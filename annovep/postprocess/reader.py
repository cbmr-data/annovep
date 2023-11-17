from __future__ import annotations

import logging
import pprint
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Union

from pydantic import BaseModel, BeforeValidator, Field, ValidationError
from typing_extensions import Annotated, TypeAlias, TypedDict

from annovep.utils import open_rb

JSON: TypeAlias = Dict[
    str, Union[str, int, float, List["JSON"], List[str], "JSON", None]
]


@dataclass
class ParsedRecord:
    vcf: VCFRecord
    vep: VEPRecord


@dataclass
class VCFRecord:
    Input: str
    Chr: str
    Pos: int
    ID: List[str]
    Ref: str
    Alts: List[str]
    Quality: Optional[float]
    Filters: List[str]
    Info: List[str]
    Samples: List[Dict[str, str]]


class VEPRecord(BaseModel):
    start: int
    strand: int
    # "most_severe_consequence": "?",
    end: int
    seq_region_name: str
    assembly_name: str
    id: str
    allele_string: str
    input: str

    transcript_consequences: List[Consequence] = Field(default_factory=list)
    intergenic_consequences: List[Consequence] = Field(default_factory=list)
    custom_annotations: Dict[str, List[CustomAnnotation]] = Field(default_factory=dict)


Consequence: TypeAlias = Dict[str, Union[str, None, List[str], int, float]]


class CustomAnnotation(BaseModel):
    # Some custom annotations use names that are interpreted as int (ClinVar)
    name: Annotated[str, BeforeValidator(lambda value: str(value))]
    allele: Optional[str] = None
    fields: Dict[str, Union[str, int, float, List[str], None]] = Field(
        default_factory=dict
    )


class MetaData(TypedDict):
    samples: List[str]


class VEPReader:
    def __init__(self, filename: Path):
        self._first_record = None

        self._log = logging.getLogger(__name__)
        self._handle = open_rb(filename)
        self.metadata = self._read_metadata()
        self.timestamp = filename.stat().st_mtime

    def _read_metadata(self) -> MetaData:
        metadata: MetaData = {"samples": []}
        for line in self._handle:
            record = self._read_record(line)

            if "AnnoVEP:Samples" in record.vcf.ID:
                metadata["samples"] = record.vcf.Info
            elif any(key.startswith("AnnoVEP:") for key in record.vcf.ID):
                self._log.warning("unexpected metadata %r", record.vcf.ID)
            else:
                self._first_record = record
                break

        return metadata

    def _read_record(
        self,
        line: bytes,
        nan_re: re.Pattern[bytes] = re.compile(rb":(-)?NaN\b", flags=re.I),
    ) -> ParsedRecord:
        # Workaround for non-standard JSON output observed in some records, where
        # an expected string value was -nan. Python accepts "NaN", but null seems
        # more reasonable for downstream compatibility
        line = nan_re.sub(b":null", line)
        try:
            vep_record = VEPRecord.model_validate_json(line, strict=True)
        except ValidationError as error:
            self._log.error("invalid record(s):\n%s", pprint.pformat(error.errors()))
            raise SystemExit(1)

        assert isinstance(vep_record.input, str)
        vcf_record = self._parse_vcf_record(vep_record.input)
        vep_record.input = vcf_record.Input

        return ParsedRecord(
            vcf=vcf_record,
            vep=vep_record,
        )

    def _parse_vcf_record(self, line: str) -> VCFRecord:
        fields = line.rstrip("\r\n").split("\t")
        chr, pos, id, ref, alt, qual, filters, info, *fmt_and_samples = fields
        chr = decode_contig_name(chr)

        # VCF record excluding any (identifying) sample information
        line = "\t".join((chr, pos, id, ref, alt, qual, filters, info))

        samples: list[dict[str, str]] = []
        if fmt_and_samples:
            fmt_keys = fmt_and_samples[0].split(":")
            for sample in fmt_and_samples[1:]:
                samples.append(dict(zip(fmt_keys, sample.split(":"))))

        return VCFRecord(
            Input=line,
            Chr=chr,
            Pos=int(pos),
            ID=[] if id == "." else id.split(";"),
            Ref=ref,
            # . is treated as a actual value, rather than an empty list. This is done so
            # that (limited) information can be retrieved for non-specific variants.
            Alts=alt.split(","),
            Quality=None if qual == "." else float(qual),
            Filters=[] if filters == "." else filters.split(";"),
            Info=[] if info == "." else info.split(";"),
            Samples=samples,
        )

    def __iter__(self) -> Iterator[ParsedRecord]:
        chrom = None
        count = 0

        if self._first_record is not None:
            yield self._first_record
            chrom = self._first_record.vcf.Chr
            count = 1

            self._first_record = None

        for line in self._handle:
            record = self._read_record(line)
            if record.vcf.Chr != chrom:
                if chrom is not None:
                    self._log.info("processed %i records on %r", count, chrom)
                chrom = record.vcf.Chr
                count = 1
            else:
                count += 1

            yield record

        if chrom is not None:
            self._log.info("processed %i records on %r", count, chrom)


def decode_contig_name(name: str) -> str:
    """Decode contig name encoded by `preprocess_vcf.py`"""
    if name.startswith("annovep_"):
        return bytes.fromhex(name[8:]).decode("utf-8")

    return name
