from __future__ import annotations

import logging
import pprint
import re
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple, Union

from pydantic import BaseModel, BeforeValidator, Field, ValidationError
from typing_extensions import Annotated, Literal, TypedDict

from annovep.utils import open_rb


class VEPRecord(TypedDict):
    Chr: str
    Pos: int
    ID: List[str]
    Ref: str
    Alts: List[str]
    Quality: Optional[float]
    Filters: List[str]
    Info: List[str]
    Samples: List[Dict[str, str]]
    VEP: VEPData


# Combined transcript / intergenic consequence
class Consequence(BaseModel):
    amino_acids: Optional[str] = None
    cdna_end: Optional[int] = None
    cdna_start: Optional[int] = None
    cds_end: Optional[int] = None
    cds_start: Optional[int] = None
    codons: Optional[str] = None
    consequence_terms: List[str] = Field(default_factory=list)
    gene_id: Optional[str] = None
    impact: Optional[str] = None
    protein_end: Optional[int] = None
    protein_start: Optional[int] = None
    strand: Optional[Literal[1, -1]] = None
    transcript_id: Optional[str] = None
    variant_allele: Optional[str] = None
    canonical: Optional[int] = None

    # Custom fields
    n_most_significant: int = Field(default=0, alias="::annovep::1")
    most_significant: Tuple[str, ...] = Field(default=(), alias="::annovep::2")
    most_significant_canonical: Optional[str] = Field(
        default=None, alias="::annovep::3"
    )
    least_significant: Optional[str] = Field(default=None, alias="::annovep::4")

    cdna_position: Optional[str] = Field(default=None, alias="::annovep::5")
    cds_position: Optional[str] = Field(default=None, alias="::annovep::6")
    protein_position: Optional[str] = Field(default=None, alias="::annovep::7")


class Custom(BaseModel):
    # Some custom annotations use names that are interpreted as int (ClinVar)
    name: Annotated[str, BeforeValidator(lambda value: str(value))]
    allele: Optional[str] = None
    fields: Dict[str, Union[str, int, float]] = Field(default_factory=dict)


class VEPData(BaseModel):
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
    custom_annotations: Dict[str, List[Custom]] = Field(default_factory=dict)


class MetaData(TypedDict):
    samples: list[str]


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
            record_id = ";".join(record["ID"])

            if record_id == "AnnoVEP:Samples":
                metadata["samples"] = record["Info"]
            elif record_id.startswith("AnnoVEP:"):
                self._log.warning("unexpected metadata %r", record_id)
            else:
                self._first_record = record
                break

        return metadata

    def _read_record(
        self,
        line: bytes,
        nan_re: re.Pattern[bytes] = re.compile(rb":(-)?NaN\b", flags=re.I),
    ) -> VEPRecord:
        # Workaround for non-standard JSON output observed in some records, where
        # an expected string value was -nan. Python accepts "NaN", but null seems
        # more reasonable for downstream compatibility
        line = nan_re.sub(b":null", line)
        try:
            data = VEPData.model_validate_json(line, strict=True)
        except ValidationError as error:
            self._log.error("invalid record(s):\n%s", pprint.pformat(error.errors()))
            raise SystemExit(1)

        vcf_record = data.input
        fields = vcf_record.rstrip("\r\n").split("\t")
        chr, pos, id, ref, alt, qual, filters, info, *fmt_and_samples = fields
        chr = decode_contig_name(chr)

        data.input = "\t".join((chr, pos, id, ref, alt, qual, filters, info))

        samples: list[dict[str, str]] = []
        if fmt_and_samples:
            fmt_keys = fmt_and_samples[0].split(":")
            for sample in fmt_and_samples[1:]:
                samples.append(dict(zip(fmt_keys, sample.split(":"))))

        return {
            "Chr": chr,
            "Pos": int(pos),
            "ID": [] if id == "." else id.split(";"),
            "Ref": ref,
            # . is treated as a actual value, rather than an empty list. This is done so
            # that (limited) information can be retrieved for non-specific variants.
            "Alts": alt.split(","),
            "Quality": None if qual == "." else float(qual),
            "Filters": [] if filters == "." else filters.split(";"),
            "Info": [] if info == "." else info.split(";"),
            "Samples": samples,
            "VEP": data,
        }

    def __iter__(self) -> Iterator[VEPRecord]:
        chrom = None
        count = 0

        if self._first_record is not None:
            yield self._first_record
            chrom = self._first_record["Chr"]
            count = 1

            self._first_record = None

        for line in self._handle:
            record = self._read_record(line)
            if record["Chr"] != chrom:
                if chrom is not None:
                    self._log.info("processed %i records on %r", count, chrom)
                chrom = record["Chr"]
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
