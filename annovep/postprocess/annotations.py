from __future__ import annotations

import collections
import functools
import re
from typing import TYPE_CHECKING, Generator

import liftover
from typing_extensions import TypedDict

from annovep.annotation import Annotation, AnnotationField
from annovep.postprocess import consequences
from annovep.postprocess.reader import (
    JSON,
    Consequence,
    CustomAnnotation,
    MetaData,
    ParsedRecord,
    VCFRecord,
    VEPRecord,
)

if TYPE_CHECKING:
    import logging

_RE_ALLELE = re.compile(r"[/|]")


class VEPAllele(TypedDict):
    start: int
    ref: str
    alt: str
    alleles: str


class Annotator:
    __slots__ = [
        "_builtin_liftover",
        "_builtin_sample_genotypes",
        "_consequence_ranks",
        "_counter",
        "_fields_collected",
        "_fields_derived",
        "_liftover_cache",
        "_liftover",
        "fields",
        "groups",
    ]

    _builtin_liftover: bool
    _builtin_sample_genotypes: list[AnnotationField]
    _consequence_ranks: dict[str, int]
    _counter: dict[str, int]
    _fields_collected: tuple[AnnotationField, ...]
    _fields_derived: tuple[AnnotationField, ...]
    _liftover_cache: str | None
    _liftover: liftover.ChainFile | None
    fields: tuple[AnnotationField, ...]
    groups: tuple[Annotation, ...]

    def __init__(
        self,
        annotations: list[Annotation],
        metadata: MetaData | None = None,
        liftover_cache: str | None = None,
    ) -> None:
        self.groups = tuple(annotations)
        self._consequence_ranks = consequences.ranks()

        self._builtin_sample_genotypes = []
        self._builtin_liftover = False
        self._liftover_cache = liftover_cache
        self._liftover = None

        if metadata is not None:
            self._apply_metadata(metadata)

        fields: list[AnnotationField] = []
        for annotation in self.groups:
            fields.extend(annotation.fields)

            if annotation.type == "builtin":
                if annotation.name == "SampleGenotypes":
                    self._builtin_sample_genotypes = annotation.fields
                elif annotation.name == "Liftover":
                    self._builtin_liftover = True

        self.fields = tuple(fields)
        self._fields_derived = tuple(it for it in self.fields if it.derived_from)
        self._fields_collected = tuple(it for it in self.fields if not it.derived_from)
        self._counter = {field.output_key: 0 for field in self.fields}

    def _apply_metadata(self, metadata: MetaData) -> None:
        # FIXME: Move
        for annotation in self.groups:
            if annotation.type == "builtin":
                if annotation.name == "SampleGenotypes":
                    annotation.fields = [
                        AnnotationField(
                            input_path=["genotypes", f"GTS_{sample}"],
                            output_key=f"GTS_{sample}",
                            type="str",
                            derived_from=[],
                            help=f"Genotypes for {sample!r}",
                        )
                        for sample in metadata["samples"]
                    ]
                elif annotation.name == "Liftover":
                    pass
                else:
                    raise NotImplementedError(
                        f"{annotation.name} not a builtin annotation"
                    )

    def annotate(self, record: ParsedRecord) -> Generator[JSON, None, None]:
        for data in self._collect_annotations(record):
            output: JSON = {}
            for field in self._fields_collected:
                source = data
                for key in field.input_path:
                    if source is not None:
                        assert isinstance(source, dict)
                        source = source.get(key)
                    else:
                        break

                if field.split_by is not None:
                    assert source is None or isinstance(source, str)
                    source = [] if source is None else source.split(field.split_by)

                if source is not None:
                    self._counter[field.output_key] += 1

                output[field.output_key] = source

            for field in self._fields_derived:
                values: list[int | float] = []
                for key in field.derived_from:
                    value = output[key]
                    if value is not None:
                        assert isinstance(value, (float, int))
                        values.append(value)

                value = None
                if field.input_path[-1] == ":min:":
                    value = min(values, default=None)
                elif field.input_path[-1] == ":max:":
                    value = max(values, default=None)

                if value is not None:
                    self._counter[field.output_key] += 1

                output[field.output_key] = value

            yield output

    def finalize(self, log: logging.Logger) -> None:
        log.info("Finalizing annotations")
        for field, count in self._counter.items():
            func = log.debug if count else log.warning
            func("  %s values collected for column %s", count, field)

    def _collect_annotations(self, record: ParsedRecord) -> Generator[JSON, None, None]:
        vep = record.vep
        vcf = record.vcf
        # Ensure that alleles are within expectations
        vcf.Alts = [self._validate_sequence(allele, "ACGTN*.") for allele in vcf.Alts]
        vcf.Ref = self._validate_sequence(vcf.Ref, "ACGTN*")

        site_depth = self._calculate_depth(vcf.Samples)
        genotype_counts = self._count_genotypes(vcf.Samples)
        frequencies = self._calculate_allele_freqs(genotype_counts)

        # Construct the cleaned up alleles / positions used by VEP
        vep_alleles = self._construct_vep_alleles(vcf)

        # Special handling of certain (optinal) annotation
        self._parse_neighbouring_genes(vep)

        liftover: JSON = {}
        if self._builtin_liftover:
            liftover = self._add_liftover_annotations(vcf)

        for allele_idx, allele in enumerate(vcf.Alts, start=1):
            # Cleaned up coordinates/sequences used by VEP
            vep_allele = vep_alleles[allele]

            # Genotype counts
            gt_00 = genotype_counts.get((0, 0), 0)
            gt_01 = genotype_counts.get((0, allele_idx), 0)
            gt_10 = genotype_counts.get((allele_idx, 0), 0)
            gt_11 = genotype_counts.get((allele_idx, allele_idx), 0)
            gt_na = genotype_counts.get((None, None), 0)

            sample_genotypes = {}
            if self._builtin_sample_genotypes:
                sample_genotypes = self._get_sample_genotypes(vcf, allele)

            yield {
                "input": vars(vcf),
                "consequence": self._get_allele_consequence(
                    vep=vep,
                    allele=vep_allele["alt"],
                ),
                "canonical_consequence": self._get_allele_consequence(
                    vep=vep,
                    allele=vep_allele["alt"],
                    canonical_only=True,
                ),
                "custom": self._get_custom_annotation(
                    vep=vep,
                    allele=vep_allele["alt"],
                ),
                "derived": {
                    # Total read depth at the site (all alleles)
                    "DP": site_depth,
                    # Current allele
                    "Alt": allele,
                    # Frequency of the current allele
                    "Freq": frequencies.get(allele_idx),
                    # The position and sequences that VEP reports for this allele
                    "VEP_allele": "{start}:{ref}:{alt}".format_map(vep_allele),
                },
                "genotypes": {
                    **sample_genotypes,
                    "GT_00": gt_00,
                    "GT_01": gt_01 + gt_10,
                    "GT_11": gt_11,
                    "GT_NA": gt_na,
                    "GT_other": (
                        sum(genotype_counts.values(), 0)
                        - gt_00
                        - gt_01
                        - gt_10
                        - gt_11
                        - gt_na
                    ),
                },
                "liftover": liftover,
            }

    def _validate_sequence(self, sequence: str, whitelist: str) -> str:
        sequence = sequence.upper()

        # Don't bother supporting old/weird VCFs
        invalid_characters = set(sequence).difference(whitelist)
        if invalid_characters:
            raise ValueError(sequence)

        return sequence

    def _calculate_depth(self, samples: list[dict[str, str]]) -> int | None:
        any_depth = False
        total_depth = 0
        for sample in samples:
            # Both missing key/value pairs and '.' values are possible
            depth = sample.get("DP", ".")
            if depth != ".":
                any_depth = True
                total_depth += int(depth)

        return total_depth if any_depth else None

    def _count_genotypes(
        self, samples: list[dict[str, str]]
    ) -> dict[tuple[None, None] | tuple[int, int], int]:
        counts: dict[
            tuple[None, None] | tuple[int, int], int
        ] = collections.defaultdict(int)
        for sample in samples:
            counts[parse_vcf_genotypes(sample.get("GT"))] += 1

        return dict(counts)

    def _calculate_allele_freqs(
        self,
        counts: dict[tuple[None, None] | tuple[int, int], int],
    ) -> dict[int, str]:
        allele_counts: dict[int, int] = collections.defaultdict(int)
        for alleles, count in counts.items():
            for allele in alleles:
                if allele is not None:
                    allele_counts[allele] += count

        frequencies: dict[int, str] = {}
        total = sum(allele_counts.values())
        for key, value in allele_counts.items():
            frequencies[key] = "{:.4g}".format(value / total)

        return frequencies

    def _construct_vep_alleles(self, record: VCFRecord) -> dict[str, VEPAllele]:
        start = record.Pos
        ref = record.Ref
        alts = record.Alts
        if alts == ["."]:
            # VEP has limited support for including records for non-variants
            return {".": {"start": start, "ref": ref, "alt": ref, "alleles": ref}}

        if any(len(ref) != len(allele) for allele in alts):
            # find out if all the alts start with the same base, ignoring "*"
            any_non_star = False
            first_bases = {ref[0]}
            for alt in alts:
                if not alt.startswith("*"):
                    any_non_star = True
                    first_bases.add(alt[0])

            if any_non_star and len(first_bases) == 1:
                start += 1
                ref = ref[1:] or "-"
                alts = list(alts)
                for idx, alt in enumerate(alts):
                    if alt.startswith("*"):
                        alts[idx] = alt
                    else:
                        alts[idx] = alt[1:] or "-"

        # String generated by VEP summarizing ref/alleles in a JSON record
        allele_string = "{}/{}".format(ref, "/".join(alts))

        return {
            vcf_alt: {
                "start": start,
                "ref": ref,
                "alt": vep_alt,
                "alleles": allele_string,
            }
            for vcf_alt, vep_alt in zip(record.Alts, alts)
        }

    def _get_allele_consequence(
        self,
        vep: VEPRecord,
        allele: str,
        canonical_only: bool = False,
    ) -> JSON:
        # The JSON record contains transcript, integenic, or no consequences
        transcript_consequences = vep.transcript_consequences
        intergenic_consequences = vep.intergenic_consequences
        assert not (transcript_consequences and intergenic_consequences), vep

        consequences: list[tuple[int, str, str | None, Consequence]] = []
        for consequence in transcript_consequences or intergenic_consequences:
            if consequence["variant_allele"] == allele:
                consequence_terms = consequence["consequence_terms"]
                assert isinstance(consequence_terms, list)
                if "NMD_transcript_variant" in consequence_terms:
                    # Consequences for NMD transcripts are not informative
                    continue

                # Gene ID will be missing for intergenetic consequences
                gene = consequence.get("gene_id")
                assert gene is None or isinstance(gene, str)

                if consequence.get("canonical") or not canonical_only:
                    for term in consequence_terms:
                        assert isinstance(term, str)
                        entry = (self._consequence_ranks[term], term, gene, consequence)
                        consequences.append(entry)

        if not consequences:
            # No consequences for non-variable sites or sites with only NMD consequences
            return {}

        consequences.sort(key=lambda it: it[0])
        # One of the most significant consequences is picked "randomly"
        _, most_significant, gene_id, consequence = consequences[0]

        n_most_significant = 0
        for _, term, _, _ in consequences:
            if term != most_significant:
                break

            n_most_significant += 1

        output: JSON = {
            **consequence,
            ":most_significant:": most_significant,
            ":most_significant_count:": n_most_significant,
            # Convert start/end coordinates into single value
            ":cdna_position:": self._format_coordinates(consequence, "cdna"),
            ":cds_position:": self._format_coordinates(consequence, "cds"),
            ":protein_position:": self._format_coordinates(consequence, "protein"),
        }

        for _, term, gene_id_, _ in reversed(consequences):
            if gene_id_ == gene_id:
                output[":least_significant:"] = term
                break

        # The explicty check for falsey values is used to catch both missing values and
        # as a workaround for bug where "aa" is -nan (converted to None in _read_record)
        if "aa" in output and not output["aa"]:
            output["aa"] = None

        return output

    def _format_coordinates(self, data: Consequence, key: str) -> str | None:
        start = data.get(f"{key}_start")
        assert start is None or isinstance(start, int)
        end = data.get(f"{key}_end")
        assert end is None or isinstance(end, int)

        if start is None:
            if end is None:
                return None

            return f"?-{end}"
        elif end is None:
            return f"{start}-?"
        elif start == end:
            return str(start)

        return f"{start}-{end}"

    ####################################################################################
    # Custom annotation

    def _get_custom_annotation(self, vep: VEPRecord, allele: str) -> JSON:
        output: JSON = {}
        for key, records in vep.custom_annotations.items():
            for record in records:
                if record.allele is None or record.allele == allele:
                    fields = record.fields

                    output[key] = {
                        "name": record.name,
                        "allele": record.allele,
                        **fields,
                    }

                    # FIXME: BED annotation requires handling multiple records (gnomAD)
                    break

        return output

    def _parse_neighbouring_genes(self, vep: VEPRecord, nnearest: int = 3) -> None:
        custom_annotations = vep.custom_annotations
        assert isinstance(custom_annotations, dict)

        neighbours = custom_annotations.get("neighbours")
        if neighbours is None:
            return

        # Start coordinate of VEP allele
        astart = vep.start
        # End coordinate of the allele. This is shared between all ALTs
        aend = vep.end

        # Alleles may cover multiple bases, so we may have multiple instances of each
        # category, some of which may also be overlapping with the allele
        neighbours_downstream: set[tuple[int, str]] = set()
        neighbours_upstream: set[tuple[int, str]] = set()
        neighbours_overlap: set[str] = set()

        def _neighbours_to_list(values: set[tuple[int, str]]) -> list[str]:
            return [
                f"{distance}:{name}" for distance, name in sorted(values)[:nnearest]
            ]

        for record in neighbours:
            values = record.name.split(";")

            # The first value (the category) is skipped
            for gene in values[1:]:
                nstart_end, name = gene.split(":")
                nstart, nend = nstart_end.split("-")
                nstart = int(nstart)
                nend = int(nend)

                if aend < nstart:
                    neighbours_downstream.add((nstart - aend, name))
                elif astart > nend:
                    neighbours_upstream.add((astart - nend, name))
                else:
                    neighbours_overlap.add(name)

        neighbours[:] = [
            CustomAnnotation(
                name="neighbours",
                fields={
                    "overlapping": sorted(neighbours_overlap),
                    "upstream": _neighbours_to_list(neighbours_upstream),
                    "downstream": _neighbours_to_list(neighbours_downstream),
                },
            )
        ]

    def _get_sample_genotypes(self, vcf: VCFRecord, allele: str) -> JSON:
        alt_genotype = str(vcf.Alts.index(allele) + 1)

        result: JSON = {}
        for field, data in zip(self._builtin_sample_genotypes, vcf.Samples):
            genotypes = data.get("GT", "./.")
            if genotypes != "./.":
                values: list[str] = []
                for value in _RE_ALLELE.split(genotypes):
                    if value == "0":
                        values.append(value)
                    elif value == alt_genotype:
                        values.append("1")
                    else:
                        values.append("x")

                # Normalize x/1 to 1/x
                if values[0] == "x":
                    values = values[::-1]

                genotypes = "/".join(values)

            result[field.output_key] = genotypes

        return result

    ####################################################################################
    # Built-in annotations

    def _add_liftover_annotations(self, vcf: VCFRecord) -> JSON:
        if self._liftover is None:
            self._liftover = liftover.get_lifter("hg38", "hg19", self._liftover_cache)

        src_chr = vcf.Chr

        # Returns list of overlapping liftover coordinates, an empty list if the
        # position does not exist in the target genome, or KeyError if unknown.
        try:
            coordinates = self._liftover.query(src_chr, vcf.Pos)
        except KeyError:
            coordinates = None

        chrom = pos = None
        if coordinates:
            # It's unclear if multiple coordinates can be returned so just use the first
            chrom, pos, _ = coordinates[0]

            if src_chr.startswith("chr") and not chrom.startswith("chr"):
                chrom = "chr" + chrom
            elif chrom.startswith("chr") and not src_chr.startswith("chr"):
                chrom = chrom[3:]

        return {
            "Hg19_chr": chrom,
            "Hg19_pos": pos,
        }


@functools.lru_cache()
def parse_vcf_genotypes(
    genotypes: str,
    _re: re.Pattern[str] = re.compile(r"[|/]"),
) -> tuple[None, None] | tuple[int, int]:
    if genotypes in (None, "./.", ".|."):
        return (None, None)

    result = tuple(int(value) for value in _re.split(genotypes))
    if len(result) != 2:
        raise ValueError(genotypes)

    return result
