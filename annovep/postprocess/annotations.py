from __future__ import annotations

import collections
import functools
import re
from typing import Any, Dict, Generator, Union

import liftover
from typing_extensions import TypeAlias, TypedDict

from annovep.annotation import Annotation, AnnotationField
from annovep.postprocess import consequences
from annovep.postprocess.reader import Consequence, MetaData, VEPData, VEPRecord

_RE_ALLELE = re.compile(r"[/|]")


class VEPAllele(TypedDict):
    start: int
    ref: str
    alt: str
    alleles: str


OutputDict: TypeAlias = Dict[str, Union[int, str, float, None]]


class Annotator:
    def __init__(
        self,
        annotations: list[Annotation],
        metadata: MetaData | None = None,
        liftover_cache: str | None = None,
    ) -> None:
        self.groups = annotations
        self._consequence_ranks = consequences.ranks()
        self._lifter = liftover.get_lifter("hg38", "hg19", liftover_cache)

        if metadata is not None:
            self._apply_metadata(metadata)

        self.fields: list[AnnotationField] = []
        for annotation in self.groups:
            self.fields.extend(annotation.fields)

    def _apply_metadata(self, metadata: MetaData) -> None:
        for annotation in self.groups:
            if annotation.type == "builtin":
                if annotation.name.lower() == "samplegenotypes":
                    annotation.fields = [
                        AnnotationField(
                            input_key=sample,
                            output_key=f"GTS_{sample}",
                            type="str",
                            help=f"Genotypes for {sample!r}",
                        )
                        for sample in metadata["samples"]
                    ]
                else:
                    raise NotImplementedError(
                        f"{annotation.name} not a builtin annotation"
                    )

    def annotate(self, record: VEPRecord) -> Generator[dict[str, Any], None, None]:
        vep = record["VEP"]
        samples = record["Samples"]

        output: dict[str, Any] = dict(record)
        output.pop("VEP")

        output["Ref"] = self._validate_sequence(record["Ref"], "ACGTN*")
        output["DP"] = self._calculate_depth(samples)

        genotype_counts = self._count_genotypes(samples)
        frequencies = self._calculate_allele_freqs(genotype_counts)

        # Construct the cleaned up alleles / positions used by VEP
        vep_alleles = self._construct_vep_alleles(record)

        for allele_idx, allele in enumerate(record["Alts"], start=1):
            allele = self._validate_sequence(allele, "ACGTN*.")

            copy = dict(output)
            copy["Alt"] = allele
            copy["Freq"] = frequencies.get(allele_idx)

            gt_00 = genotype_counts.get((0, 0), 0)
            gt_01 = genotype_counts.get((0, allele_idx), 0)
            gt_10 = genotype_counts.get((allele_idx, 0), 0)
            gt_11 = genotype_counts.get((allele_idx, allele_idx), 0)
            gt_na = genotype_counts.get((None, None), 0)

            copy["GT_00"] = gt_00
            copy["GT_01"] = gt_01 + gt_10
            copy["GT_11"] = gt_11
            copy["GT_NA"] = gt_na
            copy["GT_other"] = (
                sum(genotype_counts.values(), 0) - gt_00 - gt_01 - gt_10 - gt_11 - gt_na
            )

            # Cleaned up coordinates/sequences used by VEP
            vep_allele = vep_alleles[allele]
            copy[":vep:"] = vep_allele

            # The position and sequences that VEP reports for this allele
            copy["VEP_allele"] = "{start}:{ref}:{alt}".format(**vep_allele)

            # Add functional annotation
            consequence = self._get_allele_consequence(vep, vep_allele["alt"])

            # add custom annotation
            self._add_option_and_plugin_annotation(consequence, copy)
            self._add_custom_annotation(vep, copy)
            self._add_builtin_annotation(vep, copy)

            # Special handling of certain (optinal) annotation
            self._fix_ancestral_allele(consequence, copy)
            self._add_neighbouring_genes(vep, copy)

            # Other computed annotations
            self._add_liftover_annotations(record, copy)

            yield copy

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

    def _construct_vep_alleles(self, record: VEPRecord) -> dict[str, VEPAllele]:
        start = record["Pos"]
        ref = record["Ref"]
        alts = record["Alts"]
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
            for vcf_alt, vep_alt in zip(record["Alts"], alts)
        }

    def _get_allele_consequence(self, vep: VEPData, allele: str) -> Consequence:
        # The JSON record contains transcript, integenic, or no consequences
        transcript_consequences = vep.transcript_consequences
        intergenic_consequences = vep.intergenic_consequences
        assert not (transcript_consequences and intergenic_consequences), vep

        consequences: list[tuple[int, str, str | None, Consequence]] = []
        canonical_consequences: list[tuple[int, str, str | None, Consequence]] = []
        for consequence in transcript_consequences or intergenic_consequences:
            if consequence.variant_allele == allele:
                consequence_terms = consequence.consequence_terms
                if "NMD_transcript_variant" in consequence_terms:
                    # Consequences for NMD transcripts are not informative
                    continue

                # Gene ID will be missing for intergenetic consequences
                gene = consequence.gene_id
                is_canonical = consequence.canonical

                for term in consequence_terms:
                    entry = (self._consequence_ranks[term], term, gene, consequence)

                    consequences.append(entry)
                    if is_canonical:
                        canonical_consequences.append(entry)

        if not consequences:
            # No consequences for non-variable sites or sites with only NMD consequences
            return Consequence()

        consequences.sort(key=lambda it: it[0])
        canonical_consequences.sort(key=lambda it: it[0])

        # One of the most significant consequences is picked "randomly"
        _, most_significant, gene_id, consequence = consequences[0]

        if canonical_consequences:
            _, most_significant_canonical, _, _ = canonical_consequences[0]
        else:
            most_significant_canonical = None

        n_most_significant = 0
        for _, term, _, _ in consequences:
            if term != most_significant:
                break

            n_most_significant += 1

        consequence.most_significant_canonical = most_significant_canonical
        consequence.most_significant = (most_significant,)
        consequence.n_most_significant = n_most_significant

        for _, term, gene_id_, _ in reversed(consequences):
            if gene_id_ == gene_id:
                consequence.least_significant = term
                break

        # Convert start/end coordinates into single value
        consequence.cdna_position = self._format_coordinates(
            consequence.cdna_start, consequence.cdna_end
        )

        consequence.cdna_position = self._format_coordinates(
            consequence.cds_start, consequence.cds_end
        )

        consequence.cdna_position = self._format_coordinates(
            consequence.protein_start, consequence.protein_end
        )

        return consequence

    def _format_coordinates(self, start: int | None, end: int | None) -> str | None:
        if start is None:
            if end is None:
                return None

            return f"?-{end}"
        elif end is None:
            return f"{start}-?"
        elif start == end:
            return str(start)

        return f"{start}-{end}"

    def _fix_ancestral_allele(
        self,
        consequence: Consequence,
        dst: dict[str, object],
    ) -> None:
        # The explicty check for falsey values is used to catch both missing values and
        # as a workaround for bug where "aa" is -nan (converted to None in _read_record)
        if "Ancestral_allele" in dst and not dst["Ancestral_allele"]:
            dst["Ancestral_allele"] = None

    def _add_option_and_plugin_annotation(
        self,
        consequence: Consequence,
        copy: dict[str, Any],
    ) -> None:
        for annotation in self.groups:
            if annotation.type == "basic" or annotation.type == "plugin":
                for field in annotation.fields:
                    if field.output_key not in copy:
                        # Consequence annoation should be explicitly marked
                        copy[field.output_key] = getattr(
                            consequence, field.input_key, None
                        )

    def _add_custom_annotation(self, src: VEPData, dst: dict[str, Any]) -> None:
        for annotation in self.groups:
            if annotation.type == "custom":
                data: dict[str, str | int | float] = {}

                alt = dst[":vep:"]["alt"]
                # annotation = {"fields": ..., "name": "chr:start-end"}
                results = src.custom_annotations.get(annotation.name, [])
                for result in results:
                    if result.allele == alt:
                        # data = {"FILTER": ".", field1: value1, ...}
                        data = result.fields
                        break

                numeric_values: list[int | float] = []
                derived_values: list[tuple[str, AnnotationField]] = []

                for field in annotation.fields:
                    if field.input_key in (":min:", ":max:"):
                        derived_values.append((field.input_key, field))
                    else:
                        value = data.get(field.input_key)

                        if field.split_by is not None:
                            assert value is None or isinstance(value, str)
                            value = [] if value is None else value.split(field.split_by)
                        elif value is not None and field.type in ("int", "float"):
                            assert isinstance(value, (int, float))
                            numeric_values.append(value)

                        dst[field.output_key] = value

                for key, field in derived_values:
                    if key == ":min:":
                        value = min(numeric_values, default=None)
                    elif key == ":max:":
                        value = max(numeric_values, default=None)
                    else:
                        raise NotImplementedError(key)

                    dst[field.output_key] = value

    def _add_builtin_annotation(self, src: VEPData, dst: dict[str, Any]) -> None:
        for annotation in self.groups:
            if annotation.type == "builtin":
                if annotation.name.lower() == "samplegenotypes":
                    self._add_sample_genotypes(annotation, src, dst)
                else:
                    raise NotImplementedError(annotation.name)

    def _add_sample_genotypes(
        self,
        annotation: Annotation,
        src: VEPData,
        dst: dict[str, Any],
    ) -> None:
        alt_genotype = str(dst["Alts"].index(dst["Alt"]) + 1)

        for name, data in zip(annotation.fields, dst["Samples"]):
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

            dst[f"GTS_{name}"] = genotypes

    def _add_neighbouring_genes(
        self,
        src: VEPData,
        dst: dict[str, Any],
        nnearest: int = 3,
    ) -> None:
        # Check if pipeline was run with the 'overlapping' BED file
        if "Genes_overlapping" not in dst:
            return

        # Start coordinate of VEP allele
        astart = src.start
        # End coordinate of the allele. This is shared between all ALTs
        aend = src.end

        # Alleles may cover multiple bases, so we may have multiple instances of each
        # category, some of which may also be overlapping with the allele
        neighbours_downstream: set[tuple[int, str]] = set()
        neighbours_upstream: set[tuple[int, str]] = set()
        neighbours_overlap: set[str] = set()

        # annotation = {"fields": ..., "name": "chr:start-end"}
        annotations = src.custom_annotations.get("neighbours", [])
        for annotation in annotations:
            values = annotation.name.split(";")

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

        def _to_list(values: set[tuple[int, str]]) -> list[object]:
            result = sorted(values)[:nnearest]
            if not result:
                return []

            return [f"{distance}:{name}" for distance, name in result]

        dst["Genes_overlapping"] = sorted(neighbours_overlap)
        dst["Genes_upstream"] = _to_list(neighbours_upstream)
        dst["Genes_downstream"] = _to_list(neighbours_downstream)

    def _add_liftover_annotations(self, vep: VEPRecord, row: dict[str, Any]) -> None:
        src_chr = vep["Chr"]

        # Returns list of overlapping liftover coordinates, an empty list if the
        # position does not exist in the target genome, or KeyError if unknown.
        try:
            coordinates = self._lifter.query(src_chr, vep["Pos"])
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

        row["Hg19_chr"] = chrom
        row["Hg19_pos"] = pos


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
