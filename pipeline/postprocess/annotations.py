import collections
import functools
import re

import liftover

from . import consequences

from annotation import Builtin, Custom, Field, Option, Plugin

_RE_ALLELE = re.compile("\b")


class Annotator:
    def __init__(self, annotations, metadata=None, liftover_cache=None) -> None:
        self.groups = annotations
        self._consequence_ranks = consequences.ranks()
        self._lifter = liftover.get_lifter("hg38", "hg19", liftover_cache)

        self._apply_metadata(metadata)

        self.fields = collections.OrderedDict()
        for annotation in self.groups:
            for field in annotation.fields.values():
                if field.name is not None:
                    self.fields[field.name] = field

    def _apply_metadata(self, metadata):
        for annotation in self.groups:
            if isinstance(annotation, Builtin):
                if annotation.name.lower() == "samplegenotypes":
                    annotation.fields = {}
                    for sample in metadata["samples"]:
                        annotation.fields[sample] = Field(
                            name=sample,
                            type="str",
                            help=f"Genotypes for {sample!r}",
                            split_by=None,
                        )
                else:
                    raise NotImplementedError(
                        f"{annotation.name} not a builtin annotation"
                    )

    def annotate(self, row):
        vep = row.pop("VEP")
        samples = row["Samples"]

        row["Ref"] = self._validate_sequence(row["Ref"], "ACGTN*")
        row["DP"] = self._calculate_depth(samples)

        genotype_counts = self._count_genotypes(samples)
        frequencies = self._calculate_allele_freqs(genotype_counts)

        # Construct the cleaned up alleles / positions used by VEP
        vep_alleles = self._construct_vep_alleles(row)

        for allele_idx, allele in enumerate(row["Alts"], start=1):
            allele = self._validate_sequence(allele, "ACGTN*.")

            copy = dict(row)
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
            self._add_liftover_annotations(vep, copy)

            yield copy

    def _validate_sequence(self, sequence, whitelist):
        sequence = sequence.upper()

        # Don't bother supporting old/weird VCFs
        invalid_characters = set(sequence).difference(whitelist)
        if invalid_characters:
            raise ValueError(sequence)

        return sequence

    def _calculate_depth(self, samples):
        any_depth = False
        total_depth = 0
        for sample in samples:
            # Both missing key/value pairs and '.' values are possible
            depth = sample.get("DP", ".")
            if depth != ".":
                any_depth = True
                total_depth += int(depth)

        return total_depth if any_depth else None

    def _count_genotypes(self, samples):
        counts = collections.defaultdict(int)
        for sample in samples:
            counts[parse_vcf_genotypes(sample.get("GT"))] += 1

        return dict(counts)

    def _calculate_allele_freqs(self, counts):
        allele_counts = collections.defaultdict(int)
        for alleles, count in counts.items():
            if alleles != (None, None):
                for allele in alleles:
                    allele_counts[allele] += count

        frequencies = {}
        total = sum(allele_counts.values())
        for key, value in allele_counts.items():
            frequencies[key] = "{:.4g}".format(value / total)

        return frequencies

    def _construct_vep_alleles(self, row):
        start = row["Pos"]
        ref = row["Ref"]
        alts = row["Alts"]
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
            for vcf_alt, vep_alt in zip(row["Alts"], alts)
        }

    def _get_allele_consequence(self, vep, allele):
        # The JSON record contains transcript, integenic, or no consequences
        transcript_consequences = vep.get("transcript_consequences", ())
        intergenic_consequences = vep.get("intergenic_consequences", ())
        assert not (transcript_consequences and intergenic_consequences), vep

        consequences = []
        canonical_consequences = []
        for consequence in transcript_consequences or intergenic_consequences:
            if consequence["variant_allele"] == allele:
                consequence_terms = consequence["consequence_terms"]
                if "NMD_transcript_variant" in consequence_terms:
                    # Consequences for NMD transcripts are not informative
                    continue

                # Gene ID will be missing for intergenetic consequences
                gene = consequence.get("gene_id")
                is_canonical = consequence.get("canonical")

                for term in consequence_terms:
                    entry = (self._consequence_ranks[term], term, gene, consequence)

                    consequences.append(entry)
                    if is_canonical:
                        canonical_consequences.append(entry)

        if not consequences:
            # No consequences for non-variable sites or sites with only NMD consequences
            return {}

        consequences.sort(key=lambda it: it[0])
        canonical_consequences.sort(key=lambda it: it[0])

        # One of the most significant consequences is picked "randomly"
        _, most_significant, gene_id, consequence = consequences[0]

        if canonical_consequences:
            _, most_significant_canonical, _, _ = canonical_consequences[0]
        else:
            most_significant_canonical = None

        consequence["most_significant_canonical"] = most_significant_canonical
        consequence["most_significant"] = most_significant
        consequence["n_most_significant"] = 0
        for _, term, _, _ in consequences:
            if term != most_significant:
                break

            consequence["n_most_significant"] += 1

        for _, term, gene_id_, _ in reversed(consequences):
            if gene_id_ == gene_id:
                consequence["least_significant"] = term
                break

        # Convert start/end coordinates into single value
        for key in ("cdna", "cds", "protein"):
            consequence[f"{key}_position"] = self._format_coordinates(consequence, key)

        return consequence

    def _format_coordinates(self, consequence, key):
        start = consequence.get(f"{key}_start")
        end = consequence.get(f"{key}_end")

        if start is None:
            if end is None:
                return None

            return f"?-{end}"
        elif end is None:
            return f"{start}-?"
        elif start == end:
            return start

        return f"{start}-{end}"

    def _fix_ancestral_allele(self, consequence, dst):
        # The explicty check for falsey values is used to catch both missing values and
        # as a workaround for bug where "aa" is -nan (converted to None in _read_record)
        if "Ancestral_allele" in dst and not dst["Ancestral_allele"]:
            dst["Ancestral_allele"] = None

    def _add_option_and_plugin_annotation(self, consequence, copy):
        for annotation in self.groups:
            if isinstance(annotation, (Option, Plugin)):
                for key, field in annotation.fields.items():
                    if field.name is not None:
                        copy.setdefault(field.name, consequence.get(key))

    def _add_custom_annotation(self, src, dst):
        for annotation in self.groups:
            if isinstance(annotation, Custom):
                data = {}

                alt = dst[":vep:"]["alt"]
                # annotation = {"fields": ..., "name": "chr:start-end"}
                results = src.get("custom_annotations", {}).get(annotation.name, [])
                for result in results:
                    if result.get("allele", alt) == alt:
                        # data = {"FILTER": ".", field1: value1, ...}
                        data = result.get("fields", {})
                        break

                numeric_values = []
                derived_values = []

                for key, field in annotation.fields.items():
                    if key in (":min:", ":max:"):
                        derived_values.append((key, field))
                    elif field.name is not None:
                        value = data.get(key)

                        if field.split_by is not None:
                            value = [] if value is None else value.split(field.split_by)

                        dst[field.name] = value
                        if value is not None and field.type in ("int", "float"):
                            numeric_values.append(value)

                for key, field in derived_values:
                    if key == ":min:":
                        value = min(numeric_values, default=None)
                    elif key == ":max:":
                        value = max(numeric_values, default=None)
                    else:
                        raise NotImplementedError(key)

                    dst[field.name] = value

    def _add_builtin_annotation(self, src, dst):
        for annotation in self.groups:
            if isinstance(annotation, Builtin):
                if annotation.name == "SampleGenotypes":
                    self._add_sample_genotypes(annotation, src, dst)
                else:
                    raise NotImplementedError(annotation.name)

    def _add_sample_genotypes(self, annotation, src, dst):
        alt_idx = dst["Alts"].index(dst["Alt"]) + 1

        for name, data in zip(annotation.fields, dst["Samples"]):
            genotypes = data.get("GT", "./.")
            if genotypes != "./.":
                values = []
                for field in _RE_ALLELE.split(genotypes):
                    if field == "0":
                        values.append("0")
                    elif field == alt_idx:
                        values.append("1")
                    elif field.isdigit():
                        values.append("x")
                    elif field:
                        values.append(field)

                dst[name] = "".join(values)
            else:
                dst[name] = genotypes

    def _add_neighbouring_genes(self, src, dst, nnearest=3):
        # Check if pipeline was run with the 'overlapping' BED file
        if "Genes_overlapping" not in dst:
            return

        # Start coordinate of VEP allele
        astart = src["start"]
        # End coordinate of the allele. This is shared between all ALTs
        aend = src["end"]

        # Alleles may cover multiple bases, so we may have multiple instances of each
        # category, some of which may also be overlapping with the allele
        neighbours = {"downstream": set(), "upstream": set(), "overlap": set()}

        # annotation = {"fields": ..., "name": "chr:start-end"}
        annotations = src.get("custom_annotations", {}).get("neighbours", [])
        for annotation in annotations:
            values = annotation["name"].split(";")

            # The first value (the category) is skipped
            for gene in values[1:]:
                nstart_end, name = gene.split(":")
                nstart, nend = nstart_end.split("-")
                nstart = int(nstart)
                nend = int(nend)

                if aend < nstart:
                    neighbours["downstream"].add((nstart - aend, name))
                elif astart > nend:
                    neighbours["upstream"].add((astart - nend, name))
                else:
                    neighbours["overlap"].add(name)

        def _to_list(values):
            values = sorted(values)[:nnearest]
            if not values:
                return []

            return [f"{distance}:{name}" for distance, name in values]

        dst["Genes_overlapping"] = sorted(neighbours["overlap"])
        dst["Genes_upstream"] = _to_list(neighbours["upstream"])
        dst["Genes_downstream"] = _to_list(neighbours["downstream"])

    def _add_liftover_annotations(self, vep, row):
        src_chr = row["Chr"]

        # Returns list of overlapping liftover coordinates, an empty list if the
        # position does not exist in the target genome, or KeyError if unknown.
        try:
            coordinates = self._lifter.query(src_chr, row["Pos"])
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
def parse_vcf_genotypes(genotypes, _re=re.compile(r"[|/]")):
    if genotypes in (None, "./.", ".|."):
        return (None, None)

    result = tuple(int(value) for value in _re.split(genotypes))
    if len(result) != 2:
        raise ValueError(genotypes)

    return result
