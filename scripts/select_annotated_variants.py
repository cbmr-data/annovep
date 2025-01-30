#!/usr/bin/env python3
import argparse
import bz2
import collections
import gzip
import logging
import os
import shlex
import sys
import textwrap
from math import isnan
from os import fspath
from pathlib import Path
from typing import NoReturn


def abort(*args: object) -> NoReturn:
    logging.error(*args)
    sys.exit(1)


def quote(pathlike: str | Path) -> str:
    return shlex.quote(fspath(pathlike))


#######################################################################################
# Duplicate code to allow script to be used stand-alone


# Duplicated from pipeline.postprocess.consequences
def consequence_ranks() -> dict[str, int]:
    # https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
    consequences = [
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        # Consequences for NMD transcripts are ignored
        # "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "feature_elongation",
        "regulatory_region_variant",
        "feature_truncation",
        "intergenic_variant",
        ".",
    ]

    ranks: dict[str, int] = collections.OrderedDict()
    for rank, consequence in enumerate(reversed(consequences)):
        ranks[consequence] = rank

    return ranks


# Duplicated from pipeline.utils
def open_rb(filename: str | Path):
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(fspath(filename), "rb")
    try:
        header = handle.read(2)
        handle.seek(0)

        if header == b"\x1f\x8b":
            try:
                import isal.igzip

                handle = isal.igzip.GzipFile(mode="rb", fileobj=handle)
            except ModuleNotFoundError:
                handle = gzip.GzipFile(mode="rb", fileobj=handle)
        elif header == b"BZ":
            handle = bz2.BZ2File(handle, "rb")

        return handle
    except:
        handle.close()
        raise


#######################################################################################


class TSVReader:
    def __init__(self, filepath):
        self.path = filepath
        self.header = None
        self._handle = None
        self._header = None

    def __enter__(self):
        self._handle = open_rb(self.path)
        self.header = self._handle.readline()
        self._header = self.header.rstrip().split(b"\t")

        return self

    def __iter__(self):
        for idx, line in enumerate(self._handle, start=2):
            yield TSVRow(idx, line, self._header)

    def __exit__(self, _type, _value, _traceback):
        if self._handle:
            self._handle.close()


class TSVRow:
    def __init__(self, lineno, line, header):
        row = line.rstrip().split(b"\t")
        if len(row) != len(header):
            abort("Wrong number of columns on line %i in %s", lineno, quote(self.path))

        self.line = line
        self.lineno = lineno
        self._data = dict(zip(header, row))

    def get_float(self, key, nan=float("NaN")):
        return self._get_typed(float, key, nan)

    def get_int(self, key):
        return self._get_typed(int, key)

    def get_str(self, key, nan=None):
        return self._get_typed(lambda it: it.decode(), key, nan)

    def _get_typed(self, func, key, missing_value=None):
        value = self._data.get(key)
        if value is None:
            abort("Column in missing in annotations file: %r", key.decode())

        if value == b".":
            return missing_value

        try:
            return func(value)
        except Exception as error:
            abort(
                "Invalid %s value (%r) on row %i in annotations: %s",
                key.decode(),
                value.decode(),
                self.lineno,
                error,
            )

    def __getitem__(self, key):
        return self._data[key]


#######################################################################################


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 88)

        super().__init__(*args, **kwargs)

    def _fill_text(self, text, width, indent):
        output = []
        for segment in text.splitlines():
            output.append(
                textwrap.fill(
                    segment,
                    width,
                    initial_indent=indent,
                    subsequent_indent=indent,
                    replace_whitespace=False,
                )
            )

        return "\n".join(output)


def parse_args(argv, consequences):
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter,
        description="Helper script for post-processing of AnnoVEP annotated variants. "
        "Takes a VCF file and an annotation file (TSV) and filters these based on the "
        "selected criteria, writing the selected variants and annotations to new "
        "files. By default likely deleterious mutations and their annotations are "
        "selected. The VCF must be tabix indexed (see below)."
        "\n\n"
        "For a description of mutation consequence levels, see\n"
        "  https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html"
        "\n\n"
        "installation:\n"
        "  This script requires python module 'pysam' and (optionally) 'isal'.\n"
        "  To install, run the following command:\n"
        "    $ python3 -m pip install pysam isal"
        "\n\n"
        "example usage:\n"
        "  $ python3 %(prog)s --genotypes my_genotypes.vcf.gz \\\n"
        "      --annotations my_annotations.tsv.gz --output my_output\n"
        "  $ ls my_output.*\n"
        "  myoutput.vcf myoutput.tsv",
    )

    parser.add_argument(
        "-v",
        "--version",
        version="%(prog)s v1.0.0",
        action="version",
        help="show version number",
    )

    group = parser.add_argument_group("Input")
    group.add_argument(
        "-g",
        "--genotypes",
        required=True,
        type=Path,
        metavar="genotypes.vcf.gz",
        help="[REQUIRED] Path to the VCF file. MUST be bgzip compressed and tabix "
        "indexed",
    )
    group.add_argument(
        "-a",
        "--annotations",
        required=True,
        type=Path,
        metavar="annotations.tsv.gz",
        help="[REQUIRED] Path to the annotations. May be compressed",
    )

    group = parser.add_argument_group("Output")
    group.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        metavar="output_name",
        help="[REQUIRED] Prefix for output filenames",
    )

    group.add_argument(
        "--output-headers",
        action="store_true",
        help="Write the full VCF headers (if available) to output",
    )

    group = parser.add_argument_group("Filters")
    cosequence_choices = list(consequences)
    group.add_argument(
        "--greatest-consequence",
        default="transcript_ablation",
        metavar="NAME",
        choices=cosequence_choices,
        help="Exclude consequences more significant (as judged by VEP) than this "
        "consequence",
    )
    group.add_argument(
        "--smallest-consequence",
        default="splice_region_variant",
        metavar="NAME",
        choices=cosequence_choices,
        help="Exclude consequences less significant (as judged by VEP) than this "
        "consequence",
    )
    group.add_argument(
        "--consequence",
        metavar="NAME",
        choices=cosequence_choices,
        help="Equivalent to using --smallest-consequence and --greatest-consequence"
        "with the same consequence, i.e. only select this consequence",
    )

    group.add_argument(
        "--max-gnomAD",
        type=float,
        default=0.01,
        metavar="X",
        help="Excluded variants where the maximum gnomAD frequency is greater than X",
    )
    group.add_argument(
        "--min-gnomAD",
        type=float,
        default=0.0,
        metavar="X",
        help="Excluded variants where the maximum gnomAD frequency is less than X",
    )

    group.add_argument(
        "--ignore-filters",
        action="store_true",
        help="Set to allow FILTERS column values other than 'PASS'",
    )

    return parser.parse_args(argv)


def main(argv):
    consequences = consequence_ranks()
    args = parse_args(argv, consequences)

    if args.consequence is not None:
        args.smallest_consequence = args.consequence
        args.greatest_consequence = args.consequence

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    log = logging.getLogger(__name__)

    try:
        # ISA-L is significantly faster than the built-in gzip decompressor
        import isal.igzip
    except ModuleNotFoundError:
        log.warning("Python module 'isal' not installed; will use slow gzip reader")
        log.warning("To install, run `%s -m pip install isal`", quote(sys.executable))

    try:
        import pysam
    except ModuleNotFoundError:
        log.error("Python module 'pysam' not installed; cannot proceed!")
        log.error("To install, run `%s -m pip install pysam`", quote(sys.executable))

        return 1

    tabix_file = f"{args.genotypes}.tbi"
    if not args.genotypes.is_file():
        log.error("Genotypes do not exist/are not a file: %s", quote(args.genotypes))
        return 1
    elif not args.annotations.is_file():
        log.error(
            "Annotations do not exist/are not a file: %s", quote(args.annotations)
        )
        return 1
    elif not os.path.isfile(tabix_file):
        log.error("Genotypes are not tabix indexed!")
        log.error("To index the file, run `tabix -p vcf %s`", quote(tabix_file))
        return 1
    elif not args.genotypes.name.endswith(".vcf.gz"):
        log.warning("Genotypes do not appear to be a bgzip compressed VCF")
        log.warning("To compress the VCF, run `bgzip %s`", quote(args.genotypes))

    min_consequence = consequences[args.smallest_consequence]
    max_consequence = consequences[args.greatest_consequence]

    n_filtered_filter = 0
    n_filtered_gnomAD = 0
    n_filtered_consequence = 0
    n_selected = 0

    log.info("Reading genotypes from %s", quote(args.genotypes))
    log.info("Reading annotations from %s", quote(args.annotations))

    out_vcf_name = f"{args.output}.vcf"
    out_tsv_name = f"{args.output}.tsv"
    log.info("Writing filtered genotypes to %s", quote(out_vcf_name))
    log.info("Writing filtered annotations to %s", quote(out_tsv_name))

    with pysam.TabixFile(fspath(args.genotypes)) as vcf_handle, TSVReader(
        args.annotations
    ) as tsv_handle, open(out_vcf_name, "w") as out_vcf, open(
        out_tsv_name, "wb"
    ) as out_tsv:
        out_tsv.write(tsv_handle.header)

        if vcf_handle.header:
            vcf_header = vcf_handle.header
            if not args.output_headers:
                vcf_header = vcf_header[-1:]

            out_vcf.write("\n".join(vcf_header))
            out_vcf.write("\n")

        last_record = None
        for row in tsv_handle:
            idx = row.lineno
            if idx % 1000000 == 0 and idx:
                log.info(
                    "Processed {:,} variants, selected {}; currently at {}:{}".format(
                        idx,
                        n_selected,
                        row.get_str(b"#Chr"),
                        row.get_int(b"Pos"),
                    )
                )

            if not (args.ignore_filters or row[b"Filters"] == b"PASS"):
                n_filtered_filter += 1
                continue

            gnomAD_max = row.get_float(b"gnomAD_max")
            if not (args.min_gnomAD <= gnomAD_max <= args.max_gnomAD):
                if not isnan(gnomAD_max):
                    n_filtered_gnomAD += 1
                    continue

            consequence_str = row.get_str(b"Func_most_significant", ".")
            consequence = consequences.get(consequence_str)
            if consequence is None:
                abort("unknown consequence on line %i: %r", row.lineno, consequence_str)

            if not (min_consequence <= consequence <= max_consequence):
                n_filtered_consequence += 1
                continue

            chrom = row.get_str(b"#Chr")
            pos = row.get_int(b"Pos")
            ref = row.get_str(b"Ref")
            alt = row.get_str(b"Alt")

            n_selected += 1
            out_tsv.write(row.line)

            for record in vcf_handle.fetch(chrom, pos - 1, pos):
                vchrom, vpos, _, vref, valt, _ = record.split("\t", 5)
                if vchrom == chrom and int(vpos) == pos and vref == ref:
                    # Avoid printing duplicate variants if multiple annotations
                    # on the same position pass the filters
                    if alt in valt.split(",") and record != last_record:
                        last_record = record
                        out_vcf.write(record)
                        out_vcf.write("\n")

    total = n_filtered_consequence + n_filtered_filter + n_filtered_gnomAD + n_selected

    log.info(f"Processed {total:,} variants")
    log.info(f"{n_filtered_filter:,} variants excluded by FILTER")
    log.info(f"{n_filtered_gnomAD:,} variants excluded by gnomAD max")
    log.info(f"{n_filtered_consequence:,} variants excluded by consequence")
    log.info(f"{n_selected:,} variants selected")

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
