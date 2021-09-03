#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import bz2
import collections
import gzip
import io
import itertools
import logging
import sys
from os import fspath
from pathlib import Path
from typing import IO, Union, cast

import coloredlogs
import pysam

TEMPLATE = "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}"


class Counter:
    STEP = 1e6

    def __init__(self, log):
        self._log = log
        self._chrom = None
        self._count = 0
        self._next_print = Counter.STEP
        self._skipped = collections.defaultdict(int)

    def __call__(self, chrom, num=1):
        if chrom != self._chrom:
            self.print()
            self._chrom = chrom
            self._count = num
            self._next_print = Counter.STEP
            return

        self._count += num
        if self._count >= self._next_print:
            self._next_print += Counter.STEP
            self.print()

    def skip(self, chrom, num=1):
        self(chrom, num)
        self._skipped[chrom] += num

    def print(self):
        if self._chrom is not None:
            count = "{:,}".format(self._count)
            self._log.info("Processed %s sites on %s", count, self._chrom)

    def __enter__(self, *args, **kwargs):
        return self

    def __exit__(self, *args, **kwargs):
        self.print()

        for name, count in self._skipped.items():
            count = "{:,}".format(count)
            self._log.info("  - Skipped %s records on %r", count, name)


# Reads file like GCF_000001405.39_GRCh38.p13_assembly_report.txt
def read_assembly_report(filepath):
    mapping = {}
    with open_ro(filepath) as handle:
        # The header is the last comment
        for header, line in zip(handle, handle):
            if not line.startswith("#"):
                break
        else:
            assert False, filepath

        header = header.lstrip("#").strip().split("\t")
        header = [name.replace("-", "_").lower() for name in header]

        while line:
            row = dict(zip(header, line.rstrip().split("\t")))

            sequence_role = row["sequence_role"]
            if sequence_role == "assembled-molecule":
                template = "chr{assigned_molecule}"
            elif sequence_role == "alt-scaffold":
                template = "chr{assigned_molecule}_{genbank_accn}_alt"
            elif sequence_role == "unlocalized-scaffold":
                template = "chr{assigned_molecule}_{genbank_accn}_random"
            elif sequence_role == "unplaced-scaffold":
                template = "chrUn_{genbank_accn}"
            elif sequence_role in ("fix-patch", "novel-patch"):
                line = handle.readline()
                continue
            else:
                raise ValueError(sequence_role)

            if "na" not in (row["refseq_accn"], row["genbank_accn"]):
                row["genbank_accn"] = row["genbank_accn"].replace(".", "v")

                mapping[row["refseq_accn"]] = template.format(**row)

            line = handle.readline()

    return mapping


def refseq_names_to_chr(log, counter, args):
    log.info("Creating refseq<->chr contig name mapping from %r", args.assembly_report)

    print("RefSeq\tName")
    for refseq, name in read_assembly_report(args.assembly_report).items():
        print("{}\t{}".format(refseq, name))


# https://www.ncbi.nlm.nih.gov/snp/docs/products/vcf/redesign/#tags_with_SOterm
DBSNP_FUNCTIONS = {
    "ASS": "splice_acceptor_variant",
    "DSS": "splice_donor_variant",
    "INT": "intron_variant",
    "NSF": "frameshift_variant",
    "NSM": "missense_variant",
    "NSN": "stop_gained",
    "R3": "500B_downstream_variant",
    "R5": "2KB_upstream_variant",
    "SYN": "synonymous_variant",
    "U3": "3_prime_UTR_variant",
    "U5": "5_prime_UTR_variant",
}

DBSNP_HEADER = """\
##fileformat=VCFv4.2
##INFO=<ID=ids,Number=0,Type=String,Description="List of rsIDs for this SNP.">
##INFO=<ID=alts,Number=0,Type=String,Description="List of sets of alleles observed for a given chr:pos:ref (e.g. "A,A/T,A/T/G").">
##INFO=<ID=functions,Number=0,Type=String,Description="List of GO terms assosiated with this SNP in DBSNP.">
#CHROM  POS   ID REF  ALT   QUAL   FILTER INFO
"""


def dbsnp_to_vcf(log, counter, args):
    log.info("Creating custom DBSNP annotation from %r", args.vcf)
    with pysam.VariantFile(args.vcf) as handle:
        if args.assembly_report is not None:
            log.info("Reading genome assembly info from %r", args.assembly_report)
            mapping = read_assembly_report(args.assembly_report)
        else:
            mapping = dict(zip(handle.header.contigs, handle.header.contigs))

        print(DBSNP_HEADER, end="")

        def _grouper(record):
            return (record.contig, record.pos)

        for _, records in itertools.groupby(handle.fetch(), key=_grouper):
            records = list(records)
            contig = mapping.get(records[0].contig)
            if contig is None:
                counter.skip(records[0].contig, len(records))
                continue

            # Group by contig:pos:ref so that each comparable SNP can be annotated with
            # information from similar SNPs/Indels at the same position
            by_ref = collections.defaultdict(list)
            for record in records:
                by_ref[record.ref].append(record)

            for records in by_ref.values():
                # Group by ALT string so that identical SNPs with multiple records are
                # merged into a single entry
                by_alt_string = collections.defaultdict(list)
                for record in records:
                    by_alt_string["/".join(sorted(record.alts or "."))].append(record)

                for alt_string, records in by_alt_string.items():
                    record = records[0]
                    info_string = _dbsnp_info_string(
                        alt_string=alt_string,
                        records=records,
                        alt_string_set=by_alt_string.keys(),
                    )

                    print(
                        TEMPLATE.format(
                            chrom=contig,
                            pos=record.pos + 1,
                            id=".",
                            ref=record.ref,
                            alt=",".join(record.alts or "."),
                            qual=".",
                            filter=".",
                            info=info_string,
                        )
                    )

                    counter(contig, len(records))


def _dbsnp_info_string(alt_string, records, alt_string_set):
    info = {}

    # Other alleles at the current site
    if len(alt_string_set) > 1:
        other_alts = alt_string_set - set((alt_string,))
        info["alts"] = ",".join(sorted(other_alts))

    ids = set()
    functions = set()
    for record in records:
        ids.add(record.id)

        for key, value in DBSNP_FUNCTIONS.items():
            if record.info.get(key):
                functions.add(value)

    info["ids"] = ",".join(ids)

    if functions:
        info["functions"] = ",".join(sorted(functions))

    return ";".join(f"{key}={value}" for key, value in info.items())


THOUSAND_GENOMES_FIELDS = set(
    (
        "AN_EAS_unrel",
        "AN_AMR_unrel",
        "AN_EUR_unrel",
        "AN_SAS_unrel",
        "AN_AFR_unrel",
        "AF_EUR_unrel",
        "AF_EAS_unrel",
        "AF_AMR_unrel",
        "AF_SAS_unrel",
        "AF_AFR_unrel",
    )
)


def thousand_genomes_to_vcf(log, counter, args):
    def _repr_value(value):
        if isinstance(value, float):
            value = "{:.7f}".format(value)
            if value == "0.0000000":
                value = "0"

        return str(value)

    for idx, filepath in enumerate(args.vcfs):
        log.info("Creating custom 1k Genomes VCF from %r", filepath)
        reduce_vcf_file(
            counter=counter,
            filepath=filepath,
            fields=THOUSAND_GENOMES_FIELDS,
            repr_value=_repr_value,
            print_header=idx == 0,
        )

    return 0


GNOMAD_SITE_FIELDS = set(
    (
        "AN",
        "AF",
        "AN_ami",
        "AF_ami",
        "AN_oth",
        "AF_oth",
        "AN_afr",
        "AF_afr",
        "AN_sas",
        "AF_sas",
        "AN_raw",
        "AF_raw",
        "AN_asj",
        "AF_asj",
        "AN_fin",
        "AF_fin",
        "AN_amr",
        "AF_amr",
        "AN_nfe",
        "AF_nfe",
        "AN_eas",
        "AF_eas",
    )
)


def gnomad_sites_to_vcf(log, counter, args):
    def _repr_value(value):
        if isinstance(value, float):
            value = "{:.5e}".format(value)
            if value == "0.00000e+00":
                value = "0"

        return str(value)

    for idx, filepath in enumerate(args.vcfs):
        log.info("Creating custom gnomAD sites VCF from %r", filepath)
        reduce_vcf_file(
            counter=counter,
            filepath=filepath,
            fields=GNOMAD_SITE_FIELDS,
            repr_value=_repr_value,
            print_header=idx == 0,
        )

    return 0


def gnomad_coverage_to_vcf(log, counter, args):
    log.info("Creating custom gnomAD filters from %r", args.txt)
    with open_ro(args.txt) as handle:
        header = handle.readline().rstrip().split("\t")

        # TODO: VCF header
        for line in handle:
            fields = line.rstrip().split("\t")
            row = dict(zip(header, fields))

            info = [
                "gnomAD_mean={}".format(row["mean"]),
                "gnomAD_median={}".format(row["median_approx"]),
                "gnomAD_over_15={}".format(row["over_15"]),
                "gnomAD_over_50={}".format(row["over_50"]),
            ]

            print(
                TEMPLATE.format(
                    chrom=row["#chr"],
                    pos=row["pos"],
                    id=".",
                    ref=".",
                    alt=".",
                    qual=".",
                    filter=".",
                    info=";".join(info),
                )
            )

            counter(row["#chr"])

    return 0


def reduce_vcf_file(counter, filepath, fields, repr_value=str, print_header=False):
    def _repr_values(values):
        strings = []
        if not isinstance(values, tuple):
            values = (values,)

        for value in values:
            strings.append(repr_value(value))

        return ",".join(strings)

    with pysam.VariantFile(filepath) as handle:
        if print_header:
            print("##fileformat=VCFv4.2")
            for value in handle.header.filters.values():
                print(value.record, end="")
            for name in fields:
                value = handle.header.info[name]
                print(value.record, end="")
            for value in handle.header.contigs.values():
                # Exclude uninformative header entries
                if value.length is not None:
                    print(value.header_record, end="")
            print("#CHROM  POS   ID REF  ALT   QUAL   FILTER INFO")

        for record in handle.fetch():
            infos = []
            for key in fields:
                value = record.info.get(key)
                if value is not None:
                    infos.append("{}={}".format(key, _repr_values(value)))

            print(
                TEMPLATE.format(
                    chrom=record.contig,
                    pos=record.pos,
                    id=record.id or ".",
                    ref=record.ref,
                    alt=",".join(record.alts or "."),
                    qual=record.qual or ".",
                    filter=",".join(record.filter or "."),
                    info=(";".join(infos) or "."),
                )
            )

            counter(record.contig)


def open_ro(filename: Union[str, Path]) -> IO[str]:
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(fspath(filename), "rb")
    try:
        header = handle.read(2)
        handle.seek(0)

        if header == b"\x1f\x8b":
            handle = cast(IO[bytes], gzip.GzipFile(mode="rb", fileobj=handle))
        elif header == b"BZ":
            handle = bz2.BZ2File(handle, "rb")

        return io.TextIOWrapper(handle)
    except:
        handle.close()
        raise


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv):
    parser = argparse.ArgumentParser(formatter_class=HelpFormatter)
    parser.set_defaults(command=None)
    subparsers = parser.add_subparsers()

    sub = subparsers.add_parser("assembly")
    sub.set_defaults(command=refseq_names_to_chr)
    sub.add_argument("assembly_report", metavar="FILE")

    sub = subparsers.add_parser("dbsnp")
    sub.set_defaults(command=dbsnp_to_vcf)
    sub.add_argument("vcf", metavar="FILE")
    sub.add_argument("assembly_report", nargs="?", metavar="FILE")

    sub = subparsers.add_parser("gnomad:cov")
    sub.set_defaults(command=gnomad_coverage_to_vcf)
    sub.add_argument("txt", metavar="FILE")

    sub = subparsers.add_parser("gnomad:sites")
    sub.set_defaults(command=gnomad_sites_to_vcf)
    sub.add_argument("vcfs", nargs="+", metavar="FILE")

    sub = subparsers.add_parser("1k_genomes")
    sub.set_defaults(command=thousand_genomes_to_vcf)
    sub.add_argument("vcfs", nargs="+", metavar="FILE")

    return parser


def main(argv):
    coloredlogs.install(
        level="INFO",
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(levelname)s %(message)s",
        format="%(asctime)s %(levelname)s %(message)s",
    )

    parser = parse_args(argv)
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_usage()
        return 1

    log = logging.getLogger("__main__")
    with Counter(log) as counter:
        return args.command(log, counter, args)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))