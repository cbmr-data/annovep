#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import bz2
import collections
import fnmatch
import io
import itertools
import logging
import sys
import zipfile
from os import fspath
from pathlib import Path
from shlex import quote
from typing import IO, NamedTuple, Optional, Union, cast

import coloredlogs
import pysam
from ruamel.yaml import YAML

try:
    # ISA-L is significantly faster than the built-in gzip decompressor
    from isal.igzip import GzipFile
except ModuleNotFoundError:
    from gzip import GzipFile

VCF_HEADER = "##fileformat=VCFv4.2"
VCF_COLUMN_NAMES = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
VCF_ROW_TEMPLATE = "{chrom}\t{pos}\t{id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}"


class Counter:
    STEP = 10e6

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
##INFO=<ID=ids,Number=.,Type=String,Description="List of rsIDs for this SNP.">
##INFO=<ID=alts,Number=.,Type=String,Description="List of sets of alleles observed for a given chr:pos:ref (e.g. 'A,A/T,A/T/G').">
##INFO=<ID=functions,Number=.,Type=String,Description="List of GO terms assosiated with this SNP in DBSNP.">"""


def dbsnp_to_vcf(log, counter, args):
    log.info("Creating custom DBSNP annotation from %r", args.vcf)
    with pysam.VariantFile(args.vcf, threads=2) as handle:
        if args.assembly_report is not None:
            log.info("Reading genome assembly info from %r", args.assembly_report)
            mapping = read_assembly_report(args.assembly_report)
        else:
            mapping = dict(zip(handle.header.contigs, handle.header.contigs))

        # Some fields are never used; this saves some time
        template = VCF_ROW_TEMPLATE.format(
            chrom="{chrom}",
            pos="{pos}",
            id=".",
            ref="{ref}",
            alt="{alt}",
            qual=".",
            filter=".",
            info="{info}",
        )

        info_printed = False
        for line in str(handle.header).splitlines():
            if line.startswith("##INFO"):
                if not info_printed:
                    print(DBSNP_HEADER)
                    info_printed = True
            else:
                print(line)

        def _grouper(record):
            return (record.contig, record.pos)

        for (contig, pos), records in itertools.groupby(handle, key=_grouper):
            mapped_contig = mapping.get(contig)
            if mapped_contig is None:
                counter.skip(contig, sum(1 for _ in records))
                continue

            contig = mapped_contig
            # ALT strings for a given REF; one or more SNPs/INDELS or None
            alt_strings = collections.defaultdict(list)
            # DBSNP may contain multiple records for the same REF/ALT combination,
            # if a SNP has been registered multiple times. These are combined into
            # a single record with relevant information aggregated
            duplicate_records = collections.defaultdict(list)

            nrecords = 0
            for record in records:
                ref = record.ref
                alts = tuple(record.alts or ".")

                nrecords += 1
                alt_strings[ref].append("/".join(alts))
                duplicate_records[(ref, alts)].append(record)

            for (ref, alts), records in duplicate_records.items():
                print(
                    template.format(
                        chrom=contig,
                        pos=pos,
                        ref=ref,
                        alt=",".join(alts),
                        info=_dbsnp_info_string(
                            records=records,
                            alt_string_set=alt_strings[ref],
                        ),
                    )
                )

            counter(contig, nrecords)


def _dbsnp_info_string(records, alt_string_set):
    ids = set()
    info_keys = set()
    for record in records:
        assert record.id not in ids
        ids.add(record.id)
        info_keys.update(record.info)

    # All ALT allele strings at the current site with matching REF
    info = [
        "alts=",
        ",".join(alt_string_set),
        ";ids=",
        ",".join(ids),
    ]

    functions = info_keys & DBSNP_FUNCTIONS.keys()
    if functions:
        info.append(";functions=")
        info.append(",".join(DBSNP_FUNCTIONS[key] for key in functions))

    return "".join(info)


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

        # Some fields are never used; this saves some time
        template = VCF_ROW_TEMPLATE.format(
            chrom="{chrom}",
            pos="{pos}",
            id=".",
            ref=".",
            alt=".",
            qual=".",
            filter=".",
            info="{info}",
        )

        # TODO: VCF header
        for line in handle:
            fields = line.rstrip().split("\t")
            row = dict(zip(header, fields))

            chrom, pos = row["locus"].split(":")
            info = "mean={};median={};over_15={};over_50={}".format(
                row["mean"],
                row["median_approx"],
                row["over_15"],
                row["over_50"],
            )

            print(
                template.format(
                    chrom=chrom,
                    pos=pos,
                    info=info,
                )
            )

            counter(chrom)

    return 0


########################################################################################


class GeneRecord(NamedTuple):
    seqid: str
    start: int
    end: int
    name: Optional[str]
    id: Optional[str]

    @property
    def preferred_name(self):
        return self.name or self.id

    def __repr__(self):
        return f"{self.seqid}:{self.start}-{self.end}:{self.preferred_name}"


def neighbouring_genes_to_bed(log, counter, args):
    def _groupby(record):
        return record.seqid

    forward_and_reverse = (
        ("over", iter_nearest_genes_overlapping),
        ("up", iter_nearest_genes_upstream),
        ("down", iter_nearest_genes_downstream),
    )

    genes = read_genes_from_gff(log, args.gff)
    for seqid, genes in itertools.groupby(genes, _groupby):
        genes = list(genes)

        records = []
        for key, func in forward_and_reverse:
            for start, end, nearest in func(genes, nnearest=3):
                records.append((start, end, key, nearest))

        for start, end, key, nearest in sorted(records):
            nearest = ";".join(
                f"{it.start}-{it.end}:{it.preferred_name}" for it in nearest
            )

            # GFF uses 1-based coordiantes for both start and end, while BED uses 0
            # for the start coordinates and 1 for the end coordinates
            print(f"{seqid}\t{start - 1}\t{end}\t{key};{nearest}")

        counter(seqid, len(genes))

    return 0


def iter_nearest_genes_overlapping(genes, nnearest=3):
    for gene in genes:
        yield gene.start, gene.end, (gene,)


def iter_nearest_genes_downstream(genes, nnearest=3):
    genes.sort(key=lambda it: it.start)

    last_position = 1
    for idx, gene in enumerate(genes):
        nearest = genes[idx : idx + nnearest]

        yield last_position, gene.start - 1, nearest

        last_position = gene.start


def iter_nearest_genes_upstream(genes, nnearest=3):
    genes.sort(key=lambda it: it.end)

    last_position = genes[0].end + 1
    for idx, gene in enumerate(genes[1:]):
        nearest = genes[max(0, idx - nnearest + 1) : idx + 1]

        yield last_position, gene.end, nearest

        last_position = gene.end + 1

    # The last region must cover everything downstream of the final gene. But since we
    # don't know how big the contig is, simply use the largest value supported by tabix.
    yield last_position, 2**29 - 1, genes[-3:]


def read_genes_from_gff(log, filename):
    with open_ro(filename) as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            seqid, _, kind, start, end, _, _, _, attributes = line.split("\t")

            if kind == "gene":
                attrs = {}
                for attribute in attributes.split(";"):
                    key, value = attribute.split("=")
                    attrs[key.strip()] = value.strip()

                gene_id = attrs.get("gene_id")
                gene_name = attrs.get("Name")

                if gene_name or gene_id:
                    yield GeneRecord(
                        seqid=seqid,
                        start=int(start),
                        end=int(end),
                        name=gene_name,
                        id=gene_id,
                    )
                else:
                    log.warning("skipping %r", line)


########################################################################################


def dbnsfp4_to_vcf(log, counter, args):
    log.info("loading annotation list from %r", str(args.annotations))
    yaml = YAML(typ="safe")
    with open(args.annotations) as handle:
        data = yaml.load(handle)

    columns = []
    for settings in data.values():
        columns.extend(settings["Fields"])

    log.info("loading annotations from %r", str(args.zip))
    with zipfile.ZipFile(args.zip) as zhandle:
        print(VCF_HEADER)
        print(VCF_COLUMN_NAMES)

        for filename in fnmatch.filter(zhandle.namelist(), "dbNSFP4*_variant.chr*.gz"):
            log.info("extracting annotation file %r", str(filename))
            with zhandle.open(filename) as gzhandle:
                with GzipFile(fileobj=gzhandle) as handle:
                    handle = io.TextIOWrapper(handle)
                    header = handle.readline().rstrip().split("\t")

                    for line in handle:
                        row = dict(zip(header, line.rstrip().split("\t")))

                        info = []
                        for key in columns:
                            value = row[key]
                            if value not in (".", ".;", "./."):
                                if ";" in value:
                                    # Workaround to handle multiple values in a field:
                                    # ASCII comma is replaced with Unicode Small Comma
                                    value = ",".join(
                                        value.replace(",", "﹐").split(";")
                                    )

                                if info:
                                    info.append(";")
                                info.append(key)
                                info.append("=")
                                info.append(value)

                        info = "".join(info)

                        out = [
                            row["#chr"],
                            row["pos(1-based)"],
                            ".",
                            row["ref"],
                            row["alt"],
                            ".",
                            ".",
                            info,
                        ]

                        print("\t".join(out))

                        counter(row["#chr"])


########################################################################################


def reduce_vcf_file(counter, filepath, fields, repr_value=str, print_header=False):
    def _repr_values(values):
        strings = []
        if not isinstance(values, tuple):
            values = (values,)

        for value in values:
            strings.append(repr_value(value))

        return ",".join(strings)

    with pysam.VariantFile(filepath, threads=2) as handle:
        if print_header:
            print(VCF_HEADER)
            for value in handle.header.filters.values():
                print(value.record, end="")
            for name in fields:
                value = handle.header.info[name]
                print(value.record, end="")
            for value in handle.header.contigs.values():
                # Exclude uninformative header entries
                if value.length is not None:
                    print(value.header_record, end="")
            print(VCF_COLUMN_NAMES)

        for record in handle.fetch():
            infos = []
            for key in fields:
                value = record.info.get(key)
                if value is not None:
                    infos.append("{}={}".format(key, _repr_values(value)))

            print(
                VCF_ROW_TEMPLATE.format(
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
        header = handle.peek(2)

        if header.startswith(b"\x1f\x8b"):
            handle = cast(IO[bytes], GzipFile(mode="rb", fileobj=handle))
        elif header.startswith(b"BZ"):
            handle = bz2.BZ2File(handle, "rb")

        return io.TextIOWrapper(handle)
    except:
        handle.close()
        raise


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("width", 79)

        super().__init__(*args, **kwargs)


def parse_args(argv: list[str]) -> argparse.ArgumentParser:
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

    sub = subparsers.add_parser("neighbours")
    sub.set_defaults(command=neighbouring_genes_to_bed)
    sub.add_argument("gff", metavar="FILE")

    sub = subparsers.add_parser("dbnsfp4")
    sub.set_defaults(command=dbnsfp4_to_vcf)
    sub.add_argument("annotations", metavar="YAML")
    sub.add_argument("zip", metavar="ZIP")

    return parser


def main(argv: list[str]) -> int:
    coloredlogs.install(
        level="INFO",
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(levelname)s %(message)s",
    )

    log = logging.getLogger("__main__")

    try:
        from isal.igzip import GzipFile
    except ModuleNotFoundError:
        log.warning("Python module 'isal' not installed; will use slow gzip reader")
        log.warning("To install, run `%s -m pip install isal`", quote(sys.executable))

    parser = parse_args(argv)
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_usage()
        return 1

    # Silence log-messages from HTSLIB
    pysam.set_verbosity(0)

    with Counter(log) as counter:
        try:
            return args.command(log, counter, args)
        except BrokenPipeError:
            return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
