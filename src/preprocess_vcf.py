#!/usr/bin/env python3
# -*- coding: utf8 -*-
import argparse
import bz2
import gzip
import io
import logging
import re
import sys
from os import fspath
from pathlib import Path

import coloredlogs


_RE_CONTIG_ID = re.compile("^(##contig=<.*ID=)([^,]+)(.*>)$", re.I)


def encode_contig_name(name):
    """Reversible encoding of contig names that cause problems with VEP."""
    if ":" in name or "*" in name:
        return "annovep_{}".format(name.encode("utf-8").hex())

    return name


def fix_contig_name(line):
    match = _RE_CONTIG_ID.match(line)
    if match is not None:
        before, name, after = match.groups()
        new_name = encode_contig_name(name)

        return "".join((before, new_name, after, "\n")), name, new_name

    return line, None, None


def open_ro(filename):
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(fspath(filename), "rb")
    try:
        header = handle.read(2)
        handle.seek(0)

        if header == b"\x1f\x8b":
            handle = gzip.GzipFile(mode="rb", fileobj=handle)
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
    parser.add_argument("in_vcf", type=Path)

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    coloredlogs.install(
        fmt="%(asctime)s %(name)s %(levelname)s %(message)s",
        # Workaround for coloredlogs disabling colors in docker containers
        isatty=sys.stderr.isatty(),
    )

    n_decoys = 0
    n_bad_contigs = 0
    n_bad_contigs_vcf = 0
    n_decoy_contigs = 0
    n_records = 0

    log = logging.getLogger("preprocess_vcf")
    log.info("Reading VCF at %s", args.in_vcf)
    with open_ro(args.in_vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                line, old_name, new_name = fix_contig_name(line)
                if old_name != new_name:
                    n_bad_contigs += 1
                    if n_bad_contigs == 1:
                        log.warning("Changing bad names: %r -> %r", old_name, new_name)

                if old_name and old_name.endswith("_decoy"):
                    n_decoy_contigs += 1
                    if n_decoy_contigs == 1:
                        log.warning("Removing decoy contigs: %r", old_name)
                    continue

                print(line, end="")
                continue

            n_records += 1
            chrom, rest = line.split("\t", 1)
            if chrom.endswith("_decoy"):
                n_decoys += 1
                if n_decoys == 1:
                    log.warning("Filtering variants on decoy contigs (e.g. %r)", chrom)

            new_chrom = encode_contig_name(chrom)
            if chrom != new_chrom:
                n_bad_contigs_vcf += 1

            print(new_chrom, rest, sep="\t", end="")

            if n_records % 100_000 == 0:
                log.info("Read %s variants; at %r", "{:,}".format(n_records), chrom)

    def _fmt(value):
        return "{:,}".format(value)

    log.info("Read %s variants", _fmt(n_records))
    log.info("Renamed %s contigs with bad names", _fmt(n_bad_contigs))
    log.info("Updated %s variants on badly named contigs", _fmt(n_bad_contigs_vcf))
    log.info("Removed %s decoy contigs", _fmt(n_decoy_contigs))
    log.info("Removed %s variants on decoy contigs", _fmt(n_decoys))

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
