from __future__ import annotations

import argparse
import logging
import os
import subprocess
import sys
from typing import IO, AnyStr

from annovep.annotation import Annotation
from annovep.utils import cmd_to_str, join_procs, update_required


def exec(
    log: logging.Logger,
    command: list[str],
    stdin: None | int | IO[AnyStr] = subprocess.DEVNULL,
    stdout: None | int | IO[AnyStr] = None,
) -> subprocess.Popen[bytes]:
    log.info("Running %s", cmd_to_str(command))

    return subprocess.Popen(
        command,
        stdin=stdin,
        stdout=stdout,
        close_fds=True,
    )


def exec_self(
    log: logging.Logger,
    command: list[str],
    stdin: None | int | IO[AnyStr] = subprocess.DEVNULL,
    stdout: None | int | IO[AnyStr] = None,
) -> subprocess.Popen[bytes]:
    return exec(
        log=log,
        command=[sys.executable, "-m", "annovep", *command],
        stdin=stdin,
        stdout=stdout,
    )


def run_vep(
    args: argparse.Namespace,
    log: logging.Logger,
    annotations: list[Annotation],
) -> bool:
    command = [
        "vep",
        "--verbose",
        "--offline",
        "--cache",
        "--format",
        "vcf",
        "--json",
        "--force_overwrite",
        "--compress_output",
        "gzip",
        # Ensure that variants on unknown sequences are still written
        "--dont_skip",
        # Ensure that non-variant sites are still written
        "--allow_non_variant",
        "--dir_cache",
        args.data_cache,
        "--dir_plugins",
        args.install_plugins,
        "--assembly",
        "GRCh38",
        "--output_file",
        args.out_vep_json,
        "--stats_file",
        args.out_vep_html,
    ]

    if args.fork > 0:
        command.append("--fork")
        command.append(str(args.fork))

    if args.buffer_size > 0:
        command.append("--buffer_size")
        command.append(str(args.buffer_size))

    for annotation in annotations:
        command.extend(annotation.params)

    preproc = exec_self(
        log=log,
        command=[
            "--do",
            "pre-process",
            args.in_file,
            "--log-level",
            args.log_level,
        ],
        stdout=subprocess.PIPE,
    )

    vepproc = exec(log=log, command=command, stdin=preproc.stdout)

    return join_procs(log, [preproc, vepproc])


def run_post_proc(args: argparse.Namespace, log: logging.Logger) -> bool:
    command = [
        "--do",
        "post-process",
        args.out_vep_json,
        args.out_prefix,
        "--log-level",
        args.log_level,
        "--data-liftover",
        args.data_liftover,
        "--vcf-timestamp",
        repr(args.in_file.stat().st_mtime),
    ]

    if args.include_json:
        command.append("--include-json")

    if args.data_liftover is not None:
        command += ["--data-liftover", args.data_liftover]

    for annotation in args.annotations:
        command += ["--annotations", annotation]

    for name, enabled in args.enable.items():
        command += ["--enable" if enabled else "--disable", name]

    for fmt in args.output_format:
        command += ["--output-format", fmt]

    proc = exec_self(log=log, command=command)

    return join_procs(log, [proc])


def main(args: argparse.Namespace, annotations: list[Annotation]) -> int:
    log = logging.getLogger("annovep")

    any_errors = False
    for annotation in annotations:
        log.info("Checking files for annotation %s", annotation.name)
        for filename in annotation.files:
            if not os.path.exists(filename):
                log.error("Required %s file %r not found", annotation.name, filename)
                any_errors = True

    if any_errors:
        return 1

    args.out_vep_json = f"{args.out_prefix}.vep.json.gz"
    args.out_vep_html = f"{args.out_prefix}.vep.html"

    if update_required(
        output=args.out_vep_json,
        inputs=[args.in_file] + args.annotations,
    ):
        log.info("Running VEP")
        if not run_vep(args, log, annotations):
            return 1
    else:
        log.info("VEP annotations already up to date")

    log.info("Running post-processing")
    if not run_post_proc(args, log):
        return 1

    return 0
