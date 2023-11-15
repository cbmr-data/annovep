from __future__ import annotations

import argparse
import logging
from typing import TYPE_CHECKING, Dict

from annovep.postprocess import output
from annovep.postprocess.annotations import Annotator
from annovep.postprocess.reader import VEPReader

if TYPE_CHECKING:
    from annovep.annotation import Annotation


def main(args: argparse.Namespace, annotations: list[Annotation]) -> int:
    if not any(fmt.startswith("sql") for fmt in args.output_format):
        args.include_json = False

    log = logging.getLogger("convert_vep")
    log.info("reading VEP annotations from '%s'", args.in_file)

    output_formats = set(args.output_format)
    if not output_formats:
        output_formats = ["tsv"]

    if args.out_prefix is None and len(output_formats) > 1:
        log.error("[out_prefix] must be set when writing more than one format")
        return 1

    vep_reader = VEPReader(args.in_file)

    annotator = Annotator(
        annotations=annotations,
        metadata=vep_reader.metadata,
        liftover_cache=args.data_liftover,
    )

    writers: Dict[str, output.Output] = {}
    for key in output_formats:
        cls = output.FORMATS[key]
        writer = cls(
            annotations=annotator,
            out_prefix=args.out_prefix,
        )

        writer.set_vcf_timestamp(args.vcf_timestamp)
        writer.set_vep_timestamp(vep_reader.timestamp)
        writers[key] = writer

    try:
        for record in vep_reader:
            if args.include_json:
                for writer in writers.values():
                    writer.process_json(record)

            for row in annotator.annotate(record):
                for writer in writers.values():
                    writer.process_row(row)

        for writer in writers.values():
            writer.finalize()
    except BrokenPipeError:
        pass

    return 0
