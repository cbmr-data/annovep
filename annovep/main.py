# WORKAROUND: pydantic requies List, Union, etc, so disable ruff lints:
# ruff: noqa: UP006,UP007
from __future__ import annotations

import logging
import sys

import coloredlogs

from annovep.annotation import Annotation, AnnotationError, load_annotations
from annovep.args import parse_args
from annovep.pipeline import main as pipeline_main
from annovep.postprocess import main as postprocess_main
from annovep.preprocess import main as preprocess_main


def filter_annotations(
    log: logging.Logger,
    annotations: list[Annotation],
    enabled: dict[str, bool],
) -> bool:
    enabled = dict(enabled)
    names: set[str] = {it.name.lower() for it in annotations}
    set_all = enabled.pop("*", None)
    if set_all is not None:
        enabled.update(dict.fromkeys(names, set_all))

    result: list[Annotation] = []
    for annotation in annotations:
        name = annotation.name

        if annotation.enabled == "mandatory":
            result.append(annotation)
        elif enabled.get(name.lower(), annotation.enabled):
            log.info("   [✓] %s", name)
            result.append(annotation)
        else:
            log.warning("[☓] %s", name)

    annotations[:] = result
    unknown_annotations = enabled.keys() - names
    for unknown in unknown_annotations:
        log.error("--enable/--disable on unknown annotation: %r", unknown)

    return not unknown_annotations


def main(argv: list[str]) -> int:
    args = parse_args(argv)

    coloredlogs.install(
        level=args.log_level,
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(name)s %(levelname)s %(message)s",
        # Workaround for coloredlogs disabling colors in docker containers
        isatty=sys.stderr.isatty(),
    )

    variables = {
        # Data folders
        "data-cache": args.data_cache,
        "data-custom": args.data_custom,
        "data-plugins": args.data_plugins,
        "data-liftover": args.data_liftover,
        # Installation folders
        "install": args.install,
        "install-plugins": args.install_plugins,
    }

    log = logging.getLogger("annovep")
    try:
        annotations = load_annotations(args.annotations, variables)
    except AnnotationError:
        log.exception("error while loading annotations from %s", args.annotations)
        return 1

    if not filter_annotations(log, annotations, args.enable):
        return 1

    if args.do == "run":
        return pipeline_main(args, annotations)
    elif args.do == "pre-process":
        return preprocess_main(args, annotations)
    elif args.do == "post-process":
        return postprocess_main(args, annotations)

    raise NotImplementedError(args.do)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
