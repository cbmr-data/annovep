from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Sequence

import coloredlogs

from annovep.annotation import Annotation, AnnotationError, load_annotations
from annovep.pipeline import main as pipeline_main
from annovep.postprocess import main as postprocess_main
from annovep.preprocess import main as preprocess_main


# Enable annotation with `--enable Name`
class EnableAction(argparse.Action):
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str | Sequence[object] | None,
        option_string: str | None = None,
    ) -> None:
        assert isinstance(values, str)
        getattr(namespace, self.dest)[values.lower()] = True


# Enable annotation with `--disable Name`
class DisableAction(argparse.Action):
    def __call__(
        self,
        parser: argparse.ArgumentParser,
        namespace: argparse.Namespace,
        values: str | Sequence[object] | None,
        option_string: str | None = None,
    ) -> None:
        assert isinstance(values, str)
        getattr(namespace, self.dest)[values.lower()] = False


def filter_annotations(
    log: logging.Logger,
    annotations: list[Annotation],
    enabled: dict[str, bool],
) -> bool:
    enabled = dict(enabled)
    names: set[str] = set(it.name.lower() for it in annotations)
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


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(
        self,
        prog: str,
        indent_increment: int = 2,
        max_help_position: int = 24,
        width: int | None = 79,
    ) -> None:
        super().__init__(
            prog=prog,
            indent_increment=indent_increment,
            max_help_position=max_help_position,
            width=width,
        )


def parse_args(argv: list[str]) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter,
        prog="annovep pipeline",
    )

    # Pipeline
    parser.add_argument("in_file", type=Path)
    parser.add_argument("out_prefix", type=Path, nargs="?")

    parser.add_argument(
        "--annotations",
        metavar="FILE",
        type=Path,
        default=[],
        action="append",
        help="Optional files containing additional annotations",
    )

    parser.add_argument(
        "--enable",
        type=str.lower,
        metavar="ANNOTATION",
        default={},
        action=EnableAction,
        help="Enable annotations disabled by default",
    )

    parser.add_argument(
        "--disable",
        dest="enable",
        default={},
        type=str.lower,
        metavar="ANNOTATION",
        action=DisableAction,
        help="Disable annotations enabled by default",
    )

    parser.add_argument(
        "--do",
        type=str.lower,
        choices=("run", "pre-process", "post-process"),
        default="run",
        help=argparse.SUPPRESS,
    )

    # Allows the vcf time-stamp to be passed to post-procesing
    parser.add_argument(
        "--vcf-timestamp",
        type=float,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--root",
        metavar="DIR",
        type=Path,
        default=Path("~/annovep").expanduser(),
        help="The root location of the AnnoVEP install",
    )

    group = parser.add_argument_group("Data locations")
    group.add_argument(
        "--data-cache",
        metavar="DIR",
        type=Path,
        help="Location of VEP cache; defaults to [$root/cache]",
    )
    group.add_argument(
        "--data-custom",
        metavar="DIR",
        type=Path,
        help="Location of custom annotation files; defaults to [$root/custom]",
    )
    group.add_argument(
        "--data-plugins",
        metavar="DIR",
        type=Path,
        help="Location of plugin data; defaults to [$root/plugins]",
    )
    group.add_argument(
        "--data-liftover",
        metavar="DIR",
        type=Path,
        help="Location of liftover cache; defaults to [$root/liftover]",
    )

    group = parser.add_argument_group("Installation locations")
    group.add_argument(
        "--install",
        metavar="DIR",
        type=Path,
        help="Installation folder; defaults to [$root/install]",
    )
    group.add_argument(
        "--install-plugins",
        metavar="DIR",
        type=Path,
        help="Installation folder for plugins; defaults to [$install/vep-plugins/Plugins]",
    )

    group = parser.add_argument_group("Output")
    group.add_argument(
        "--output-format",
        default=[],
        action="append",
        type=str.lower,
        choices=("tsv", "json", "sql", "sqlite3"),
        help="Output format for aggregated annotations. Maybe be specified zero or "
        "more times. Defaults to TSV if not specified",
    )

    group.add_argument(
        "--include-json",
        action="store_true",
        help="Include JSON data in SQL output, excluding sample specific information",
    )

    group = parser.add_argument_group("VEP options")
    group.add_argument(
        "--fork",
        metavar="N",
        type=int,
        default=0,
        help="Use forking to improve VEP runtime",
    )
    group.add_argument(
        "--buffer-size",
        metavar="N",
        default=100_000,
        type=int,
        help="Number of VCF records read by VEP per cycle",
    )

    group = parser.add_argument_group("Logging")
    group.add_argument(
        "--log-level",
        default="info",
        choices=("debug", "info", "warning", "error"),
        type=str.lower,
        help="Log messages at the specified level. This option applies to the "
        "`--log-file` option and to log messages printed to the terminal.",
    )

    return parser


def main(argv: list[str]) -> int:
    parser = parse_args(argv)
    args = parser.parse_args(argv)

    coloredlogs.install(
        level=args.log_level,
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(name)s %(levelname)s %(message)s",
        # Workaround for coloredlogs disabling colors in docker containers
        isatty=sys.stderr.isatty(),
    )

    if args.out_prefix is None:
        args.out_prefix = args.in_file

    if args.data_cache is None:
        args.data_cache = args.root / "cache"
    if args.data_custom is None:
        args.data_custom = args.root / "custom"
    if args.data_plugins is None:
        args.data_plugins = args.root / "plugins"
    if args.data_liftover is None:
        args.data_liftover = args.root / "liftover"

    if args.install is None:
        args.install = args.root / "install"
    if args.install_plugins is None:
        args.install_plugins = args.install / "vep-plugins" / "Plugins"

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
    except AnnotationError as error:
        log.error("error while loading annotations: %s", error)
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
