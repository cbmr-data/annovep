# WORKAROUND: pydantic requires List, Union, etc., so disable ruff lints:
# ruff: noqa: UP006,UP007,TC002
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from pydantic import BaseModel
from typing_extensions import Literal


class Args(BaseModel):
    in_file: Path
    out_prefix: Path
    transcript_strategy: Literal["canonical", "most-significant"]
    annotations: List[Path]
    annotations_list: bool
    enable: Dict[str, bool]
    do: Literal["run", "pre-process", "post-process"]
    vcf_timestamp: Optional[float]
    root: Path
    data_cache: Path
    data_custom: Path
    data_plugins: Path
    data_liftover: Path
    install: Path
    install_plugins: Path
    output_format: List[Literal["tsv", "json", "sql", "sqlite3"]]
    include_json: bool
    fork: int
    buffer_size: int
    log_level: Literal["debug", "info", "warning", "error"]


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter,
        prog="annovep pipeline",
    )

    # Pipeline
    parser.add_argument("in_file", type=Path)
    parser.add_argument("out_prefix", type=Path, nargs="?")

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

    group = parser.add_argument_group("Annotations")
    group.add_argument(
        "--annotations",
        metavar="FILE",
        type=Path,
        default=[],
        action="append",
        help="Optional files containing additional annotations",
    )

    group.add_argument(
        "--annotations-list",
        action="store_true",
        help="Quit after listing available annotation sources",
    )

    group.add_argument(
        "--enable",
        type=str.lower,
        metavar="NAME",
        default={},
        action=EnableAction,
        help="Enable annotations disabled by default; run annovep with "
        "--annotations-list to view available annotation sources",
    )

    group.add_argument(
        "--disable",
        dest="enable",
        default={},
        type=str.lower,
        metavar="NAME",
        action=DisableAction,
        help="Disable annotations enabled by default; run annovep with "
        "--annotations-list to view available annotation sources",
    )

    group.add_argument(
        "--transcript-strategy",
        metavar="STRATEGY",
        choices=("canonical", "most-significant"),
        default="most-significant",
        help="Strategy for selecting what transcript to report results for: Either the "
        "canonical transcript or the transcript for which the most significant "
        "consequence was observed",
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
        help="Installation folder for plugins; "
        "defaults to [$install/vep-plugins/Plugins]",
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


def parse_args(argv: list[str]) -> Args:
    parser = build_parser()
    args = parser.parse_args(argv)

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

    return Args.model_validate(args, from_attributes=True, strict=True)
