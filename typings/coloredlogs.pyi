import logging
from typing import IO, TypedDict

class Style(TypedDict, total=False):
    color: str
    bold: bool

def install(
    level: int | str | None = None,
    logger: logging.Logger = ...,
    fmt: str = ...,
    datefmt: str = ...,
    milliseconds: bool = ...,
    level_styles: dict[str, Style] = ...,
    field_styles: dict[str, Style] = ...,
    stream: IO[str] = ...,
    isatty: bool = ...,
    reconfigure: bool = ...,
    use_chroot: bool = ...,
    programname: str | None = None,
    syslog: bool = ...,
) -> None: ...
