from __future__ import annotations

import bz2
import io
import os
import shlex
import signal
import subprocess
import sys
import time
from os import fspath
from typing import IO, TYPE_CHECKING, AnyStr, Iterable, cast

if TYPE_CHECKING:
    import logging
    from pathlib import Path

try:
    # Isal provides a significantly faster GZip implementation
    from isal.igzip import GzipFile
except ModuleNotFoundError:
    print(
        "WARNING: isal is not install; using built-in (slow) gzip reader\n"
        "         to install isal, run `python3 -m pip install isal`",
        file=sys.stderr,
    )
    from gzip import GzipFile


def cmd_to_str(command: Iterable[object], max_length: int | None = None) -> str:
    text = " ".join(shlex.quote(str(value)) for value in command)
    if max_length is not None and len(text) > max(0, max_length - 3):
        return text[: max_length - 3] + "..."

    return text


def update_required(output: str | Path, inputs: Iterable[str | Path]) -> bool:
    if not os.path.exists(output):
        return True

    mtime = os.stat(output).st_mtime
    return any(os.stat(dependency).st_mtime > mtime for dependency in inputs)


def open_rb(filename: str | Path) -> IO[bytes]:
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle: IO[bytes] = open(fspath(filename), "rb")

    try:
        header = handle.peek(2)

        if header.startswith(b"\x1f\x8b"):
            # ISA-l or builting GzipFile (see above)
            handle = cast(IO[bytes], GzipFile(mode="rb", fileobj=handle))
        elif header.startswith(b"BZ"):
            handle = bz2.BZ2File(handle, "rb")
    except:
        handle.close()
        raise
    else:
        return handle


def open_ro(filename: str | Path) -> IO[str]:
    return io.TextIOWrapper(open_rb(filename))


def join_procs(log: logging.Logger, procs: Iterable[subprocess.Popen[AnyStr]]) -> bool:
    start_time = time.time()
    sleep_time = 0.05
    commands = list(enumerate(procs))
    return_codes: list[int | None] = [None] * len(commands)

    assert all(hasattr(cmd, "args") for (_, cmd) in commands)

    while commands and not any(return_codes):
        try:
            # Wait for arbitrary command
            commands[0][1].wait(sleep_time if len(commands) > 1 else None)
        except subprocess.TimeoutExpired:
            sleep_time = min(1, sleep_time * 2)

        for index, command in list(commands):
            if command.poll() is not None:
                return_code = command.wait()
                return_codes[index] = return_code
                commands.remove((index, command))
                sleep_time = 0.05

                if return_code < 0:
                    return_code = signal.Signals(-return_code).name

                assert isinstance(command.args, Iterable)
                log.info("Command finished: %s", cmd_to_str(command.args, 80))
                log.info("  Runtime:        %.1fs", time.time() - start_time)
                log.info("  Return-code:    %s", return_code)

    if any(return_codes):
        for index, command in commands:
            assert isinstance(command.args, Iterable)
            log.warning("Terminating command: %s", cmd_to_str(command.args, 80))
            command.terminate()
            return_codes[index] = command.wait()

        log.error("Errors occured during processing!")

    return not any(return_codes)
