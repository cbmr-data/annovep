#!/usr/bin/env python3
# -*- coding: utf8 -*-
import logging
import os
import pwd
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Union

import coloredlogs

OK = "[✓]"
ERR = "[☓]"

MARKER = "not_mounted"

USER_ROOT = Path("/") / "data" / "user"
CACHE_ROOT = Path("/") / "data" / "cache"
ANNOVEP_ROOT = Path("/opt/annovep")

VEP_ROOT = Path("/opt/vep/src/ensembl-vep")
VEP_PLUGINS = Path("/opt/vep-plugins/Plugins")


COMMANDS: Dict[str, List[Union[Path, str]]] = {
    "bash": ["/bin/bash"],
    # VEP setup and direct execution
    "vep": ["perl", VEP_ROOT / "vep"],
    "setup": ["bash", ANNOVEP_ROOT / "scripts" / "setup_vep.sh"],
    # Main pipeline
    "pipeline": [
        "python3",
        ANNOVEP_ROOT / "pipeline" / "main.py",
        "--root",
        CACHE_ROOT,
        "--install-plugins",
        VEP_PLUGINS,
        "--annotations",
        ANNOVEP_ROOT / "annotations",
    ],
}


def quote_command(args):
    return [shlex.quote(str(value)) for value in args]


def user_exists(uid):
    try:
        pwd.getpwuid(uid)
    except KeyError:
        return False

    return True


def user_add(uid, name):
    return not subprocess.call(["useradd", "-u", str(uid), name])


def check_container(log):
    log.info("Checking container ..")

    folders = [
        ["cache", CACHE_ROOT],
        ["user data", USER_ROOT],
    ]

    for name, root in folders:
        if not root.is_dir():
            log.error(" %s annovep %s folder does not exist: %s", ERR, name, root)
            return False
        log.info("  %s annovep %s folder exists ..", OK, name)

        if (root / MARKER).exists():
            log.error(" %s annovep %s folder not mounted in podman!", ERR, name)
            return False
        log.info("  %s annovep %s folder mounted ..", OK, name)

    return True


def check_databases(log):
    log.info("VEP cache accessible at '%s'", CACHE_ROOT)

    return True


def main(argv):
    coloredlogs.install(
        level="INFO",
        datefmt="%Y-%m-%d %H:%M:%S",
        fmt="%(asctime)s %(name)s %(levelname)s %(message)s",
    )

    log = logging.getLogger("annovep")
    if not argv:
        log.error("No command specified. Available commands are ")
        for nth, key in enumerate(COMMANDS, start=1):
            log.error("  %i. %s", nth, key)

        return 1

    command, *argv = argv
    commandline = COMMANDS.get(command.lower())
    if commandline is None:
        log.error("unknown command %r", command)
        return 1

    if not check_container(log):
        if command != "bash":
            return 1

    if not check_databases(log):
        if command != "bash":
            return 1

    # Set home to user data folder
    env = dict(os.environ)

    env["HOME"] = "/home/dummy"

    env["ANNOVEP_CACHE"] = str(CACHE_ROOT)

    env["VEP_ROOT"] = str(VEP_ROOT)

    log.info("Running %s", " ".join(quote_command(commandline + argv)))

    return subprocess.call(
        commandline + argv,
        stdin=sys.stdin,
        cwd=USER_ROOT,
        env=env,
    )


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
