import nox

nox.options.sessions = [
    "style",
    "lints",
    "typing",
]


SOURCES = (
    "annovep",
    "bin/annovep",
    "noxfile.py",
    # "scripts",
)


RUFF_REQUIREMENT = "ruff==0.9.3"


@nox.session
def style(session: nox.Session) -> None:
    session.install(RUFF_REQUIREMENT)
    # Replaces `black --check`
    session.run("ruff", "format", "--check", *SOURCES)
    # Replaces `isort --check-only`
    session.run("ruff", "check", "--select", "I", *SOURCES)


@nox.session
def lints(session: nox.Session) -> None:
    session.install(RUFF_REQUIREMENT)
    session.run("ruff", "check", *SOURCES)


@nox.session()
def typing(session: nox.Session) -> None:
    session.install(".")
    session.install("nox~=2023.4.22")
    session.install("pyright==v1.1.393")
    session.run("pyright", *SOURCES)
