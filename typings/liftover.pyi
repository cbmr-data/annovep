from pathlib import Path

class ChainFile:
    def query(self, chrom: str, pos: int) -> list[tuple[str, int, bool]]: ...

def get_lifter(
    target: str,
    query: str,
    cache: str | Path | None = None,
) -> ChainFile: ...
