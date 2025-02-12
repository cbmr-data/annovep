from __future__ import annotations

import collections
import datetime
import json
import logging
import sqlite3
import sys
import zlib
from typing import TYPE_CHECKING, ClassVar, List, Sequence, cast

from typing_extensions import Literal, TypeAlias, TypedDict, override

from annovep._version import VERSION
from annovep.postprocess import consequences

if TYPE_CHECKING:
    from pathlib import Path

    from annovep.postprocess.annotations import Annotator
    from annovep.postprocess.reader import JSON

ConsequenceColumns: TypeAlias = Literal[
    "Func_most_significant",
    "Func_least_significant",
    "Func_most_significant_canonical",
]

# Columns that contain consequence terms (see `consequences.ranks()`)
CONSEQUENCE_COLUMNS = (
    "Func_most_significant",
    "Func_least_significant",
    "Func_most_significant_canonical",
)

UTC = datetime.timezone.utc


class GeneInfo(TypedDict):
    Chr: str
    MinPos: int
    MaxPos: int
    Variants: int
    Most_significant: int | None
    Most_significant_canonical: int | None


class Output:
    def __init__(
        self,
        annotations: Annotator,
        out_prefix: str | Path | None,
        extension: str,
    ) -> None:
        self.annotations = annotations
        self.fields = annotations.fields
        self._handle = sys.stdout
        if out_prefix is not None:
            self._handle = open(f"{out_prefix}{extension}", "w")

        self._date_vcf = "Unspecified"
        self._date_vep = "Unspecified"
        self._date_post = (
            datetime.datetime.now(tz=UTC).replace(microsecond=0).isoformat()
        )

    def set_vcf_timestamp(self, timestamp: float | None) -> None:
        if timestamp is not None:
            value = datetime.datetime.fromtimestamp(timestamp, tz=UTC)
            self._date_vcf = value.replace(microsecond=0).isoformat()
        else:
            self._date_vcf = "Unspecified"

    def set_vep_timestamp(self, timestamp: float) -> None:
        value = datetime.datetime.fromtimestamp(timestamp, tz=UTC)
        self._date_vep = value.replace(microsecond=0).isoformat()

    def finalize(self) -> None:
        if self._handle is not sys.stdout:
            self._handle.close()

    def process_json(self, data: JSON) -> None:
        pass

    def process_row(self, data: JSON) -> None:
        raise NotImplementedError

    def _print(self, line: str = "", *args: object) -> None:
        if args:
            line = line.format(*args)

        print(line, file=self._handle)


class JSONOutput(Output):
    def __init__(self, annotations: Annotator, out_prefix: str | Path | None) -> None:
        super().__init__(annotations, out_prefix, ".json")

    @override
    def process_row(self, data: JSON) -> None:
        json.dump(
            {field.output_key: data[field.output_key] for field in self.fields},
            self._handle,
        )
        self._handle.write("\n")


class TSVOutput(Output):
    def __init__(self, annotations: Annotator, out_prefix: str | Path | None) -> None:
        super().__init__(annotations, out_prefix, ".tsv")

        self._print("#{}", "\t".join([field.output_key for field in self.fields]))

        if out_prefix is not None:
            with open(f"{out_prefix}.tsv.columns", "w") as handle:
                print("Name\tDescription", file=handle)

                for field in self.fields:
                    print(field.output_key, field.help, sep="\t", file=handle)

    @override
    def process_row(self, data: JSON) -> None:
        row = [self._to_string(data[field.output_key]) for field in self.fields]

        self._print("\t".join(row))

    @staticmethod
    def _to_string(value: object) -> str:
        if isinstance(value, (tuple, list)):
            return ";".join(map(str, cast(List[object], value) or "."))
        elif value is None:
            return "."

        return str(value)


class SQLOutput(Output):
    # Columns that are renamed for the DB/shiny interface
    COLUMN_MAPPING: ClassVar[dict[str, str]] = {
        "Chr": "Hg38_chr",
        "Pos": "Hg38_pos",
    }

    def __init__(self, annotations: Annotator, out_prefix: str | Path | None) -> None:
        super().__init__(annotations, out_prefix, ".sql")

        self._consequence_ranks = self._build_consequence_ranks()
        self._contigs: dict[str, dict[str, int]] = {
            "hg19": collections.defaultdict(int),
            "hg38": collections.defaultdict(int),
        }
        self._genes: dict[str, GeneInfo] = {}
        self._n_overlap = 0
        self._n_row = 0
        self._n_json = 0

        self._print("PRAGMA TEMP_STORE=MEMORY;")
        self._print("PRAGMA JOURNAL_MODE=OFF;")
        self._print("PRAGMA SYNCHRONOUS=OFF;")
        self._print("PRAGMA LOCKING_MODE=EXCLUSIVE;")

        self._print("BEGIN;")
        self._print()
        self._print_descriptions()
        self._print()
        self._print_consequence_terms()
        self._print()
        self._print_gene_tables()
        self._print()

        for table in ("Annotations",):
            self._print("DROP TABLE IF EXISTS [{}];", table)
        self._print()

        query = [
            "CREATE TABLE [Annotations] (\n",
            "    [pk] INTEGER PRIMARY KEY ASC",
        ]

        for field in self.fields:
            key = field.output_key
            datatype = "TEXT"
            if key in CONSEQUENCE_COLUMNS:
                key = f"{key}_id"
                datatype = "INTEGER REFERENCES [Consequenes]([pk])"

            # Rename columns for SQL output only
            key = self.COLUMN_MAPPING.get(key, key)

            if field.type == "int":
                datatype = "INTEGER"
            elif field.type == "float":
                datatype = "REAL"
            else:
                assert field.type == "str", field.type

            query.append(f",\n    [{key}] {datatype}")
        query.append("\n);")
        self._print("".join(query))

        self._print()
        self._print_json_table()
        self._print()
        self._print("END;")
        self._print()
        self._print("BEGIN;")

    def finalize(self) -> None:
        self._print("END;")
        self._print()
        self._print("BEGIN;")

        for idx, (gene, info) in enumerate(sorted(self._genes.items())):
            self._print(
                "INSERT INTO [Genes] VALUES ({}, {}, {}, {}, {}, {}, {}, {});",
                idx,
                self._to_string(gene),
                self._to_string(info["Chr"]),
                self._to_string(info["MinPos"]),
                self._to_string(info["MaxPos"]),
                self._to_string(info["Variants"]),
                self._to_string(info["Most_significant"]),
                self._to_string(info["Most_significant_canonical"]),
            )

        self._print()
        self._print_contig_names()
        self._print()
        self._print_meta_data()
        self._print()
        self._print("CREATE INDEX IPositions_hg38 ON Annotations (Hg38_chr, Hg38_pos);")
        self._print("CREATE INDEX IPositions_hg19 ON Annotations (Hg19_chr, Hg19_pos);")
        self._print("CREATE INDEX IPositions_json ON JSON (Hg38_chr, Hg38_pos);")
        self._print("END;")
        self._print()

        # Collect information for the query planner
        self._print("ANALYZE;")

        self._handle.close()

    def process_json(self, data: JSON) -> None:
        data = dict(data)
        assert isinstance(data["input"], str)
        chrom, pos, _ = data["input"].split("\t", 2)

        # Keys are sorted to improve compression ratio
        blob = zlib.compress(json.dumps(data, sort_keys=True).encode("utf-8")).hex()

        self._print(
            "INSERT INTO [JSON] VALUES ({}, {}, {}, {});",
            self._n_json,
            self._to_string(chrom),
            int(pos),
            f"X'{blob}'",
        )

        self._n_json += 1

    @override
    def process_row(self, data: JSON) -> None:
        self._n_row += 1

        data = dict(data)
        assert isinstance(data["Chr"], str)
        assert isinstance(data["Pos"], int)
        self._contigs["hg38"][data["Chr"]] += 1
        if isinstance(data["Hg19_chr"], str):
            self._contigs["hg19"][data["Hg19_chr"]] += 1

        # Convert VEP consequence terms to ranks/numeric keys
        consequences: dict[ConsequenceColumns, int | None] = {
            key: self._consequence_ranks[data[key]] for key in CONSEQUENCE_COLUMNS
        }
        data.update(consequences.items())

        values = [str(self._n_row)]
        values.extend(self._to_string(data[field.output_key]) for field in self.fields)

        self._print("INSERT INTO [Annotations] VALUES ({});", ", ".join(values))

        overlapping_genes = data.get("Genes_overlapping")
        if overlapping_genes is not None:
            assert isinstance(overlapping_genes, list)

            for gene in overlapping_genes:
                assert isinstance(gene, str)

                gene_info = self._genes.get(gene)
                if gene_info is None:
                    self._genes[gene] = GeneInfo(
                        Chr=data["Chr"],
                        MinPos=data["Pos"],
                        MaxPos=data["Pos"],
                        Variants=1,
                        Most_significant=consequences["Func_most_significant"],
                        Most_significant_canonical=consequences[
                            "Func_most_significant_canonical"
                        ],
                    )
                elif gene_info["Chr"] != data["Chr"]:
                    raise ValueError(f"gene {gene!r} found on multiple contigs")
                else:
                    gene_info["MinPos"] = min(data["Pos"], gene_info["MinPos"])
                    gene_info["MaxPos"] = max(data["Pos"], gene_info["MaxPos"])
                    assert isinstance(gene_info["Variants"], int)
                    gene_info["Variants"] += 1

                    gene_info["Most_significant"] = self._worst_consequence(
                        gene_info["Most_significant"],
                        consequences["Func_most_significant"],
                    )

                    gene_info["Most_significant_canonical"] = self._worst_consequence(
                        gene_info["Most_significant_canonical"],
                        consequences["Func_most_significant_canonical"],
                    )

    @staticmethod
    def _worst_consequence(
        consequence_a: int | None,
        consequence_b: int | None,
    ) -> int | None:
        if consequence_a is None:
            return consequence_b
        elif consequence_b is None:
            return consequence_a

        return max(consequence_a, consequence_b)

    def _print_descriptions(self) -> None:
        self._print("DROP TABLE IF EXISTS [Columns];")
        self._print(
            """
            CREATE TABLE [Columns] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Name] TEXT,
              [Table] TEXT,
              [Column] TEXT,
              [Description] TEXT,
              [Type] TEXT,
              [ThousandsSep] TEXT,
              [Digits] INTEGER
            );
            """
        )

        for pk, field in enumerate(self.fields):
            # Rename columns for SQL output only
            key = self.COLUMN_MAPPING.get(field.output_key, field.output_key)

            table = "Annotations"
            column = key

            if key in CONSEQUENCE_COLUMNS:
                table = "Consequences"
                column = f"{key}_id"

            self._print(
                "INSERT INTO [Columns] VALUES ({}, {}, {}, {}, {}, {}, {}, {});",
                pk,
                self._to_string(key),
                self._to_string(table),
                self._to_string(column),
                self._to_string(field.help),
                self._to_string(field.type),
                self._to_string("," if field.thousands_sep else ""),
                field.digits,
            )

    def _print_consequence_terms(self) -> None:
        self._print("DROP TABLE IF EXISTS [Consequences];")
        self._print(
            """
            CREATE TABLE [Consequences] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Name] TEXT
            );
            """
        )

        for name, pk in self._consequence_ranks.items():
            if name is not None and pk is not None:
                self._print(
                    "INSERT INTO [Consequences] VALUES ({}, {});",
                    pk,
                    self._to_string(name),
                )

    def _print_gene_tables(self) -> None:
        self._print("DROP TABLE IF EXISTS [Genes];")
        self._print(
            """
            CREATE TABLE [Genes] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Name] TEXT,
              [Hg38_chr] TEXT,
              [Hg38_start] INTEGER,
              [Hg38_end] INTEGER,
              [Variants] INTEGER,
              [Most_significant] INTEGER REFERENCES [Consequenes]([pk]),
              [Most_significant_canonical] INTEGER REFERENCES [Consequenes]([pk])
            );
            """
        )

    def _print_json_table(self) -> None:
        self._print("DROP TABLE IF EXISTS [JSON];")
        self._print(
            """
            CREATE TABLE [JSON] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Hg38_chr] TEXT,
              [Hg38_pos] INTEGER,
              [Data] BINARY
            );
            """
        )
        self._print()

    def _print_contig_names(self) -> None:
        self._print("DROP TABLE IF EXISTS [Contigs];")
        self._print(
            """
            CREATE TABLE [Contigs] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Build] TEXT,
              [Name] TEXT,
              [Variants] INTEGER
            );
            """
        )

        contigs: list[tuple[str, str, int]] = []
        overlap = self._contigs["hg19"].keys() & self._contigs["hg38"].keys()

        # Collect hg38 contigs. These should be in the proper order
        for name, variants in self._contigs["hg38"].items():
            contigs.append(("hg38", name, variants))
            if name in overlap:
                contigs.append(("hg19", name, variants))

        # Collect hg19 only contigs; unmapped variants are ignored
        for name in sorted(self._contigs["hg19"].keys() - overlap - {None}):
            contigs.append(("hg19", name, self._contigs["hg19"][name]))

        for primary_key, (build, name, variants) in enumerate(contigs):
            self._print(
                "INSERT INTO [Contigs] VALUES ({}, {}, {}, {});",
                primary_key,
                self._to_string(build),
                self._to_string(name),
                self._to_string(variants),
            )

    def _print_meta_data(self) -> None:
        self._print("DROP TABLE IF EXISTS [Meta];")
        self._print(
            """
            CREATE TABLE [Meta] (
              [pk] INTEGER PRIMARY KEY ASC,
              [Label] TEXT,
              [Text] TEXT
            );
            """
        )

        metadata: list[tuple[str, int | str]] = [
            ("Version", VERSION),
            ("Date[VCF]", self._date_vcf),
            ("Date[Annotation]", self._date_vep),
            ("Date[Post-processing]", self._date_post),
            ("Variants", self._n_row),
        ]

        for idx, group in enumerate(self.annotations.groups, start=1):
            metadata.append((f"Annotation[{idx}]", group.name))

        for pk, (label, text) in enumerate(metadata):
            self._print(
                "INSERT INTO [Meta] VALUES ({}, {}, {});",
                pk,
                self._to_string(label),
                self._to_string(text),
            )

    @staticmethod
    def _build_consequence_ranks() -> dict[object, int | None]:
        """Returns consequences with a human friendly ranking: bad > insignificant."""
        human_friendly_ranks: dict[object, int | None] = collections.OrderedDict(
            (k, v) for v, k in enumerate(reversed(list(consequences.ranks())))
        )

        human_friendly_ranks[None] = None
        return human_friendly_ranks

    @staticmethod
    def _to_string(value: object) -> str:
        if isinstance(value, (int, float)):
            return repr(value)
        elif isinstance(value, (tuple, list)):
            if not value:
                return "NULL"

            value = ";".join(map(str, cast(Sequence[object], value)))
        elif value is None:
            return "NULL"
        else:
            value = str(value)

        return "'{}'".format(value.replace("'", "''"))


class SQLite3Output(SQLOutput):
    def __init__(self, annotations: Annotator, out_prefix: str | Path | None) -> None:
        self._conn = sqlite3.connect(f"{out_prefix}.sqlite3")
        self._curs = self._conn.cursor()

        super().__init__(annotations, None)

    def finalize(self) -> None:
        super().finalize()

        self._conn.commit()
        self._conn.close()

    def _print(self, line: str = "", *args: object) -> None:
        if args:
            line = line.format(*args)

        query = line.strip()
        if query:
            try:
                self._curs.execute(query)
            except sqlite3.OperationalError as error:
                pretty_query = " ".join(line.strip() for line in query.split("\n"))

                log = logging.getLogger(__name__)
                log.error("Error executing query: %s", error)
                log.error("  Query = %r", pretty_query)
                sys.exit(1)


FORMATS = {
    "json": JSONOutput,
    "tsv": TSVOutput,
    "sql": SQLOutput,
    "sqlite3": SQLite3Output,
}
