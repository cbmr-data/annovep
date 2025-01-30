#!/usr/bin/env python3
import argparse
import collections
import io
import sys
import warnings
from dataclasses import dataclass
from os import fspath
from pathlib import Path
from typing import (
    IO,
    Any,
    Dict,
    FrozenSet,
    List,
    NamedTuple,
    NoReturn,
    Set,
    Tuple,
    Union,
    cast,
)

import ruamel.yaml
import ruamel.yaml.error
from dataclass_wizard import JSONWizard
from dataclass_wizard.errors import ParseError

try:
    from tqdm import tqdm
except ImportError:

    def tqdm(iterable: Any, **_kwargs: Any) -> Any:
        return iterable


try:
    from isal.igzip import GzipFile
except ModuleNotFoundError:
    warnings.warn("Module isal not found; using slow gzip reader")
    from gzip import GzipFile


SUMMARY = """
Scans annovep annotions for variants matching one or more models.

Models are specified in a YAML file, grouped by family. Variants matching a given model
is written to a per-family file, one row per matching model.

$ cat models.yaml
Family1:
  Members:
    Mother:   "Sample1"
    Father:   "Sample2"
    Sibling:  "Sample3"
    Proband:  "Sample4"

  Models:
    Recesssive:
      "1/1": Proband
      "0/1": [Mother, Father]
      "!1/1": Sibling

    Denovo:
      "0/1": Proband
      "0/0": [Mother, Father, Sibling]

$ annovep pipeline genotypes.vcf.gz annotations \
    --enable SampleGenotypes
$ python3 select_candidate_variants.py \
    --annotations annotations.tsv \
    --models models.yaml \
    --output output_prefix
"""
#######################################################################################


def abort(msg: str, *args: Any) -> NoReturn:
    sys.stderr.write(msg.format(*args) + "\n")
    sys.exit(1)


def load_yaml(filepath: Path) -> object:
    yaml = ruamel.yaml.YAML(typ="safe", pure=True)
    yaml.version = (1, 1)  # type: ignore

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ruamel.yaml.error.MantissaNoDotYAML1_1Warning)

        with filepath.open() as handle:
            return yaml.load(handle)


#######################################################################################
# Duplicate code to allow script to be used stand-alone


# Duplicated from pipeline.postprocess.consequences
def consequence_ranks() -> Dict[str, int]:
    # https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html
    consequences = [
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "splice_region_variant",
        "incomplete_terminal_codon_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        # Consequences for NMD transcripts are ignored
        # "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "feature_elongation",
        "regulatory_region_variant",
        "feature_truncation",
        "intergenic_variant",
        ".",
    ]

    ranks: collections.OrderedDict[str, int] = collections.OrderedDict()
    for rank, consequence in enumerate(reversed(consequences)):
        ranks[consequence] = rank

    return ranks


# Duplicated from pipeline.utils
def open_ro(filename: Path):
    """Opens a file for reading, transparently handling
    GZip and BZip2 compressed files. Returns a file handle."""
    handle = open(fspath(filename), "rb")
    try:
        header = handle.peek(2)
        if header == b"\x1f\x8b":
            handle = GzipFile(mode="rb", fileobj=handle)

        return io.TextIOWrapper(cast(IO[bytes], handle))
    except:
        handle.close()
        raise


#######################################################################################


@dataclass
class RawFamily:
    members: Dict[str, str]
    models: Dict[str, Dict[str, Union[str, List[str]]]]


@dataclass
class RawProject(JSONWizard):
    data: Dict[str, RawFamily]


# Supported genotypes
GENOTYPES: FrozenSet[str] = frozenset(("0/0", "0/1", "1/1"))


class Genotype(NamedTuple):
    genotypes: FrozenSet[str]
    samples: Tuple[str, ...]

    @staticmethod
    def load(genotype: str, samples: Tuple[str, ...]) -> "Genotype":
        invert: bool = False
        allow_uncalled: bool = False
        while genotype:
            if genotype.startswith("!"):
                genotype = genotype[1:]
                invert = True
            elif genotype.startswith("?"):
                genotype = genotype[1:]
                allow_uncalled = True
            else:
                break

        if genotype not in GENOTYPES:
            abort("invalid genotype {!r}", genotype)
        genotypes = frozenset((genotype,))
        if invert:
            genotypes = GENOTYPES - genotypes
        if allow_uncalled:
            genotypes = genotypes | set(("./.",))

        return Genotype(
            genotypes=genotypes,
            samples=tuple(f"GTS_{name}" for name in samples),
        )


class Model(NamedTuple):
    name: str
    requirements: Tuple[Genotype, ...]

    @staticmethod
    def load(
        name: str,
        data: Dict[str, Union[str, List[str]]],
        mapping: Dict[str, str],
    ) -> "Model":
        requirements: list[Genotype] = []
        for genotype, members in data.items():
            if isinstance(members, str):
                members = [members]

            samples = tuple(mapping[name] for name in members)
            requirements.append(Genotype.load(genotype, samples))

        if not requirements:
            abort(f"Model {name} has no requirements")

        return Model(
            name=name,
            requirements=tuple(requirements),
        )

    def matches(self, row: Dict[str, str]) -> bool:
        return all(
            (row[sample] in req.genotypes)
            for req in self.requirements
            for sample in req.samples
        )


class Family(NamedTuple):
    name: str
    models: Tuple[Model, ...]
    members: Dict[str, str]

    @staticmethod
    def load(name: str, data: RawFamily) -> "Family":
        if not data.members:
            abort("Family {} has no 'Members'", name)
        elif not data.models:
            abort("Family {} has no 'Models'", name)

        unique_samples = collections.Counter(data.members.values())
        for member, count in unique_samples.items():
            if count > 1:
                abort("Sample {} has been used multiple times", member)

        requirement_members: Set[str] = set()
        for model in data.models.values():
            for members in model.values():
                if isinstance(members, str):
                    requirement_members.add(members)
                else:
                    requirement_members.update(members)

        missing_members = requirement_members - data.members.keys()
        if missing_members:
            missing_members = ", ".join(repr(member) for member in missing_members)
            abort("Unknown family member(s) in family {}: {}", name, missing_members)

        family = Family(
            name=name,
            members={
                f"GTS_{label}": f"GTS_{sample}"
                for label, sample in data.members.items()
            },
            models=tuple(
                Model.load(name=model, data=model_data, mapping=data.members)
                for model, model_data in data.models.items()
            ),
        )

        return family


def load_project(filepath: Path) -> List[Family]:
    data = load_yaml(filepath)
    if isinstance(data, dict):
        try:
            project = RawProject.from_dict({"data": data})
        except ParseError as error:
            error._json_object = None
            error.kwargs = {}
            abort(str(error))

        return [Family.load(key, value) for key, value in project.data.items()]

    abort("YAML file {} does not contain dictionary", filepath)


#######################################################################################


class TableReader:
    def __init__(self, filepath: Path, consequences: FrozenSet[str]) -> None:
        self._filepath = filepath
        self._consequences = consequences

        self._handle = open_ro(self._filepath)
        self.header = tuple(
            self._handle.readline().lstrip("#").rstrip("\r\n").split("\t")
        )

    def __enter__(self):
        return self

    def __iter__(self):
        handle = tqdm(self._handle, unit_scale=True, unit=" rows")

        for lineno, line in enumerate(handle, start=1):
            row = line.rstrip("\r\n").split("\t")
            if len(row) != len(self.header):
                abort(f"Wrong number of columns at {self._filepath}:{lineno}")

            row = dict(zip(self.header, row))
            if row["Func_most_significant"] in self._consequences:
                yield row

    def __exit__(self, type: object, value: object, traceback: object) -> bool:
        self._handle.close()
        return False


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args: Any, **kwargs: Any):
        kwargs.setdefault("width", 88)

        super().__init__(*args, **kwargs)


def parse_args(argv: List[str], consequences: Dict[str, int]):
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter,
    )

    group = parser.add_argument_group("Input/Output")
    group.add_argument(
        "--models",
        type=Path,
        required=True,
        help="YAML file describing families and models",
    )
    group.add_argument(
        "--annotations",
        type=Path,
        required=True,
        help="Table containing annotations and genotypes",
    )
    group.add_argument(
        "--output",
        required=True,
        help="Prefix for output tables",
    )

    group = parser.add_argument_group("Filters")
    consequence_choices = list(consequences)
    group.add_argument(
        "--greatest-consequence",
        default="transcript_ablation",
        metavar="NAME",
        choices=consequence_choices,
        help="Exclude consequences more significant (as judged by VEP) than this "
        "consequence",
    )
    group.add_argument(
        "--smallest-consequence",
        default="splice_region_variant",
        metavar="NAME",
        choices=consequence_choices,
        help="Exclude consequences less significant (as judged by VEP) than this "
        "consequence",
    )
    group.add_argument(
        "--consequence",
        metavar="NAME",
        choices=consequence_choices,
        help="Equivalent to using --smallest-consequence and --greatest-consequence"
        "with the same consequence, i.e. only select this consequence",
    )

    return parser.parse_args(argv)


def main(argv: List[str]):
    consequences = consequence_ranks()
    args = parse_args(argv, consequences)

    if args.consequence is not None:
        args.smallest_consequence = args.consequence
        args.greatest_consequence = args.consequence

    min_consequence = consequences[args.smallest_consequence]
    max_consequence = consequences[args.greatest_consequence]
    selected_consequences = frozenset(
        key
        for key, value in consequences.items()
        if min_consequence <= value <= max_consequence
    )

    project = load_project(args.models)

    with TableReader(args.annotations, selected_consequences) as reader:
        sample_names: set[str] = set()
        for family in project:
            sample_names.update(family.members.values())

        missing_samples: set[str] = sample_names - set(reader.header)
        if missing_samples:
            abort("Samples from {} not found: {}", args.models, missing_samples)

        common_columns = [key for key in reader.header if not key.startswith("GTS_")]

        handles: Dict[str, IO[str]] = {}
        headers: Dict[str, List[str]] = {}
        for family in project:
            header = headers[family.name] = [
                "Family",
                "Model",
            ]

            for label in family.members:
                header.append(label)

            header.extend(common_columns)

            filename = f"{args.output}.{family.name}.tsv"
            handle = open(filename, "w")
            handles[family.name] = handle
            print(*header, sep="\t", file=handle)

        for row in reader:
            for family in project:
                for model in family.models:
                    if model.matches(row):
                        handle = handles[family.name]
                        header = headers[family.name]

                        row["Family"] = family.name
                        row["Model"] = model.name
                        for label, sample in family.members.items():
                            row[label] = row[sample]

                        print(*(row[key] for key in header), sep="\t", file=handle)

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
