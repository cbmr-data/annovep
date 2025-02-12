# WORKAROUND: pydantic requies List, Union, etc, so disable ruff lints:
# ruff: noqa: UP006,UP007
from __future__ import annotations

import logging
from dataclasses import dataclass
from fnmatch import fnmatch
from pathlib import Path
from typing import Dict, List, Optional, Union

import pydantic
import ruamel.yaml
from pydantic import ConfigDict, Field, RootModel, ValidationError
from typing_extensions import Annotated, Literal, TypeAlias, override

from annovep.resources import access_resources

# Built-in annotations derived from input
_BUILTINS = {
    "samplegenotypes": "SampleGenotypes",
    "liftover": "Liftover",
}

DataSource: TypeAlias = Literal[
    "canonical_consequence",
    "consequence",
    "derived",
    "genotypes",
    "input",
    "liftover",
]


FieldType: TypeAlias = Literal["str", "int", "float"]

VariablesType: TypeAlias = Dict[str, Union[str, Path]]


class AnnotationError(Exception):
    pass


########################################################################################
# Annotation parsing and validation


class BaseModel(pydantic.BaseModel):
    model_config = ConfigDict(extra="forbid")


class _AnnotationFieldModel(BaseModel):
    input_key: str = Field(alias="Input")
    output_key: Optional[str] = Field(alias="Output", default=None)
    fieldtype: Optional[FieldType] = Field(alias="FieldType", default=None)
    help: Optional[str] = Field(alias="Help", default=None)
    # Normalize lists of items using this separator
    split_by: Optional[str] = Field(alias="SplitBy", default=None)
    # If split_by is set, also sort the resulting list
    sort: bool = Field(alias="Sort", default=False)
    # Enable thousands separator
    thousands_sep: bool = Field(alias="ThousandsSep", default=False)
    # Floating point precision (int/float only)
    digits: Optional[int] = Field(alias="Digits", default=None)
    # Optional data source
    source: Optional[DataSource] = Field(alias="Source", default=None)
    # Optional glob for output keys
    derived_from: Optional[Union[str, List[str]]] = Field(
        alias="DerivedFrom", default=None
    )


class _AnnotationBaseModel(BaseModel):
    rank: int = Field(alias="Rank", default=0)
    fieldtype: FieldType = Field(alias="FieldType", default="str")
    digits: int = Field(alias="Digits", default=-1)
    fields: List[_AnnotationFieldModel] = Field(alias="Fields", default_factory=list)
    enabled: Literal[True, False, "mandatory"] = Field(alias="Enabled", default=True)
    options: List[str] = Field(alias="Options", default_factory=list)
    # Optional data source
    source: Optional[DataSource] = Field(alias="Source", default=None)

    def to_annotation(
        self,
        *,
        name: str,
        variables: VariablesType,
    ) -> Annotation:
        return Annotation(
            type=self._get_type(),
            name=name,
            rank=self.rank,
            fields=self._get_fields(name=name),
            enabled=self.enabled,
            files=self._get_files(variables=variables),
            params=self._get_options(name=name, variables=variables),
        )

    def _get_fields(self, name: str) -> list[AnnotationField]:
        derived_fields: list[tuple[AnnotationField, str | list[str]]] = []

        log = logging.getLogger(__name__)
        fields: list[AnnotationField] = []
        for it in self.fields:
            if it.input_key == it.output_key:
                log.warning("Redundant Output key in %s:%s", name, it.input_key)
            elif it.derived_from and it.input_key not in (":min:", ":max:"):
                log.error("DerivedFrom specified but Input is %r", it.input_key)
                raise AnnotationError(it)

            field = AnnotationField(
                input_path=self._get_path(
                    group=name,
                    field=it.input_key,
                    source=self.source if it.source is None else it.source,
                ),
                output_key=it.input_key if it.output_key is None else it.output_key,
                type=self.fieldtype if it.fieldtype is None else it.fieldtype,
                derived_from=[],
                help=it.help,
                split_by=it.split_by,
                sort=it.sort,
                thousands_sep=it.thousands_sep,
                digits=self.digits if it.digits is None else it.digits,
            )

            fields.append(field)
            if it.derived_from:
                derived_fields.append((field, it.derived_from))

        for it, derived_from in derived_fields:
            if isinstance(derived_from, str):
                for other in fields:
                    if other is not it and fnmatch(other.output_key, derived_from):
                        it.derived_from.append(other.output_key)
            else:
                names = {other.output_key for other in fields}
                unknown = set(derived_from) - names
                if unknown:
                    log.error("Unknown field names for %r: %s", it.output_key, unknown)
                    raise AnnotationError(it)

                it.derived_from = derived_from

            log.debug("%r derived from fields %r", it.output_key, it.derived_from)

        return fields

    def _get_type(self) -> AnnotationTypes:
        raise NotImplementedError

    def _get_files(self, *, variables: VariablesType) -> list[str]:
        return []

    def _get_options(self, *, name: str, variables: VariablesType) -> list[str]:
        return self.options

    def _get_path(self, *, group: str, field: str, source: str | None) -> list[str]:
        if source is None:
            raise AnnotationError(f"no source specified for {group}")

        return [source, field]


class _BasicAnnotationModel(_AnnotationBaseModel):
    type: Literal["Basic"] = Field(..., alias="Type")

    @override
    def _get_type(self) -> AnnotationTypes:
        return "basic"


class _PluginModel(_AnnotationBaseModel):
    type: Literal["Plugin"] = Field(..., alias="Type")
    files: List[str] = Field(alias="Files")
    parameters: List[str] = Field(alias="Parameters")
    variables: Dict[str, str] = Field(alias="Variables", default_factory=dict)

    @override
    def _get_type(self) -> AnnotationTypes:
        return "plugin"

    @override
    def _get_files(self, *, variables: VariablesType) -> list[str]:
        variables = dict(variables)
        variables.update(self.variables)

        return _apply_variables(self.files, variables)

    @override
    def _get_options(self, *, name: str, variables: VariablesType) -> list[str]:
        variables = dict(variables)
        variables.update(self.variables)
        parameters = _apply_variables(self.parameters, variables)

        return [*self.options, "--plugin", ",".join([name, *parameters])]

    @override
    def _get_path(self, *, group: str, field: str, source: str | None) -> list[str]:
        if source is not None:
            raise AnnotationError(f"source must not be specified for {group}")

        return ["consequence", field]


class _CustomModel(_AnnotationBaseModel):
    type: Literal["BED", "VCF"] = Field(..., alias="Type")
    mode: Literal["exact", "overlap"] = Field(..., alias="Mode")
    file: str = Field(..., alias="File")
    variables: Dict[str, str] = Field(alias="Variables", default_factory=dict)

    @override
    def _get_type(self) -> AnnotationTypes:
        return "custom"

    @override
    def _get_files(self, *, variables: VariablesType) -> list[str]:
        variables = dict(variables)
        variables.update(self.variables)

        return _apply_variables([self.file, f"{self.file}.tbi"], variables)

    @override
    def _get_options(self, *, name: str, variables: VariablesType) -> list[str]:
        file, _ = self._get_files(variables=variables)
        options = [*self.options, "--custom"]
        params = [file, name, self.type, self.mode, "0"]
        for field in self.fields:
            key = field.input_key
            if not (key.startswith(":") and key.endswith(":")):
                params.append(key)

        options.append(",".join(params))
        return options

    @override
    def _get_path(self, *, group: str, field: str, source: str | None) -> list[str]:
        if source is not None:
            raise AnnotationError(f"source must not be specified for {group}")

        return ["custom", group, field]


class _BuiltinModel(_AnnotationBaseModel):
    type: Literal["Builtin"] = Field(..., alias="Type")

    @override
    def to_annotation(
        self,
        *,
        name: str,
        variables: VariablesType,
    ) -> Annotation:
        builtin_name = _BUILTINS.get(name.lower())
        if builtin_name is None:
            raise ValueError(f"unknown built-in annotation {name!r}")

        return super().to_annotation(name=name, variables=variables)

    @override
    def _get_type(self) -> AnnotationTypes:
        return "builtin"


_Annotations = Annotated[
    Union[_BasicAnnotationModel, _PluginModel, _CustomModel, _BuiltinModel],
    Field(discriminator="type"),
]

_AnnotationsDict = Dict[str, _Annotations]


class _Root(RootModel[_AnnotationsDict]):
    root: _AnnotationsDict


########################################################################################
# Simplified user-facing models


@dataclass
class AnnotationField:
    input_path: list[str]
    output_key: str
    type: FieldType
    help: Optional[str]
    derived_from: list[str]
    # Normalize lists of items using this separator
    split_by: Optional[str] = None
    # If split_by is set, also sort the resulting list
    sort: bool = False
    # Enable thousands separator
    thousands_sep: bool = False
    # Floating point precision (int/float only)
    digits: int = -1


AnnotationTypes: TypeAlias = Literal["basic", "custom", "builtin", "plugin"]


@dataclass
class Annotation:
    type: AnnotationTypes
    name: str
    rank: int
    fields: List[AnnotationField]
    enabled: Literal[True, False, "mandatory"]
    files: list[str]
    params: List[str]


########################################################################################


def _collect_yaml_files(filepaths: List[Path]) -> List[Path]:
    result: List[Path] = []
    for filepath in filepaths:
        if filepath.is_dir():
            result.extend(filepath.glob("*.yaml"))
            result.extend(filepath.glob("*.yml"))
        else:
            result.append(filepath)

    result.sort(key=lambda it: it.name)
    return result


def _apply_variables(
    values: list[str],
    variables: Dict[str, Union[str, Path]],
) -> list[str]:
    result: list[str] = []
    for value in values:
        last_value = ""
        while value != last_value:
            last_value = value
            value = value.format_map(variables)
        result.append(value)

    return result


def load_annotations(
    filepaths: List[Path],
    variables: Dict[str, Path],
) -> List[Annotation]:
    log = logging.getLogger(__name__)
    yaml = ruamel.yaml.YAML(typ="safe", pure=True)
    yaml.version = (1, 1)

    annotations: Dict[str, Annotation] = {}
    with access_resources("annotations") as built_in:
        for filepath in _collect_yaml_files([built_in, *filepaths]):
            log.debug("reading annotation settings from %s", filepath)
            with filepath.open("rt") as handle:
                data: object = yaml.load(handle)

            if data is None:
                log.warning("file is empty")
                continue

            try:
                result = _Root.model_validate(data, strict=True)
            except ValidationError as error:
                for err in error.errors():
                    raise AnnotationError(
                        "error at {loc}: {msg}: {input!r}".format(
                            loc=".".join(map(str, err["loc"])),
                            input=err["input"],
                            msg=err["msg"],
                        )
                    ) from error

                raise AssertionError("should not happen") from error

            for name, obj in result.root.items():
                annotation = obj.to_annotation(name=name, variables=dict(variables))
                if annotation.name in annotations:
                    log.warning("Overriding settings for annotations %r", name)

                annotations[name] = annotation

    input_paths: set[tuple[str, ...]] = set()
    output_keys: set[str] = set()
    for annotation in annotations.values():
        for field in annotation.fields:
            input_path = tuple(field.input_path)
            if input_path in input_paths:
                log.warning("Multiple values read from %r", input_path)
            if field.output_key in output_keys:
                log.error("Multiple values written to %r", field.output_key)

            input_paths.add(input_path)
            output_keys.add(field.output_key)

    return sorted(annotations.values(), key=lambda it: it.rank)
