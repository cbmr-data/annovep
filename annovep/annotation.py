from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Union

import pydantic
import ruamel.yaml
from pydantic import ConfigDict, Field, RootModel, ValidationError, model_validator
from typing_extensions import Annotated, Literal, override

from annovep.resources import access_resources

# Built-in annotations derived from input
_BUILTINS = {
    "samplegenotypes": "SampleGenotypes",
}


class AnnotationError(Exception):
    pass


class BaseModel(pydantic.BaseModel):
    model_config = ConfigDict(extra="forbid")


class Column(BaseModel):
    name: str = Field(alias="Name")
    fieldtype: Optional[str] = Field(alias="FieldType", default=None)
    help: Optional[str] = Field(alias="Help", default=None)
    # Normalize lists of items using this separator
    split_by: Optional[str] = Field(alias="SplitBy", default=None)
    # Enable thousands separator
    thousands_sep: Optional[bool] = Field(alias="ThousandsSep", default=None)
    # Floating point precision (int/float only)
    digits: Optional[int] = Field(alias="Digits", default=None)


class AnnotationBase(BaseModel):
    _name: Optional[str] = None
    rank: int = Field(alias="Rank", default=0)
    fieldtype: str = Field(alias="FieldType", default="str")
    thousands_sep: bool = Field(alias="ThousandsSep", default=False)
    digits: int = Field(alias="Digits", default=-1)
    fields: Dict[str, Column] = Field(alias="Fields", default_factory=dict)
    enabled: bool = Field(alias="Enabled", default=True)

    def apply_variables(self, variables: Dict[str, Union[str, Path]]) -> None:
        ...

    @property
    def name(self) -> str:
        if self._name is None:
            raise AssertionError(self)
        return self._name

    def set_name(self, name: str) -> None:
        self._name = name

    @property
    def params(self) -> List[str]:
        return []

    @model_validator(mode="after")
    def set_column_defaults(self) -> AnnotationBase:
        for field in self.fields.values():
            if field.fieldtype is None:
                field.fieldtype = self.fieldtype
            if field.thousands_sep is None:
                field.thousands_sep = self.thousands_sep
            if field.digits is None:
                field.digits = self.digits

        return self


class Option(AnnotationBase):
    type: Literal["Option"] = Field(..., alias="Type")
    command: List[str] = Field(alias="Command")

    @property
    def files(self) -> list[str]:
        return []


class Plugin(AnnotationBase):
    type: Literal["Plugin"] = Field(..., alias="Type")
    files: List[str] = Field(alias="Files")
    parameters: List[str] = Field(alias="Parameters")
    variables: Dict[str, str] = Field(alias="Variables", default_factory=dict)

    @override
    def apply_variables(self, variables: Dict[str, Union[str, Path]]) -> None:
        variables = dict(variables)
        variables.update(self.variables)

        self.files = _apply_variables(self.files, variables)
        self.parameters = _apply_variables(self.parameters, variables)

    @property
    def params(self) -> List[str]:
        params: List[str] = [self.name]
        params.extend(self.parameters)

        return ["--plugin", ",".join(params)]


class Custom(AnnotationBase):
    type: Literal["BED", "VCF"] = Field(..., alias="Type")
    mode: Literal["exact", "overlap"] = Field(..., alias="Mode")
    file: str = Field(..., alias="File")
    variables: Dict[str, str] = Field(alias="Variables", default_factory=dict)

    @property
    def files(self) -> List[str]:
        return [self.file, f"{self.file}.tbi"]

    @override
    def apply_variables(self, variables: Dict[str, Union[str, Path]]) -> None:
        variables = dict(variables)
        variables.update(self.variables)

        (self.file,) = _apply_variables([self.file], variables)

    @property
    def params(self) -> List[str]:
        params = [self.file, self.name, self.type, self.mode, "0"]
        for name in self.fields:
            if not (name.startswith(":") and name.endswith(":")):
                params.append(name)

        return ["--custom", ",".join(params)]


class Builtin(AnnotationBase):
    type: Literal["Builtin"] = Field(..., alias="Type")

    @property
    def files(self) -> list[str]:
        return []

    @override
    def set_name(self, name: str) -> None:
        builtin_name = _BUILTINS.get(name.lower())
        if builtin_name is None:
            raise ValueError(f"unknown built-in annotation {name!r}")

        super().set_name(builtin_name)


Annotations = Annotated[
    Union[Option, Plugin, Custom, Builtin], Field(discriminator="type")
]

AnnotationsDict = Dict[str, Annotations]


class _Root(RootModel[AnnotationsDict]):
    root: AnnotationsDict


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
    last_value: str | None = None
    for value in values:
        while value != last_value:
            last_value = value
            value = value.format_map(variables)
        result.append(value)

    return result


def load_annotations(
    log: logging.Logger,
    filepaths: List[Path],
    variables: Dict[str, Union[str, Path]],
) -> List[Annotations]:
    yaml = ruamel.yaml.YAML(typ="safe", pure=True)
    yaml.version = (1, 1)

    annotations: Dict[str, Annotations] = {}
    with access_resources("annotations") as built_in:
        for filepath in _collect_yaml_files([built_in, *filepaths]):
            log.info("reading annotation settings from %s", filepath)
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
                    )

                raise AssertionError("should not happen")

            for name, obj in result.root.items():
                if name in annotations:
                    log.warning("Overriding settings for annotations %r", name)

                obj.set_name(name)
                obj.apply_variables(variables)

                annotations[name] = obj

    return sorted(annotations.values(), key=lambda it: it.rank)
