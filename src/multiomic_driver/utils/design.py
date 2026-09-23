"""Single validated representation of the experimental design."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ExperimentalDesign:
    """Treatment labels and biological replicates used by the current pipeline."""

    control: str
    treated: str
    replicates: tuple[int, ...]

    @property
    def conditions(self) -> tuple[str, str]:
        return self.control, self.treated

    @classmethod
    def from_config(cls, config: dict) -> ExperimentalDesign:
        conditions = config["conditions"]
        design = cls(
            control=str(conditions["control"]),
            treated=str(conditions["treated"]),
            replicates=tuple(int(value) for value in config["replicates"]),
        )
        if design.control == design.treated or len(set(design.replicates)) != len(
            design.replicates
        ):
            raise ValueError("Experimental conditions and replicates must be unique")
        if len(design.replicates) != 2:
            raise ValueError("The current replicate-aware workflow requires two replicates")
        return design


DEFAULT_DESIGN = ExperimentalDesign("DMSO", "Dasatinib", (1, 2))
