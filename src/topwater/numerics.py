"""Controls and exact dense-operator storage preflight for the Laplace prototype.

This module does not assemble operators or claim numerical verification. The
budget covers S and K only, not geometry, temporary arrays or a factorization.
"""

from dataclasses import asdict, dataclass, field
import json
import math
from numbers import Integral, Real


def _integer(value, name, minimum=1):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _mapping(text):
    try:
        data = json.loads(text)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid numerical record: {error}") from error
    if not isinstance(data, dict):
        raise ValueError("numerical record must be a JSON object")
    if type(data.get("format_version")) is not int or data["format_version"] != 1:
        raise ValueError("unsupported numerical format_version; expected 1")
    return data


@dataclass(frozen=True)
class LaplaceControls:
    """Error requests are on S/panel_diameter and dimensionless K.

    Tolerances are requests, not certified integration-error bounds. Depth caps
    bound fallback work; failure to meet accuracy must raise, not return success.
    Block size is a tunable performance control, not part of the mathematics.
    """

    quadrature_rtol: float = 1e-8
    quadrature_atol: float = 1e-10
    max_refinement: int = 8
    target_block_size: int = 128
    max_operator_bytes: int = 512 * 1024**2

    def __post_init__(self):
        for name in ("quadrature_rtol", "quadrature_atol"):
            value = getattr(self, name)
            if (isinstance(value, bool) or not isinstance(value, Real)
                    or not math.isfinite(value) or not 0 < value < 1):
                raise ValueError(f"{name} must be finite and between 0 and 1")
            object.__setattr__(self, name, float(value))
        for name in ("target_block_size", "max_operator_bytes", "max_refinement"):
            object.__setattr__(self, name, _integer(
                getattr(self, name), name, minimum=0 if name == "max_refinement" else 1,
            ))

    def to_json(self):
        return json.dumps({"format_version": 1, **asdict(self)}, allow_nan=False, indent=2)

    @classmethod
    def from_json(cls, text):
        data = _mapping(text)
        data.pop("format_version")
        try:
            return cls(**data)
        except TypeError as error:
            raise ValueError(f"Invalid Laplace controls: {error}") from error


@dataclass(frozen=True)
class DenseOperatorPlan:
    """Preflight two square float64 matrices before numerical allocation.

    This is a computed storage record, not an assembly-result or timing record.
    A successful preflight does not guarantee sufficient total process memory.
    """

    panel_count: int
    controls: LaplaceControls = field(default_factory=LaplaceControls)

    def __post_init__(self):
        object.__setattr__(self, "panel_count", _integer(self.panel_count, "panel_count"))
        if not isinstance(self.controls, LaplaceControls):
            raise ValueError("controls must be LaplaceControls")
        if self.operator_bytes > self.controls.max_operator_bytes:
            raise MemoryError(
                f"S and K need {self.operator_bytes} bytes, exceeding max_operator_bytes="
                f"{self.controls.max_operator_bytes}; this excludes solve/scratch memory"
            )

    @property
    def operator_bytes(self):
        return 2 * 8 * self.panel_count**2

    @property
    def block_rows(self):
        return min(self.controls.target_block_size, self.panel_count)

    def to_json(self):
        return json.dumps({
            "format_version": 1, "panel_count": self.panel_count,
            "operator_bytes": self.operator_bytes, "block_rows": self.block_rows,
            "controls": asdict(self.controls),
        }, allow_nan=False, indent=2)

    @classmethod
    def from_json(cls, text):
        data = _mapping(text)
        if set(data) != {"format_version", "panel_count", "operator_bytes", "block_rows", "controls"}:
            raise ValueError("unexpected or missing dense-operator plan fields")
        try:
            result = cls(data["panel_count"], LaplaceControls(**data["controls"]))
        except TypeError as error:
            raise ValueError(f"Invalid dense-operator plan: {error}") from error
        for name in ("operator_bytes", "block_rows"):
            if type(data[name]) is not int or data[name] != getattr(result, name):
                raise ValueError(f"{name} does not match the computed plan")
        return result
