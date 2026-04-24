from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import yaml


@dataclass
class FilterConfig:
    # Fraction of (r_cov_i + r_cov_j). Structures with any pair closer than
    # this are rejected as MLIP collapses.
    min_pair_distance_factor: float = 0.7

    # eV/atom above the lowest-energy converged structure in the batch.
    max_energy_per_atom_above_min: float = 2.0

    # Angstrom. Largest per-atom displacement from the `undistorted` reference
    # in the defectool DB. Set to null/None to skip.
    max_displacement: Optional[float] = 3.0

    # Require `is_converged=True` in the defectool DB.
    require_converged: bool = True


@dataclass
class QEConfig:
    pseudopotentials: Dict[str, str]
    pseudo_dir: str

    # Launch command for pw.x. Examples: "pw.x", "mpirun -np 8 pw.x",
    # "srun pw.x". Parsed with shlex; no shell metacharacters expanded.
    command: str = "pw.x"

    ecutwfc: float = 80.0
    ecutrho: Optional[float] = None
    kpoints: Optional[List[int]] = None
    kspacing: Optional[float] = 0.04
    smearing: str = "cold"
    degauss: float = 0.01
    spin_polarized: bool = False
    mixing_beta: float = 0.3
    electron_maxstep: int = 200
    conv_thr: float = 1.0e-8

    # Escape hatch. Deep-merged into the QE input_data namelists.
    # Use e.g. {"system": {"tot_charge": -1}} to do a charged calculation.
    extra_input_data: Dict[str, Dict[str, Any]] = field(default_factory=dict)


@dataclass
class OutputConfig:
    dataset_path: str = "defect_dataset.extxyz"
    rejected_path: str = "defect_rejected.extxyz"
    log_path: str = "defectset.log"


@dataclass
class Config:
    defectool_config: str
    qe: QEConfig
    filters: FilterConfig = field(default_factory=FilterConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    skip_defectool: bool = False

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "Config":
        with open(path) as f:
            data = yaml.safe_load(f)
        if not isinstance(data, dict):
            raise ValueError(f"Config root must be a mapping, got {type(data).__name__}")

        qe_data = data.pop("qe", None)
        if qe_data is None:
            raise ValueError("Config must contain a `qe` section.")
        filters_data = data.pop("filters", {}) or {}
        output_data = data.pop("output", {}) or {}

        return cls(
            qe=QEConfig(**qe_data),
            filters=FilterConfig(**filters_data),
            output=OutputConfig(**output_data),
            **data,
        )

    def defectool_output_dir(self) -> Path:
        with open(self.defectool_config) as f:
            d = yaml.safe_load(f) or {}
        return Path(d.get("output_dir", "defectool_output"))
