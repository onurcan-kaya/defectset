from pathlib import Path
from typing import Any, Dict, Optional

from ase import Atoms
from ase.calculators.espresso import Espresso, EspressoProfile


def _deep_update(base: Dict[str, Any], extra: Dict[str, Dict[str, Any]]) -> None:
    for section, vals in extra.items():
        base.setdefault(section, {}).update(vals)


def run_qe_single_point(
    atoms: Atoms,
    pseudopotentials: Dict[str, str],
    pseudo_dir: str,
    command: str = "pw.x",
    ecutwfc: float = 80.0,
    ecutrho: Optional[float] = None,
    kpoints=None,
    kspacing: Optional[float] = 0.04,
    smearing: str = "cold",
    degauss: float = 0.01,
    spin_polarized: bool = False,
    mixing_beta: float = 0.3,
    electron_maxstep: int = 200,
    conv_thr: float = 1.0e-8,
    extra_input_data: Optional[Dict[str, Dict[str, Any]]] = None,
    work_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Run one SCF calculation. Returns a dict with energy, forces, status.

    The launch command is passed explicitly via `command` and turned into an
    EspressoProfile, so no ASE env vars or ~/.config/ase/config.ini are read.
    """

    input_data: Dict[str, Dict[str, Any]] = {
        "control": {
            "calculation": "scf",
            "tstress": False,
            "tprnfor": True,
        },
        "system": {
            "ecutwfc": ecutwfc,
            "occupations": "smearing",
            "smearing": smearing,
            "degauss": degauss,
        },
        "electrons": {
            "conv_thr": conv_thr,
            "mixing_beta": mixing_beta,
            "electron_maxstep": electron_maxstep,
        },
    }

    if ecutrho is not None:
        input_data["system"]["ecutrho"] = ecutrho

    if spin_polarized:
        input_data["system"]["nspin"] = 2
        for i, _sym in enumerate(sorted(set(atoms.get_chemical_symbols())), start=1):
            input_data["system"][f"starting_magnetization({i})"] = 0.1

    if extra_input_data:
        _deep_update(input_data, extra_input_data)

    profile = EspressoProfile(command=command, pseudo_dir=str(pseudo_dir))

    calc_kwargs: Dict[str, Any] = dict(
        profile=profile,
        input_data=input_data,
        pseudopotentials=pseudopotentials,
    )
    if kpoints is not None:
        calc_kwargs["kpts"] = tuple(int(x) for x in kpoints)
    elif kspacing is not None:
        calc_kwargs["kspacing"] = kspacing

    if work_dir is not None:
        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        calc_kwargs["directory"] = str(work_dir)

    atoms_copy = atoms.copy()
    try:
        atoms_copy.calc = Espresso(**calc_kwargs)
        energy = atoms_copy.get_potential_energy()
        forces = atoms_copy.get_forces()
        return {
            "status": "ok",
            "energy": float(energy),
            "forces": forces,
            "atoms": atoms_copy,
        }
    except Exception as e:
        return {
            "status": "failed",
            "error": f"{type(e).__name__}: {e}",
            "atoms": atoms_copy,
        }
