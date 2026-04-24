from typing import Tuple

import numpy as np
from ase import Atoms
from ase.data import covalent_radii

FilterResult = Tuple[bool, str]


def check_converged(is_converged: bool) -> FilterResult:
    if not is_converged:
        return False, "not_converged"
    return True, "ok"


def check_min_pair_distance(atoms: Atoms, factor: float = 0.7) -> FilterResult:
    n = len(atoms)
    if n < 2:
        return True, "ok"
    distances = atoms.get_all_distances(mic=True)
    numbers = atoms.get_atomic_numbers()
    radii = np.array([covalent_radii[z] for z in numbers])
    min_allowed = factor * (radii[:, None] + radii[None, :])
    np.fill_diagonal(distances, np.inf)
    violations = distances < min_allowed
    if violations.any():
        i, j = [int(x) for x in np.argwhere(violations)[0]]
        return (
            False,
            f"pair_too_close: {atoms.symbols[i]}-{atoms.symbols[j]} "
            f"at {distances[i, j]:.2f} A",
        )
    return True, "ok"


def check_energy_per_atom(
    energy: float,
    n_atoms: int,
    min_energy_per_atom: float,
    max_above: float = 2.0,
) -> FilterResult:
    epa = energy / n_atoms
    delta = epa - min_energy_per_atom
    if delta > max_above:
        return False, f"energy_too_high: {delta:.3f} eV/atom above batch min"
    return True, "ok"


def check_max_displacement(
    atoms: Atoms,
    reference: Atoms,
    max_disp: float,
) -> FilterResult:
    if len(atoms) != len(reference):
        return False, "atom_count_mismatch_vs_reference"
    cell = atoms.get_cell()
    disps = atoms.get_positions() - reference.get_positions()
    # MIC: wrap fractional displacements to [-0.5, 0.5)
    frac = np.linalg.solve(np.array(cell).T, disps.T).T
    frac -= np.round(frac)
    disps = frac @ np.array(cell)
    max_d = float(np.max(np.linalg.norm(disps, axis=1)))
    if max_d > max_disp:
        return False, f"displacement_too_large: {max_d:.2f} A"
    return True, "ok"
