import logging
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from ase import Atoms
from ase.db import connect
from ase.io import read, write

from .collector import write_datasets
from .config import Config, FilterConfig, QEConfig
from .filters import (
    check_converged,
    check_energy_per_atom,
    check_max_displacement,
    check_min_pair_distance,
)
from .qe import run_qe_single_point

logger = logging.getLogger(__name__)


def run_defectool(defectool_config: str) -> None:
    for stage in ("generate", "relax"):
        logger.info("Running defectool %s on %s", stage, defectool_config)
        result = subprocess.run(
            ["defectool", stage, defectool_config],
            capture_output=True,
            text=True,
        )
        if result.stdout:
            logger.debug("defectool %s stdout:\n%s", stage, result.stdout)
        if result.returncode != 0:
            logger.error("defectool %s failed:\n%s", stage, result.stderr)
            raise RuntimeError(f"defectool {stage} failed (exit {result.returncode})")


def load_relaxed_from_db(db_path: Path) -> List[Tuple[Atoms, Dict[str, Any]]]:
    db = connect(str(db_path))
    rows: List[Tuple[Atoms, Dict[str, Any]]] = []
    for row in db.select():
        atoms = row.toatoms()
        meta = {
            "label": str(getattr(row, "label", f"id{row.id}")),
            "is_converged": bool(getattr(row, "is_converged", False)),
            "total_energy": _maybe_float(getattr(row, "total_energy", None)),
            "final_fmax": _maybe_float(getattr(row, "final_fmax", None)),
            "defect_type": str(getattr(row, "defect_type", "") or ""),
            "distortion_type": str(getattr(row, "distortion_type", "") or ""),
            "distortion_mag": _maybe_float(getattr(row, "distortion_mag", None)),
        }
        rows.append((atoms, meta))
    return rows


def _maybe_float(x) -> Optional[float]:
    if x is None:
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _find_reference(
    rows: List[Tuple[Atoms, Dict[str, Any]]],
) -> Optional[Atoms]:
    for atoms, meta in rows:
        if meta.get("label") == "undistorted" and meta.get("is_converged"):
            return atoms
    return None


def apply_filters(
    rows: List[Tuple[Atoms, Dict[str, Any]]],
    filter_config: FilterConfig,
) -> Tuple[
    List[Tuple[Atoms, Dict[str, Any]]],
    List[Tuple[Atoms, Dict[str, Any], str]],
]:
    survivors: List[Tuple[Atoms, Dict[str, Any]]] = []
    rejected: List[Tuple[Atoms, Dict[str, Any], str]] = []
    reference = _find_reference(rows)

    # Stage 1: per-structure filters (convergence, pair distance).
    stage1: List[Tuple[Atoms, Dict[str, Any]]] = []
    for atoms, meta in rows:
        if filter_config.require_converged:
            ok, reason = check_converged(meta["is_converged"])
            if not ok:
                rejected.append((atoms, meta, reason))
                continue

        ok, reason = check_min_pair_distance(atoms, filter_config.min_pair_distance_factor)
        if not ok:
            rejected.append((atoms, meta, reason))
            continue

        stage1.append((atoms, meta))

    # Stage 2: batch-relative filters (energy vs min, displacement vs reference).
    min_epa: Optional[float] = None
    epas = [
        (m["total_energy"] / len(a))
        for a, m in stage1
        if m["total_energy"] is not None and len(a) > 0
    ]
    if epas:
        min_epa = min(epas)

    for atoms, meta in stage1:
        if min_epa is not None and meta["total_energy"] is not None:
            ok, reason = check_energy_per_atom(
                meta["total_energy"],
                len(atoms),
                min_epa,
                filter_config.max_energy_per_atom_above_min,
            )
            if not ok:
                rejected.append((atoms, meta, reason))
                continue

        if filter_config.max_displacement is not None and reference is not None:
            ok, reason = check_max_displacement(
                atoms, reference, filter_config.max_displacement
            )
            if not ok:
                rejected.append((atoms, meta, reason))
                continue

        survivors.append((atoms, meta))

    if filter_config.max_displacement is not None and reference is None:
        logger.warning(
            "max_displacement filter requested but no converged `undistorted` "
            "reference in DB; displacement check skipped."
        )

    return survivors, rejected


def _cache_get(cache_dir: Path, label: str) -> Optional[Dict[str, Any]]:
    path = cache_dir / f"{label}.extxyz"
    if not path.exists():
        return None
    try:
        a = read(str(path))
        if "REF_energy" in a.info and "REF_forces" in a.arrays:
            return {
                "status": "ok",
                "energy": float(a.info["REF_energy"]),
                "forces": a.arrays["REF_forces"],
                "atoms": a,
            }
    except Exception as e:
        logger.warning("QE cache read failed for %s: %s", label, e)
    return None


def _cache_put(cache_dir: Path, label: str, result: Dict[str, Any]) -> None:
    if result.get("status") != "ok":
        return
    cache_dir.mkdir(parents=True, exist_ok=True)
    a = result["atoms"].copy()
    a.info["REF_energy"] = float(result["energy"])
    a.arrays["REF_forces"] = result["forces"]
    write(str(cache_dir / f"{label}.extxyz"), a, format="extxyz")


def label_with_qe(
    survivors: List[Tuple[Atoms, Dict[str, Any]]],
    qe_config: QEConfig,
    work_dir: Path,
) -> List[Tuple[Atoms, Dict[str, Any], Dict[str, Any]]]:
    cache_dir = work_dir / "qe_cache"
    results: List[Tuple[Atoms, Dict[str, Any], Dict[str, Any]]] = []
    n = len(survivors)
    for i, (atoms, meta) in enumerate(survivors, start=1):
        label = meta["label"]
        cached = _cache_get(cache_dir, label)
        if cached is not None:
            logger.info("[%d/%d] %s: cached", i, n, label)
            results.append((atoms, meta, cached))
            continue

        job_dir = work_dir / "qe" / label
        logger.info("[%d/%d] %s: running QE", i, n, label)
        result = run_qe_single_point(
            atoms,
            pseudopotentials=qe_config.pseudopotentials,
            pseudo_dir=qe_config.pseudo_dir,
            command=qe_config.command,
            ecutwfc=qe_config.ecutwfc,
            ecutrho=qe_config.ecutrho,
            kpoints=qe_config.kpoints,
            kspacing=qe_config.kspacing,
            smearing=qe_config.smearing,
            degauss=qe_config.degauss,
            spin_polarized=qe_config.spin_polarized,
            mixing_beta=qe_config.mixing_beta,
            electron_maxstep=qe_config.electron_maxstep,
            conv_thr=qe_config.conv_thr,
            extra_input_data=qe_config.extra_input_data,
            work_dir=job_dir,
        )
        if result["status"] == "ok":
            _cache_put(cache_dir, label, result)
            logger.info("[%d/%d] %s: E = %.6f eV", i, n, label, result["energy"])
        else:
            logger.warning("[%d/%d] %s: QE FAILED (%s)", i, n, label, result.get("error"))
        results.append((atoms, meta, result))
    return results


def _setup_logging(log_file: Path) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    # Only configure once.
    if getattr(root, "_defectset_configured", False):
        return
    root.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    fh = logging.FileHandler(log_file)
    fh.setFormatter(fmt)
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    root.addHandler(fh)
    root.addHandler(sh)
    root._defectset_configured = True  # type: ignore[attr-defined]


def run_pipeline(config: Config, work_dir: Path = Path("defectset_work")) -> Dict[str, int]:
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    _setup_logging(work_dir / config.output.log_path)

    if not config.skip_defectool:
        run_defectool(config.defectool_config)
    else:
        logger.info("skip_defectool=True -> reading existing DB only")

    db_path = config.defectool_output_dir() / "defectool.db"
    if not db_path.exists():
        raise FileNotFoundError(f"defectool database not found at {db_path}")

    rows = load_relaxed_from_db(db_path)
    logger.info("Loaded %d structures from %s", len(rows), db_path)

    survivors, rejected = apply_filters(rows, config.filters)
    logger.info(
        "After filters: %d survivors, %d rejected", len(survivors), len(rejected)
    )
    for atoms, meta, reason in rejected:
        logger.info("  rejected %s: %s", meta.get("label", "?"), reason)

    qe_results = label_with_qe(survivors, config.qe, work_dir)

    counts = write_datasets(
        qe_results,
        rejected,
        dataset_path=work_dir / config.output.dataset_path,
        rejected_path=work_dir / config.output.rejected_path,
    )
    logger.info(
        "Done. DFT-labelled: %d, QE-failed: %d, pre-filtered: %d",
        counts["success"],
        counts["qe_failed"],
        counts["filtered"],
    )
    return counts
