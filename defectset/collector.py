from pathlib import Path
from typing import Any, Dict, List, Tuple

from ase import Atoms
from ase.io import write


def _serializable_info(meta: Dict[str, Any]) -> Dict[str, Any]:
    """extxyz writer can't handle None; drop those keys."""
    out = {}
    for k, v in meta.items():
        if v is None:
            continue
        out[k] = v
    return out


def write_datasets(
    qe_results: List[Tuple[Atoms, Dict[str, Any], Dict[str, Any]]],
    rejected: List[Tuple[Atoms, Dict[str, Any], str]],
    dataset_path: Path,
    rejected_path: Path,
) -> Dict[str, int]:
    """Write survivors with DFT labels (REF_energy, REF_forces) and a rejected log.

    Returns a counts dict: {"success": n, "qe_failed": m, "filtered": k}.
    """

    dataset_path = Path(dataset_path)
    rejected_path = Path(rejected_path)
    dataset_path.parent.mkdir(parents=True, exist_ok=True)
    rejected_path.parent.mkdir(parents=True, exist_ok=True)

    success_atoms: List[Atoms] = []
    rejected_atoms: List[Atoms] = []

    qe_failed = 0
    for atoms, meta, result in qe_results:
        if result["status"] == "ok":
            a = result["atoms"].copy()
            a.info["REF_energy"] = float(result["energy"])
            a.arrays["REF_forces"] = result["forces"]
            info = _serializable_info({
                "label": meta.get("label", ""),
                "defect_type": meta.get("defect_type", ""),
                "distortion_type": meta.get("distortion_type", ""),
                "distortion_mag": meta.get("distortion_mag"),
                "mlip_energy": meta.get("total_energy"),
            })
            a.info.update(info)
            success_atoms.append(a)
        else:
            qe_failed += 1
            a = atoms.copy()
            info = _serializable_info({
                "label": meta.get("label", ""),
                "reject_reason": f"qe_failed: {result.get('error', 'unknown')}",
                "mlip_energy": meta.get("total_energy"),
                "defect_type": meta.get("defect_type", ""),
                "distortion_type": meta.get("distortion_type", ""),
                "distortion_mag": meta.get("distortion_mag"),
            })
            a.info.update(info)
            rejected_atoms.append(a)

    for atoms, meta, reason in rejected:
        a = atoms.copy()
        info = _serializable_info({
            "label": meta.get("label", ""),
            "reject_reason": reason,
            "mlip_energy": meta.get("total_energy"),
            "defect_type": meta.get("defect_type", ""),
            "distortion_type": meta.get("distortion_type", ""),
            "distortion_mag": meta.get("distortion_mag"),
        })
        a.info.update(info)
        rejected_atoms.append(a)

    if success_atoms:
        write(str(dataset_path), success_atoms, format="extxyz")
    if rejected_atoms:
        write(str(rejected_path), rejected_atoms, format="extxyz")

    return {
        "success": len(success_atoms),
        "qe_failed": qe_failed,
        "filtered": len(rejected),
    }
