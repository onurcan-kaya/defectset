# defectset

Close the loop from [`defectool`](https://github.com/onurcan-kaya/defectool)-generated defect structures to a DFT-labelled training set.

`defectset` runs `defectool` (MLIP-based distort → rattle → relax), filters out the broken geometries the MLIP produced, runs a single-point Quantum ESPRESSO SCF on each survivor, and writes an `extxyz` dataset with DFT energies and forces. The MLIP pre-filter is the point: structures that an MLIP collapses, fails to converge on, or relaxes into some bizarre basin are exactly the ones that would waste hours of QE wall-time and produce low-quality labels.

## What it does

```
defectool yaml
      │
      ├── defectool generate   (build defect supercells + distortions + rattles)
      ├── defectool relax      (MLIP relaxation, default MACE-MP-0 via ASE)
      │
      ▼
defectool.db  ──►  filter  ──►  QE single-point  ──►  extxyz with REF_energy/REF_forces
                     │
                     └── rejected.extxyz (reason logged per structure)
```

Four filters run before QE:

| Filter | Rejects |
|---|---|
| `require_converged` | structures defectool marked `is_converged=False` (force divergence, >max_steps, etc.) |
| `min_pair_distance_factor` | pairs closer than `factor × (r_cov_i + r_cov_j)` — catches MLIP atom collapses |
| `max_energy_per_atom_above_min` | structures with MLIP energy per atom more than X eV above the batch minimum |
| `max_displacement` | atoms displaced more than X Å (MIC) from the converged `undistorted` reference |

All thresholds are configurable. Each reject is written to `defect_rejected.extxyz` with a `reject_reason` in the info dict so you can audit.

QE calculations are cached per label in `<work_dir>/qe_cache/`. Re-running the pipeline reuses successful QE outputs; failed ones are re-attempted.

## Prerequisites

- **Python ≥ 3.10**
- **ASE ≥ 3.23** (installed automatically as a dependency). `defectset` uses ASE's `EspressoProfile` API, which requires ASE 3.23+. It does **not** use the legacy `ASE_ESPRESSO_COMMAND` env var or `~/.config/ase/config.ini`.
- **[`defectool`](https://github.com/onurcan-kaya/defectool)** on your `PATH`. `defectset` invokes it via the `defectool` CLI. Install with:
  ```bash
  git clone https://github.com/onurcan-kaya/defectool.git
  cd defectool && pip install -e ".[mace]"
  ```
  Install defectool and defectset into the **same** Python environment, otherwise the `defectool` command will not be found when `defectset` shells out.
- **An MLIP calculator usable by defectool.** The default is MACE-MP-0, pulled in by the `[mace]` extra above. Any defectool-supported calculator works — the pre-filter quality just depends on how well it relaxes your defects.
- **Quantum ESPRESSO `pw.x`** on your `PATH`, and pseudopotentials for every element in your system. The path to pseudos goes into `qe.pseudo_dir` in the defectset config; the launch command goes into `qe.command`.

`defectset` is intentionally local-only: no scheduler integration, no job arrays. One `pw.x` SCF per structure, run in-process. Wrap it in your own batch script if you need parallel labelling.

## Install

```bash
git clone <this-repo> defectset
cd defectset
pip install -e .
```

`defectool` is **not** declared as a Python dependency because it is not on PyPI. Install it yourself as above, and make sure `which defectool` resolves.

## Usage

### CLI

```bash
defectset run config.yaml
defectset run config.yaml --work-dir ./my_run
defectset run config.yaml --skip-defectool    # reuse existing defectool.db
```

### Python

```python
from pathlib import Path
from defectset import Config, run_pipeline

cfg = Config.from_yaml("config.yaml")
counts = run_pipeline(cfg, work_dir=Path("my_run"))
# counts -> {"success": N, "qe_failed": M, "filtered": K}
```

## Config reference

The config is a YAML file with four sections: `defectool_config`, `filters`, `qe`, `output`. Only `defectool_config` and `qe` are required. See [`example_config.yaml`](./example_config.yaml) for an annotated template.

### `defectool_config`
Path to an existing defectool YAML. `defectset` reads its `output_dir` to find the resulting `defectool.db`. Defaults match defectool (`defectool_output/`).

### `filters`

| Key | Default | Meaning |
|---|---|---|
| `min_pair_distance_factor` | `0.7` | Reject if any pair distance < `factor × (r_cov_i + r_cov_j)` |
| `max_energy_per_atom_above_min` | `2.0` | eV/atom above batch minimum MLIP energy |
| `max_displacement` | `3.0` | Å from converged `undistorted` reference. `null` to skip |
| `require_converged` | `true` | Require `is_converged=True` in defectool DB |

### `qe`

| Key | Default | Meaning |
|---|---|---|
| `pseudopotentials` | **required** | Mapping `element: UPF filename` |
| `pseudo_dir` | **required** | Directory containing UPFs. Use an absolute path |
| `command` | `"pw.x"` | Full launch command. e.g. `"pw.x"`, `"mpirun -np 8 pw.x"`, `"srun pw.x"`. Parsed with `shlex` — no shell metacharacters expanded |
| `ecutwfc` | `80.0` | Ry |
| `ecutrho` | `null` (= 4·ecutwfc) | Ry |
| `kpoints` | `null` | Explicit `[nx, ny, nz]` (overrides kspacing) |
| `kspacing` | `0.04` | Å⁻¹ |
| `smearing` | `cold` | QE smearing type |
| `degauss` | `0.01` | Ry |
| `spin_polarized` | `false` | Set `nspin=2` + small starting magnetisation |
| `mixing_beta` | `0.3` | |
| `electron_maxstep` | `200` | |
| `conv_thr` | `1.0e-8` | |
| `extra_input_data` | `{}` | Deep-merged into QE namelists. Escape hatch for anything above |

### `output`

| Key | Default | Meaning |
|---|---|---|
| `dataset_path` | `defect_dataset.extxyz` | Final DFT-labelled dataset |
| `rejected_path` | `defect_rejected.extxyz` | All rejects (filtered + QE-failed) with `reject_reason` |
| `log_path` | `defectset.log` | |

All output paths are resolved relative to `--work-dir`.

## Output format

`defect_dataset.extxyz` is extended XYZ. Each frame:

- `info["REF_energy"]` — DFT total energy (eV)
- `arrays["REF_forces"]` — DFT forces (eV/Å, shape N×3)
- `info["label"]` — defectool label (e.g. `bond_-0.02`, `rattle_3`, `undistorted`)
- `info["defect_type"]`, `info["distortion_type"]`, `info["distortion_mag"]` — passthrough from defectool
- `info["mlip_energy"]` — the MLIP total energy that defectool produced, for comparison

Compatible out of the box with MLIP training frameworks that auto-detect `REF_energy`/`REF_forces` keys (MACE, NequIP/Allegro, AutoMLIP).

## Working directory layout

```
<work_dir>/
    defectset.log
    defect_dataset.extxyz         # DFT-labelled survivors
    defect_rejected.extxyz        # filtered + QE-failed, with reject_reason
    qe/<label>/                   # raw pw.x working dirs (pwi, pwo, tmp/)
    qe_cache/<label>.extxyz       # per-label successful QE result (used for restart)
```

Delete `qe_cache/` to force re-labelling. Delete `qe/` to save space after the dataset is built.

## Known limitations

- **Neutral cells only by default.** `spin_polarized: true` works; for charged defects you need to set `tot_charge` via `extra_input_data.system` *and* accept that a naive plane-wave charged-cell energy needs a finite-size correction (FNV / eFNV) that this tool does not do. Use `doped` for that step.
- **No parallel QE scheduling.** One SCF at a time, in-process. For >~50 structures, wrap `defectset run` calls in your own job script, or split the defectool output into chunks.
- **Reference-based displacement filter assumes the `undistorted` label exists.** If you ran defectool with only distortions or rattles and no undistorted reference, this filter is silently skipped (warning in log). Set `max_displacement: null` to disable explicitly.
- **MLIP pre-filter quality is only as good as the MLIP.** MACE-MP-0 is a reasonable default for main-group chemistries but misbehaves on charged, polaronic, or strongly correlated systems. Inspect `defect_rejected.extxyz` to audit.

## License

MIT.
