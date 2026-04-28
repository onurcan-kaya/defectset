# defectset

Takes defect structures from [`defectool`](https://github.com/onurcan-kaya/defectool), filters out the junk, runs Quantum ESPRESSO single-point SCFs on what survives, and writes an extxyz dataset with DFT energies and forces. You get a training set for MLIPs without wasting QE time on structures the MLIP mangled.

## Pipeline

```
defectool generate + relax
            |
            v
      defectool.db
            |
      [stage 1 filters]  -- per-structure: convergence, pair distance
            |
      [stage 2 filters]  -- batch-relative: energy/atom, displacement
            |
      QE single-point SCF (cached per label)
            |
    +-------+-------+
    |               |
    v               v
defect_dataset    defect_rejected
  .extxyz           .extxyz
```

`defectset` shells out to `defectool generate` then `defectool relax` (unless you pass `--skip-defectool`), reads the resulting ASE database, applies filters in two stages, runs `pw.x` on each survivor, and writes the output.

## Filters

Filters run in order. A structure rejected in stage 1 never reaches stage 2.

**Stage 1 -- per-structure:**

- `require_converged` (default: true) -- rejects anything defectool marked `is_converged=False`.
- `min_pair_distance_factor` (default: 0.7) -- rejects if any pair distance is less than `factor * (r_cov_i + r_cov_j)`, using ASE covalent radii and MIC distances. Catches MLIP atom collapses.

**Stage 2 -- batch-relative (runs on stage 1 survivors only):**

- `max_energy_per_atom_above_min` (default: 2.0 eV/atom) -- computes MLIP energy/atom for all stage 1 survivors, finds the minimum, rejects anything more than X eV/atom above it.
- `max_displacement` (default: 3.0 Ang) -- computes MIC displacement of each atom from the converged `undistorted` structure in the DB. Rejects if any atom moved more than X Ang. If there is no converged `undistorted` entry in the DB, this filter is skipped (warning in log). Set to `null` to disable.

Every reject is written to `defect_rejected.extxyz` with a `reject_reason` string in the info dict.

## QE caching

QE results are cached in `<work_dir>/qe_cache/<label>.extxyz`. Re-running the pipeline skips any label that already has a cached result. Failed QE runs are not cached, so they get retried.

**Watch out:** the cache key is the defectool `label` string. If your defectool config produces duplicate labels, the second structure with the same label will read the first one's cached result. This is a real problem if you have multiple defect types producing identically-named distortions.

## Prerequisites

- Python >= 3.10
- ASE >= 3.23 (installed as a dependency). The code uses `EspressoProfile` directly -- it does not use `ASE_ESPRESSO_COMMAND` or `~/.config/ase/config.ini`.
- [`defectool`](https://github.com/onurcan-kaya/defectool) on your PATH. Not on PyPI, install it yourself:

```
git clone https://github.com/onurcan-kaya/defectool.git
cd defectool && pip install -e ".[mace]"
```

Must be in the same Python environment as defectset, otherwise the subprocess call to `defectool` will not find it.

- Quantum ESPRESSO `pw.x` on your PATH.
- Pseudopotentials for every element in your system.

No scheduler integration. One `pw.x` call at a time, in-process. Wrap it yourself if you need parallelism.

## Install

```
git clone https://github.com/onurcan-kaya/defectset.git
cd defectset
pip install -e .
```

## Usage

### CLI

```
defectset run config.yaml                       # work dir defaults to ./defectset_work
defectset run config.yaml --work-dir ./my_run   # custom work dir
defectset run config.yaml --skip-defectool      # reuse existing defectool.db
```

### Python

```python
from pathlib import Path
from defectset import Config, run_pipeline

cfg = Config.from_yaml("config.yaml")
counts = run_pipeline(cfg, work_dir=Path("my_run"))
# counts = {"success": N, "qe_failed": M, "filtered": K}
```

## Config

YAML with four sections. Only `defectool_config` and `qe` are required. Full annotated example in [`example_config.yaml`](example_config.yaml).

### `defectool_config`

Path to an existing defectool YAML. `defectset` reads that file's `output_dir` field to locate `defectool.db`. If the defectool YAML has no `output_dir`, it defaults to `defectool_output/`.

### `filters`

| Key | Default | What it does |
| --- | --- | --- |
| `require_converged` | `true` | Reject `is_converged=False` |
| `min_pair_distance_factor` | `0.7` | Reject if any pair < `factor * (r_cov_i + r_cov_j)` |
| `max_energy_per_atom_above_min` | `2.0` | eV/atom above batch min MLIP energy |
| `max_displacement` | `3.0` | Ang from converged `undistorted` ref. `null` to skip |

### `qe`

| Key | Default | What it does |
| --- | --- | --- |
| `pseudopotentials` | required | Map of `element: UPF_filename` |
| `pseudo_dir` | required | Absolute path to UPF directory |
| `command` | `"pw.x"` | Launch command, parsed with `shlex`. Example: `"mpirun -np 8 pw.x"` |
| `ecutwfc` | `80.0` | Ry |
| `ecutrho` | `null` (= 4 * ecutwfc) | Ry |
| `kpoints` | `null` | Explicit `[nx, ny, nz]`. Overrides kspacing |
| `kspacing` | `0.04` | Reciprocal-space spacing in 1/Ang |
| `smearing` | `cold` | QE smearing type |
| `degauss` | `0.01` | Ry |
| `spin_polarized` | `false` | Sets `nspin=2` and `starting_magnetization(i)=0.1` for each unique element |
| `mixing_beta` | `0.3` | |
| `electron_maxstep` | `200` | |
| `conv_thr` | `1.0e-8` | |
| `extra_input_data` | `{}` | Dict of dicts, deep-merged into QE namelists. E.g. `{system: {tot_charge: -1}}` |

### `output`

| Key | Default |
| --- | --- |
| `dataset_path` | `defect_dataset.extxyz` |
| `rejected_path` | `defect_rejected.extxyz` |
| `log_path` | `defectset.log` |

All paths resolved relative to `--work-dir`.

## Output format

Each frame in `defect_dataset.extxyz`:

- `info["REF_energy"]` -- DFT total energy in eV
- `arrays["REF_forces"]` -- DFT forces in eV/Ang, shape (N, 3)
- `info["label"]` -- defectool label, e.g. `bond_-0.02`, `rattle_3`, `undistorted`
- `info["defect_type"]`, `info["distortion_type"]`, `info["distortion_mag"]` -- from defectool
- `info["mlip_energy"]` -- MLIP total energy from defectool

Works with MACE, NequIP/Allegro, and anything else that reads `REF_energy`/`REF_forces` keys.

Each frame in `defect_rejected.extxyz` has the same metadata plus `info["reject_reason"]`.

## Work directory layout

```
<work_dir>/
    defectset.log
    defect_dataset.extxyz
    defect_rejected.extxyz
    qe/<label>/                   # pw.x working dirs (pwi, pwo, tmp/)
    qe_cache/<label>.extxyz       # per-label cached QE result
```

Delete `qe_cache/` to force re-labelling. Delete `qe/` to reclaim disk after the dataset is built.

## License

MIT.
