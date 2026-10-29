# counterfeit — Ground Factory Pipeline

Iteratively optimize a crystal structure model to near-perfect geometry while
maintaining a reasonable fit to X-ray data.

**Strategy:** Run `phenix.refine` with a heavily down-weighted X-ray term
(`wxc_scale=0.1`) so geometry restraints dominate.  Each cycle: validate with
MolProbity, update harmonic hold-back restraints on atoms that are getting
worse, then add or tighten one geometry-correction restraint per outlier type.
Repeat until the weighted energy score converges (target wE < 9).

---

## Pipeline scripts

### `ground_factory_runme.com`

Main orchestrator.  Runs the iterative phenix loop and calls all sub-scripts.

```
ground_factory_runme.com  model.pdb  data.mtz  [options]

  maxcycles=N        max iterations (default 300)
  maxadd=N           new restraints added per outlier type per cycle (default 1)
  wxc=0.1            X-ray coordinate weight scale
  wxu=0.05           X-ray ADP weight scale
  prefix=start       phenix output prefix  (start_001.pdb, start_002.pdb …)
  srun=1             submit phenix jobs via "srun -p debug" (1 = yes)
  propagate=0        if 1, each iteration refines from previous output
  snap_interval=0    if >0, run a full-weight snap iteration every N cycles
  snap_wxc=1.0       wxc_scale for snap iterations
  fudgeify=0         if 1, use fudgeify.py instead of reweight_opts for
                     phenix_fudge.eff (recommended for new runs)
  overrides=file     TORSION OVERRIDE file — residues with genuine outliers
                     that must not be tightened
```

**Outputs**

| File | Description |
|------|-------------|
| `start_NNN.pdb` | phenix model each iteration |
| `itr_vs_energy.txt` | weighted energy vs iteration number |
| `phenix_holdme.eff` | harmonic positional restraints (good atoms held back) |
| `phenix_fudge.eff` | geometry-correction restraints (problem atoms) |
| `refmac_opts.txt` | same restraints in REFMAC external format |
| `fudge_table.txt` | fudgeify state: sigma history across cycles |
| `fmodelme.pdb` | final model after REFMAC B-factor convergence |
| `nonoise.mtz` / `simplenoise.mtz` | ideal / ~3%-noisy structure factors |

**Loop per iteration**

1. `phenix.refine` with `wxc_scale=0.1`
2. `molprobify_runme.com` → weighted energy, fullgeo file
3. `redo_harmonics_runme.com` → update `phenix_holdme.eff`
4. `reweight_opts_runme.com` → update `refmac_opts.txt`; writes `phenix_fudge.eff`
   (overwritten by fudgeify when `fudgeify=1`)
5. `fudgeify.py` (if `fudgeify=1`) → update `fudge_table.txt`, overwrite `phenix_fudge.eff`
6. 1 % sigma escalation on all entries in `refmac_opts.txt`

---

### `molprobify_runme.com`

Runs the full MolProbity validation suite and computes the weighted energy score:

```
energy = Σ (deviation/sigma)²
```

weighted by a Gaussian-population correction so that rare large deviates count
more than routine ones.  Writes a `*_fullgeo.txt` file — one line per outlier
in the format:

```
TYPE  energy  delta  obs  ideal  sigma | atom1 altid resname chain resi - atom2 …
```

Types: `NONBOND`, `BOND`, `ANGLE`, `TORSION`, `PLANE`.

---

### `redo_harmonics_runme.com`

Compares the starting model to the current iteration output.  Atoms that became
worse (appeared in bad-atom lists for the first time) get harmonic positional
restraints (`sigma = 0.1 Å`) in `phenix_holdme.eff` so the X-ray term cannot
pull them further from ideal geometry.  Waters that drifted > 0.1 Å are also
restrained.

---

### `reweight_opts_runme.com`

Reads `refmac_opts_debump.txt` / `refmac_opts_fix*.txt` (REFMAC external
restraint proposals from molprobify) and builds the working `refmac_opts.txt`.

Sigma-tightening logic:
- Each restraint that still appears as an outlier is tightened by ÷1.5.
- Each restraint that has resolved is escalated by ×1.01 (and eventually dropped
  once sigma ≥ 5 for BOND/NONBOND).
- Sigma floors prevent runaway tightening:

| Type | REFMAC floor | phenix_fudge.eff floor |
|------|-------------|----------------------|
| NONBOND (interval) | 0.1 | — |
| BOND (dist) | 0.001 Å | 0.001 Å |
| ANGLE | 0.01 ° | 0.01 ° |
| TORSION | 0.01 ° | 0.01 ° |
| PLANE | **0.001** (fixed from 1e-7) | 0.001 |
| CBdev | 2.0 ° | 2.0 ° |

Also assembles `phenix_fudge.eff` from individual `phenix_opts_*.txt` files.
When `fudgeify=1`, this file is immediately overwritten by `fudgeify.py`.

---

### `fudgeify.py`

Python replacement for the tcsh/awk fudge-generation logic in
`reweight_opts_runme.com`.  Recommended for new runs (`fudgeify=1`).

Reads the current iteration's `*_fullgeo.txt`; accumulates sigma state in
`fudge_table.txt`; writes `phenix_fudge.eff`.

**Sigma behaviour**

| Event | Action |
|-------|--------|
| Restraint appears in new outlier list | Add at natural sigma (NONBOND: `0.4/(|delta|+1.1)`) |
| Restraint still in outlier list | Tighten: `sigma /= 1.5` (with type floor) |
| Restraint no longer an outlier | Escalate: `sigma *= 1.01` |
| NONBOND or BOND sigma ≥ 5 | Drop (restraint resolved) |

**Sigma floors on tightening** (prevents numerical blowup in phenix)

| Type | Floor |
|------|-------|
| BOND | 0.001 Å |
| NONBOND | 0.001 Å |
| ANGLE | **0.5 °** |
| TORSION | **5.0 °** |
| PLANE | 0.001 (not yet implemented) |

**Usage**

```
python3 fudgeify.py  start_006_fullgeo.txt  [maxadd=1]  [outfile=phenix_fudge.eff]
                     [table=fudge_table.txt]  [overridefile=overrides.txt]
```

Key-value arguments (`key=val`) are accepted alongside normal flags.

---

### `make_ss_link.com`

Detects CYS–CYS pairs with SG–SG < 2.6 Å and writes `ss_link.cif` /
`ss_link.eff` using the NEW_SS geometry (CB–SG–SG = 100.2°, not Phenix's
default 104.2°).  Uses `apply_cif_link` because
`geometry_restraints.edits` cannot override SS angles in Phenix 2.1.

---

## Key files created at runtime

| File | Created by | Purpose |
|------|-----------|---------|
| `noocc.eff` | ground_factory (first run) | Disable occupancy refinement; generate R-free flags |
| `ss_link.cif/eff` | make_ss_link | SS-bond geometry (if SS bonds present) |
| `default.geo` | phenix iteration 001 | Reference geometry for redo_harmonics |
| `fudge_table.txt` | fudgeify | Sigma state across cycles |
| `overrides.txt` | user-supplied or empty | TORSION OVERRIDE: genuine outliers to skip |

---

## Typical run

```tcsh
cd /my/target/dir
/path/to/counterfeit/ground_factory_runme.com \
    model.pdb  data.mtz \
    fudgeify=1  maxcycles=300  srun=1
```

Resume after interruption: re-run the same command; the pipeline skips
iterations whose `start_NNN.pdb` already exists and picks up from the last
completed one.

---

## Notes

- All scripts must be invoked directly (`./script.com`), **not** via
  `tcsh script.com`, so the shebang is honoured.
- MTZ files without R-free flags are handled automatically
  (`data_manager.fmodel.xray_data.r_free_flags.generate = True` in `noocc.eff`).
- Phenix planarity restraints via `geometry_restraints.edits` are unreliable
  for peptide planes; plane corrections go through REFMAC
  `extern plane` entries instead.
- OVERRIDE policy: add a residue to `overrides.txt` only when Fo−Fc difference
  density appears near a tightened restraint — not based on wE score alone.
