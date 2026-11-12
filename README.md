# counterfeit — Ground Factory Pipeline

Iteratively optimize a crystal structure model to near-perfect geometry while
maintaining a reasonable fit to X-ray data.

**Strategy:** Run `phenix.refine` with a heavily down-weighted X-ray term
(`wxc_scale=0.1`) so geometry restraints dominate.  Each cycle: validate with
MolProbity, update harmonic hold-back restraints on atoms that are getting
worse, then add or tighten one geometry-correction restraint per outlier type.
Repeat until the weighted energy score converges (target wE < 9).

Every cycle refines the same starting model, and sigmas change gradually
(÷1.5 for the worst outliers, ×1.01 for everything else), because the system
is prone to instability.  Changes are never undone: if wE goes up, something
else needs more restraint, and the next cycle gives it that.

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
  crowbars=N         with fudgeify, clashes and contacts crowbarred per cycle,
                     each (default maxadd; 0 = no crowbars)
  bonds=N  angles=N  torsions=N  planes=N
                     with fudgeify, outliers of that kind restrained per cycle
                     (default maxadd; torsions include CB deviations, planes
                     include the impropers that stand in for peptide planes)
  cdl=True           phenix conformation-dependent library in molprobify scoring
  cdl_interp=1       CDL ideals interpolated between phi/psi bins
                     (cdl_interp.py); 0 = binned, as the Challenge scorer
  refine_cdl=False   CDL in phenix.refine (CDL undoes backbone bond edits, so
                     fudgeify passes the scorer's targets on instead)
  fixweights=1       fix phenix's X-ray weights at iteration 001's values
                     (see below); 0 = phenix recomputes them every macro-cycle
  fix_wxc=X          with fixweights=1, pass X as fix_wxc instead of iteration
                     001's wxc × wxc_scale (snap iterations: X × snap_wxc / wxc)
  omegas=N           with fudgeify, omegas restrained per cycle (default maxadd)
  omega_esd=X        phenix's sigma for every omega, degrees (phenix 5;
                     molprobify scores at 4.01)
  overrides=file     TORSION OVERRIDE file — genuine outliers that must not
                     be tightened (see OVERRIDE policy below)
```

**Outputs**

| File | Description |
|------|-------------|
| `start_NNN.pdb` | phenix model each iteration |
| `itr_vs_energy.txt` | weighted energy vs iteration number (`snap`, `snap_reverted` tags) |
| `phenix_holdme.eff` | harmonic positional restraints (good atoms held back) |
| `phenix_fudge.eff` | geometry-correction restraints (problem atoms) |
| `refmac_opts.txt` | same restraints in REFMAC external format |
| `fudge_table.txt` | fudgeify state: sigma history across cycles |
| `start_NNN_crowbars.txt` | contacts each crowbar switches off or weakens in phenix |
| `ground_factory.lock` | host and PID of the run using this directory |
| `fixed_weights.txt` | wxc and wxu from iteration 001, unscaled |
| `fmodelme.pdb` | final model after REFMAC B-factor convergence |
| `nonoise.mtz` / `simplenoise.mtz` | ideal / ~3%-noisy structure factors, at the resolution of the input mtz |

**Loop per iteration**

1. `phenix.refine` with `wxc_scale=0.1`
2. `molprobify_runme.com` → weighted energy, fullgeo file, restraint
   proposals (`refmac_opts_*.txt`, read by steps 3 and 4)
3. `redo_harmonics_runme.com` → update `phenix_holdme.eff`
4. `reweight_opts_runme.com` → update `refmac_opts.txt`; writes `phenix_fudge.eff`
   (overwritten by fudgeify when `fudgeify=1`)
5. `fudgeify.py` (if `fudgeify=1`) → update `fudge_table.txt`, overwrite `phenix_fudge.eff`
6. 1 % sigma escalation on all entries in `refmac_opts.txt`

The loop stops at `maxcycles`, or when phenix, molprobify or fudgeify fails.
The final REFMAC step starts from the last iteration.

**X-ray weights.**  phenix normally recomputes its X-ray weights (wxc, wxu)
every macro-cycle from gradient norms, and a new or tightened restraint can
change them more than tenfold from one iteration to the next (7KR0: wxc 23 →
334), which undoes the down-weighting.  With `fixweights=1` they are fixed
at their averages over iteration 001's macro-cycles (in `fixed_weights.txt`),
and every later iteration passes `fix_wxc = wxc × wxc_scale` (phenix ignores
`wxc_scale` once wxc is fixed) and `fix_wxu = wxu` (phenix still applies
`wxu_scale`).  phenix's default target is least-squares for the first two
macro-cycles and ml after, with weights ~10× apart, so every iteration uses
`main.target=ml`.

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

Ideal values come from `phenix.geometry_minimization`, with CDL on unless
`cdl=False` is given (ground_factory passes its own setting).

Types: `NONBOND`, `BOND`, `ANGLE`, `TORSION`, `PLANE`.  A contact with a
symmetry mate ends in `symop=<op>`: the `.geo` `sym.op.`, applied to the
second atom.  Every restraint made from such a contact carries the operator
(`symmetry_operation` for phenix, `symm y` for REFMAC, with the explicit
operator in a trailing `! symop` comment).

---

### `redo_harmonics_runme.com`

Compares the starting model to the current iteration output.  Atoms that became
worse (appeared in bad-atom lists for the first time) get harmonic positional
restraints (`sigma = 0.1 Å`) in `phenix_holdme.eff` so the X-ray term cannot
pull them further from ideal geometry.  Waters that drifted > 0.1 Å are also
restrained.  When both atoms of a crowbar (a `NONBOND` entry in
`fudge_table.txt`) are held, both are released, since the crowbar could not
move either.

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

Reads a `*_fullgeo.txt` and the matching `*_clashes.txt`; accumulates sigma
state in `fudge_table.txt`; writes `phenix_fudge.eff`.  Each call is one
cycle:

| Event | Action |
|-------|--------|
| One of the worst `maxadd` outliers of its type, not yet in the table | Add: a crowbar at sigma 0.4; a plane at phenix's sigma (it is a second restraint); BOND/ANGLE/TORSION at natural sigma/1.5 (they change phenix's own restraint, so natural sigma would change nothing) |
| One of the worst `maxadd` outliers of its type, already in the table | Tighten: `sigma /= 1.5`, down to the floors (crowbar 0.01, BOND 0.001, ANGLE and TORSION 0.01) |
| A crowbar chosen again at its floor | Lengthen: `dmin += 0.1` |
| Anything else in the table | Relax: `sigma *= 1.01` |
| NONBOND, BOND or PLANE sigma ≥ 5 | Drop (restraint resolved) |

**Crowbars** (NONBOND entries) push contacts apart: a phenix bond that only
pushes, flat from `dmin` out to `dmin+2` (`distance_ideal = dmin+1`,
`slack = 1`).  Probe clashes and scored contacts (fullgeo `NONBOND` lines)
are two separate lists, each giving its worst `crowbars` per cycle (default
`maxadd`).  A
hydrogen is replaced by its parent heavy atom, since crowbarring a hydrogen
only bends it, and `dmin` is the parents' contact distance plus the overlap
(plus 0.1 Å per hydrogen replaced, for scored contacts).  Parents that
phenix does not treat as a contact (bonded or angle neighbours) get no
crowbar.

A crowbar goes into phenix's bond graph, so phenix switches off the contacts
between each crowbar atom and the other atom's bonded neighbours (1-3) and
weakens the next shell (1-4, vdW × 0.8).  They are still scored.  Each cycle
`start_NNN_crowbars.txt` lists them under their crowbar, with the crowbar's
current distance and the contact distance, and the log counts how many are
now too close.  Lengthening has no limit, so a contact that never resolves
keeps getting longer; that file shows it.

**Planes.**  fullgeo lists plane deviations per atom without saying which
plane, so fudgeify measures every phenix plane (from `default.geo`) in the
model (`<prefix>.pdb`); a plane is out when an atom is out of it both there
and in the fullgeo score (so overrides apply).  The worst `maxadd` get a
second planarity restraint on the plane's heavy atoms.  phenix refuses a
plane with the same atoms as one of its own, so such a copy drops its
best-fitting atom, and a 4-atom one (a peptide's CA, C, O, N) is restrained
instead by an improper dihedral on the same 4 atoms (0 or 180°, sigma =
the plane's sigma over the 4th atom's distance from the dihedral axis, in
degrees), which is the only way 4 atoms leave a plane.

molprobify scores with phenix's conformation-dependent library on
(`cdl=True`), as the Challenge scorer does, and by default interpolated
(`cdl_interp=1`, below).  Refinement runs with CDL off (`refine_cdl=False`):
CDL sets the C–N, N–CA, CA–CB, CA–C and C–O bonds from phi/psi after the
edits are applied, and again every macro-cycle, so with CDL on in refinement
an edit to one of those bonds does nothing.  Instead fudgeify passes the
scorer's targets on: every cycle each BOND and ANGLE entry takes the scorer's
current ideal from fullgeo, and every other bond and angle whose scorer ideal
or sigma differs from phenix's own (`default.geo`) gets an edit with the
scorer's values, rebuilt each cycle.  So refinement aims at the CDL targets
of the model's current phi/psi.  Runs before 12 Nov 2026 scored with
`cdl=False`, so their wE is not the Challenge scorer's.  Bonds between
symmetry mates are not changed.

phenix takes CDL ideals from the phi/psi bin a residue is in, so they jump at
bin edges.  `mmtbx.conformation_dependent_library.update_restraints` can
interpolate between bins instead (`interpolate=True`), but phenix exposes no
parameter for it, so `cdl_interp.py` turns it on and then runs the program:
`phenix.python cdl_interp.py phenix.command_line.phenix_refine …` (or
`mmtbx.command_line.geometry_minimization`).  ground_factory passes
`cdl_interp=1` to molprobify, which uses it for geometry_minimization, and
uses it for phenix.refine only with `refine_cdl=True`.  molprobify on its
own defaults to `cdl_interp=0`, as the Challenge scorer
(`~/Develop/molprobify_runme.com`) does.

Crowbars between symmetry mates keep the operator (`sym=` in the table,
`symmetry_operation` in the eff).  A torsion phenix already restrains
(found in `default.geo`) is changed in place rather than added, because phenix
refuses duplicate dihedrals; fudge torsions are harmonic (`periodicity = 0`;
phenix's default of 1 means sinusoidal).  Torsions named in `overrides.txt`
are skipped (exact atom match, either direction).

**CB deviations.**  ground_factory turns phenix's `c_beta_restraints` off,
so fudgeify restrains them itself, from `<prefix>_cbetadev.log` (passed as
`cbetadev=`).  Each deviation is scored as (dev/0.05)², as molprobify scores
it, and ranked with the fullgeo outliers; the worst `torsions` per cycle get
the pair of pseudo-dihedrals phenix would use (N–C–CA–CB and C–N–CA–CB, ideals
from phenix `c_beta.py`, sigma 2.5), added as new restraints and tightened
like the others.  The N/C–CA–CB torsions in fullgeo are skipped, since these
pairs replace them.

**Omegas.**  The CA–C–N–CA torsions are their own category (`omegas` per
cycle, default maxadd), ranked by molprobify's omega energy: its sigma is 4.01°
where phenix restrains omega at 5°, so an omega phenix calls fine still scores,
and the scored average over all omegas sits at about (5/4.01)² ≈ 1.5.  One per
cycle cannot move that average; ground_factory `omega_esd=X` sets phenix's
sigma for every omega (`pdb_interpretation.peptide_link.omega_esd_override_value`).
Cis peptides are left out.

**Usage**

```
python3 fudgeify.py  start_006_fullgeo.txt  [maxadd=1]  [crowbars=1]  [outfile=phenix_fudge.eff]
                     [bonds=1]  [angles=1]  [torsions=1]  [planes=1]  [omegas=1]
                     [table=fudge_table.txt]  [overridefile=overrides.txt]
                     [clashes=start_006_clashes.txt]  [watch=start_006_crowbars.txt]
                     [pdb=start_006.pdb]  [geo=default.geo]
                     [cbetadev=start_006_cbetadev.log]
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
| `noocc.eff` | ground_factory (first run) | Disable occupancy refinement; isotropic B for every non-hydrogen atom (phenix refuses riding H in the selection); generate R-free flags |
| `ss_link.cif/eff` | make_ss_link | SS-bond geometry (if SS bonds present) |
| `default.geo` | phenix iteration 001 | Reference geometry for redo_harmonics |
| `fudge_table.txt` | fudgeify | Sigma state across cycles |
| `overrides.txt` | user-supplied before the run, or empty | TORSION OVERRIDE: genuine outliers to skip |

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
completed one.  A second run in the same directory refuses to start while
`ground_factory.lock` names a live process; delete the lock if that run is
dead on another host.

---

## Notes

- All scripts must be invoked directly (`./script.com`), **not** via
  `tcsh script.com`, so the shebang is honoured.
- MTZ files without R-free flags are handled automatically
  (`data_manager.fmodel.xray_data.r_free_flags.generate = True` in `noocc.eff`).
- phenix refuses a planarity edit with the same atoms as one of its own
  planes (such as the peptide's CA, C, O, N).  reweight_opts skips those;
  fudgeify restrains them with an improper dihedral instead.
- OVERRIDE policy: `overrides.txt` is populated only at the beginning, from
  obvious outliers in obviously good density.  Nothing adds to it during a
  run.
- tcsh has no `2>/dev/null`: inside backticks it either swallows the output
  or fails with "Ambiguous output redirect".  Test for the file first.
