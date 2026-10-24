#! /bin/tcsh -f
#
#   ground_factory_runme.com  - James Holton  6-27-26
#
#   Iteratively optimize a crystal structure model to have near-perfect
#   geometry while maintaining reasonable fit to X-ray data.
#
#   Strategy: phenix.refine with heavily downweighted X-ray term (wxc=0.1)
#   so geometry restraints dominate.  Each cycle: validate with molprobify,
#   update harmonic holdback restraints on good atoms, add/tighten one geometry
#   correction restraint per outlier type.  Run until energy converges.
#
#   Usage:
#     ground_factory_runme.com model.pdb data.mtz [options]
#
#   Options:
#     maxcycles=N        max iterations (default 300)
#     maxadd=N           new restraints per type per cycle (default 1)
#     wxc=0.1            X-ray coordinate weight scale
#     wxu=0.05           X-ray ADP weight scale
#     overrides=file     TORSION OVERRIDE file for genuine outliers
#     prefix=start       phenix output prefix
#     srun=1             submit phenix jobs via "srun -p debug" (default 1)
#     propagate=0        if 1, each iteration starts from previous output (default 0)
#     snap_interval=0    if >0, run a snap iteration every N cycles (default 0=off)
#     snap_wxc=1.0       wxc_scale used during snap iterations (default 1.0)
#
#   Outputs:
#     start_NNN.pdb    phenix-refined model each iteration
#     fmodelme.pdb     final model after REFMAC B-factor convergence
#     fmodel.mtz       ideal structure factors (phenix.fmodel)
#     nonoise.mtz      clean Fcalc + free-R flags
#     simplenoise.mtz  ~3% noise added
#     itr_vs_energy.txt  energy per iteration

set pdbfile       = ""
set mtzfile       = ""
set maxcycles     = 300
set maxadd        = 1
set wxc           = 0.1
set wxu           = 0.05
set overrides     = ""
set prefix        = start
set srun          = 1
set propagate     = 0
set snap_interval = 0
set snap_wxc      = 1.0

foreach Arg ( $* )
  set key = `echo $Arg | awk -F= '{print $1}'`
  set val = `echo $Arg | awk -F= '{print $2}'`
  if("$Arg" =~ *.pdb) set pdbfile = $Arg
  if("$Arg" =~ *.mtz) set mtzfile = $Arg
  if("$key" == "maxcycles")     set maxcycles     = $val
  if("$key" == "maxadd")        set maxadd        = $val
  if("$key" == "wxc")           set wxc           = $val
  if("$key" == "wxu")           set wxu           = $val
  if("$key" == "overrides")     set overrides     = $val
  if("$key" == "prefix")        set prefix        = $val
  if("$key" == "srun")          set srun          = $val
  if("$key" == "propagate")     set propagate     = $val
  if("$key" == "snap_interval") set snap_interval = $val
  if("$key" == "snap_wxc")      set snap_wxc      = $val
end

if("$pdbfile" == "" || "$mtzfile" == "") then
  echo "Usage: ground_factory_runme.com model.pdb data.mtz [options]"
  exit 1
endif

set scriptdir = `dirname $0`

set PHENIXRUN = "phenix.refine"
if($srun) set PHENIXRUN = "srun -p debug phenix.refine"

# -----------------------------------------------------------------------
# Step 1: SS bonds
# -----------------------------------------------------------------------
if(! -e ss_link.cif || ! -e ss_link.eff) then
  echo "=== Detecting SS bonds ==="
  ${scriptdir}/make_ss_link.com $pdbfile
endif

set ssargs = ""
if(-s ss_link.cif && -s ss_link.eff) set ssargs = "ss_link.cif ss_link.eff"

# -----------------------------------------------------------------------
# Step 2: Static control files
# -----------------------------------------------------------------------
if(! -e noocc.eff) then
cat >! noocc.eff << EOF
refinement {
  pdb_interpretation {
    c_beta_restraints = True
    chir_volume_esd = 0.1
  }
  refine {
    strategy = *individual_sites individual_sites_real_space rigid_body \
               *individual_adp group_adp tls occupancies group_anomalous den
    occupancies { remove_selection = All }
  }
  geometry_restraints.edits {
    excessive_bond_distance_limit = 100
  }
}
EOF
endif

if(! -e phenix_holdme.eff) then
  echo 'refine {}' >! phenix_holdme.eff
endif
if(! -e phenix_fudge.eff) then
  echo 'refine {}' >! phenix_fudge.eff
endif
if(! -e refmac_opts.txt)   touch refmac_opts.txt
if(! -e overrides.txt && "$overrides" != "") cp $overrides overrides.txt
if(! -e overrides.txt) touch overrides.txt

set overridearg = ""
if(-s overrides.txt) set overridearg = "overridefile=overrides.txt"

# -----------------------------------------------------------------------
# Step 3: Validate starting model, save reference geometry
# -----------------------------------------------------------------------
echo "=== Validating starting model ==="
mkdir -p this start
if(! -e start/refmac_opts_debump.txt) then
  ${scriptdir}/molprobify_runme.com $pdbfile keepgeo=1 writefudge=1 $overridearg
  cp refmac_opts_de*.txt refmac_opts_fix*.txt start/
endif

if(! -e default.geo) then
  echo "=== Running first phenix cycle to get reference .geo ==="
  $PHENIXRUN $pdbfile $mtzfile \
    prefix=${prefix} serial=001 \
    $ssargs noocc.eff \
    wxc_scale=${wxc} wxu_scale=${wxu} \
    main.number_of_macro_cycles=5 \
    write_geo_file=True \
    phenix_holdme.eff phenix_fudge.eff >! phenix_001.log
  if($status) then
    echo "ERROR: initial phenix.refine failed — see phenix_001.log"
    exit 1
  endif
  cp ${prefix}_001.geo default.geo
endif

# -----------------------------------------------------------------------
# Step 4: Iterative geometry-correction loop
# -----------------------------------------------------------------------
set num = `ls -1rt ${prefix}_???.pdb | tail -n 1 | awk -F_ '{print substr($NF,1,3)+0}'`
if("$num" == "") set num = 0

# For forward propagation, resume from the last completed iteration's output
set startpdb = $pdbfile
if($propagate) then
  set _last = `ls -1rt ${prefix}_???.pdb 2>/dev/null | tail -1`
  if("$_last" != "") set startpdb = $_last
endif

echo "=== Starting iterative loop from iteration $num ==="
if($num == 0) echo "" >! itr_vs_energy.txt

foreach x ( `seq 1 $maxcycles` )

  set num = `echo $num | awk '{printf("%03d",$1+1)}'`

  # Determine wxc for this iteration (snap if interval is set and num is a multiple)
  set wxc_now = $wxc
  set is_snap = `echo $num $snap_interval | awk '{print ($2>0 && $1%$2==0)}'`
  if($is_snap) then
    set wxc_now = $snap_wxc
    echo "=== SNAP iteration $num : wxc_scale=$snap_wxc ==="
  endif

  # Skip if this iteration's output already exists
  if(-e ${prefix}_${num}.pdb) goto already_done

  $PHENIXRUN $startpdb $mtzfile \
    prefix=${prefix} serial=${num} \
    $ssargs noocc.eff \
    wxc_scale=${wxc_now} wxu_scale=${wxu} \
    main.number_of_macro_cycles=5 \
    write_geo_file=True \
    phenix_holdme.eff phenix_fudge.eff >! phenix_${num}.log
  if($status) then
    echo "phenix.refine failed at iteration $num — stopping"
    break
  endif
  # NFS cache may lag after srun returns; wait up to 60s for the file to appear
  set _wait = 0
  while(! -e ${prefix}_${num}.pdb && $_wait < 60)
    sleep 1
    @ _wait++
  end
  if(! -e ${prefix}_${num}.pdb) break

  already_done:

  # Propagate: next iteration starts from this iteration's output
  if($propagate) set startpdb = ${prefix}_${num}.pdb

  ${scriptdir}/molprobify_runme.com ${prefix}_${num}.pdb keepgeo=1 writefudge=1 debug=1 \
    tempfile=/dev/shm/${USER}/gf_${num}_ \
    $overridearg | tee molprobify_${num}.log

  # Extract energy and log it
  set energy = `grep "weighted energy" molprobify_${num}.log | awk '{print $NF}'`
  if($is_snap) then
    echo "$num $energy snap" >> itr_vs_energy.txt
  else
    echo "$num $energy" >> itr_vs_energy.txt
  endif
  echo "Iteration $num : energy = $energy"

  mkdir -p this
  cp refmac_opts_de*.txt refmac_opts_fix*.txt this/

  ${scriptdir}/redo_harmonics_runme.com startpdb=$startpdb prefix=$prefix
  cp refmac_opts_new.txt refmac_opts.txt

  ${scriptdir}/reweight_opts_runme.com maxadd=$maxadd
  cp refmac_opts_new.txt refmac_opts.txt

  # Slow sigma escalation: tighten all external restraints by 1% per cycle,
  # drop distance restraints that went slack (sigma >= 5)
  cat refmac_opts.txt |\
  awk '{sc=0}\
       /^extern inter/{sc=26}\
       /^extern dist/{sc=24}\
       /^extern plane/{sc=NF}\
       sc{$sc*=1.01}\
       /^extern dist/ && $sc>=5{next}\
       /^extern inter/ && $sc>=5{next}\
       {print}' >! refmac_opts_new.txt
  mv refmac_opts_new.txt refmac_opts.txt

end

# -----------------------------------------------------------------------
# Step 5: Final B-factor convergence with REFMAC
# -----------------------------------------------------------------------
set lastpdb = `ls -1rt ${prefix}_???.pdb | tail -n 1`
echo "=== Final REFMAC B-factor convergence from $lastpdb ==="

cat << EOF >! refmac_opts_final.txt
damp 0 0.5 0
bfactor 5
solvent vdwprobe 1.1 ionprobe 1.1 rashrink 0.9
make hout Y
make hydr Y
EOF

# Fix H1 naming on N-terminal residue (VAL/ALA etc.) if present
cat $lastpdb |\
  awk '/H1 .VAL .   1 /{gsub("  H1 ","  H  ")} {print}' |\
  egrep -v LINK >! refmacme.pdb

cp refmac_opts.txt refmac_opts_beforefinal.txt
cp refmac_opts_final.txt refmac_opts.txt
converge_refmac.com $mtzfile refmacme.pdb >&! converge_B.log
cp refmac_opts_beforefinal.txt refmac_opts.txt
cp refmacout.pdb fmodelme.pdb

echo "=== Final validation ==="
${scriptdir}/molprobify_runme.com fmodelme.pdb keepgeo=1 $overridearg | tee molprobify_fmodelme.log

# -----------------------------------------------------------------------
# Step 6: Ideal structure factors
# -----------------------------------------------------------------------
set hires = `grep "_cell.length_a\|CRYST1\|Resolution" $mtzfile \
  phenix_001.log 2>/dev/null | \
  awk '/High resolution limit/{print $NF;exit} \
       /high_resolution/{print $NF;exit}' 2>/dev/null`
if("$hires" == "") set hires = `grep "high_res\|resolution" phenix_001.log | \
  awk '/high_resolution/{print $NF;exit}'`
if("$hires" == "") set hires = 1.5

echo "=== phenix.fmodel at ${hires} A ==="
rm -f fmodel.mtz
phenix.fmodel fmodelme.pdb k_sol=0.42 b_sol=45 high_res=${hires} \
  solvent_radius=1.1 shrink_truncation_radius=0.9 \
  grid_resolution_factor=1/5 file_name=fmodel.mtz | tee fmodel.log

chltofom -mtzin fmodel.mtz -mtzout hlified.mtz \
  -colin-phifom '/*/*/[PHIFMODEL,]'

rm -f freeme.mtz nonoise.mtz simplenoise.mtz
sftools << EOF
read hlified.mtz
set labels
F
PHI
HLA
HLB
HLC
HLD
calc Q col SIGF = 0.1
calc W col FOM = 0.99
write freeme.mtz col F SIGF PHI FOM HLA HLB HLC HLD
quit
EOF
FreeRer.com freeme.mtz > /dev/null
mv FreeRed.mtz nonoise.mtz

# Add ~3% noise: sigma_noise = F/30, added in quadrature with SIGF
sftools << EOF
read nonoise.mtz
calc col frac  = col F 30 /
calc col frac  = col frac 2 **
calc col SIGF  = col SIGF 2 **
calc col SIGF  = col SIGF col frac +
calc col SIGF  = col SIGF 0.5 **
calc col noise = col SIGF ran_g *
calc F col Fnew = col F col noise +
select col Fnew < 0.05
calc col noise = 0 col noise -
calc F col Fnew = col F col noise +
select all
calc F col F = col Fnew
delete col Fnew noise frac
write simplenoise.mtz
quit
EOF

echo ""
echo "=== Done ==="
echo "  Final model:         fmodelme.pdb"
echo "  Clean data:          nonoise.mtz"
echo "  Noisy data:          simplenoise.mtz"
echo "  Energy vs iteration: itr_vs_energy.txt"
