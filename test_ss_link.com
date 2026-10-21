#! /bin/tcsh -f
#
#   test_ss_link.com   - James Holton  6-20-26
#
#   Test whether the current Phenix version correctly applies the observed
#   CB-SG-SG disulfide angle (100.2 deg) automatically, or whether the
#   ss_link.cif / apply_cif_link workaround is still needed.
#
#   Background:
#     Phenix built-in data_link_SS uses CB-SG-SG = 104.2 deg (Engh & Huber).
#     The NEW_SS link in ss_link.cif uses 100.2 deg (better matches observed structures).
#     Molprobify cannot distinguish the two cases (it uses its own 104.2 deg reference);
#     the raw phenix .geo file is the authoritative source.
#
#   Tested result: Phenix 2.1rc2-6037 still needs ss_link.cif.
#
#   Usage:
#     test_ss_link.com model.pdb data.mtz [ss_link.cif] [ss_link.eff]

set pdbfile = ""
set mtzfile = ""
set ciflinkfile = ""
set effinkfile = ""

foreach arg ( $* )
  if("$arg" =~ *.pdb) set pdbfile = $arg
  if("$arg" =~ *.mtz) set mtzfile = $arg
  if("$arg" =~ *link*.cif || "$arg" =~ *ss*.cif) set ciflinkfile = $arg
  if("$arg" =~ *link*.eff || "$arg" =~ *ss*.eff) set effinkfile = $arg
end

if("$pdbfile" == "" || "$mtzfile" == "") then
  echo "Usage: test_ss_link.com model.pdb data.mtz [ss_link.cif] [ss_link.eff]"
  exit 1
endif
if("$ciflinkfile" == "" && -e ss_link.cif) set ciflinkfile = ss_link.cif
if("$effinkfile"  == "" && -e ss_link.eff) set effinkfile  = ss_link.eff

# Check the model actually has SS bonds
set nSS = `awk '/^ATOM|^HETAT/{if($3=="SG"&&$4=="CYS"){x[$5$6]=$7;y[$5$6]=$8;z[$5$6]=$9}}\
  END{n=0;for(i in x){for(j in x){if(i>=j)next;\
    d=sqrt((x[i]-x[j])^2+(y[i]-y[j])^2+(z[i]-z[j])^2);\
    if(d>1.8&&d<2.6)n++}};print n}' $pdbfile`
if($nSS == 0) then
  echo "ERROR: no SS bonds found in $pdbfile (SG-SG 1.8-2.6 A range)"
  exit 1
endif
echo "Found $nSS SS bond(s) in $pdbfile"

# Write noocc.eff if not present
if(! -e noocc.eff) then
cat >! noocc.eff << EOF
refinement {
  pdb_interpretation { c_beta_restraints = True; chir_volume_esd = 0.1 }
  refine {
    strategy = *individual_sites *individual_adp occupancies
    occupancies { remove_selection = All }
  }
  geometry_restraints.edits { excessive_bond_distance_limit = 100 }
}
EOF
endif

echo ""
echo "=== Run 1: WITHOUT ss_link.cif ==="
srun -p debug phenix.refine $pdbfile $mtzfile \
  prefix=nolink serial=001 \
  main.number_of_macro_cycles=1 \
  write_geo_file=True wxc_scale=0.1 wxu_scale=0.05 \
  noocc.eff >& nolink.log &

echo "=== Run 2: WITH ss_link.cif + ss_link.eff ==="
if("$ciflinkfile" != "" && "$effinkfile" != "") then
  srun -p debug phenix.refine $pdbfile $mtzfile \
    $ciflinkfile $effinkfile \
    prefix=withlink serial=002 \
    main.number_of_macro_cycles=1 \
    write_geo_file=True wxc_scale=0.1 wxu_scale=0.05 \
    noocc.eff >& withlink.log &
else
  echo "  (no ss_link.cif/eff — skipping run 2)"
endif

wait
echo ""
echo "=== CB-SG-SG ideal angles in raw phenix .geo files ==="
foreach geo ( nolink_001.geo withlink_002.geo )
  if(! -e $geo) continue
  echo "--- $geo ---"
  awk '/^angle pdb=/{a1=$0;getline;a2=$0;getline;a3=$0;getline;getline;\
    if((a1~/CB.*CYS/&&a2~/SG.*CYS/&&a3~/SG.*CYS/)||\
       (a1~/SG.*CYS/&&a2~/SG.*CYS/&&a3~/CB.*CYS/))\
      printf "  ideal=%s sigma=%s | %s | %s | %s\n",$1,$4,a1,a2,a3}' $geo | head -3
end

echo ""
echo "=== Summary ==="
set nolink_angle = `awk '/^angle pdb=/{a1=$0;getline;a2=$0;getline;a3=$0;getline;getline;\
  if((a1~/CB.*CYS/&&a2~/SG.*CYS/&&a3~/SG.*CYS/)||\
     (a1~/SG.*CYS/&&a2~/SG.*CYS/&&a3~/CB.*CYS/))\
    {print $1;exit}}' nolink_001.geo 2>/dev/null`

if("$nolink_angle" == "") then
  echo "No CB-SG-SG angle found without ss_link.cif — Phenix not applying SS geometry."
  echo "CONCLUSION: ss_link.cif + ss_link.eff are REQUIRED."
else if(`echo $nolink_angle | awk '{print ($1>103&&$1<105)}'`) then
  echo "Without ss_link.cif: CB-SG-SG ideal = $nolink_angle deg (Phenix default 104.2)"
  echo "CONCLUSION: ss_link.cif + ss_link.eff STILL NEEDED to get 100.2 deg."
else if(`echo $nolink_angle | awk '{print ($1>99&&$1<101)}'`) then
  echo "Without ss_link.cif: CB-SG-SG ideal = $nolink_angle deg — already 100.2."
  echo "CONCLUSION: ss_link.cif + ss_link.eff NOT NEEDED in this Phenix version."
else
  echo "Without ss_link.cif: CB-SG-SG ideal = $nolink_angle deg (unexpected — inspect manually)"
endif

echo ""
echo "Phenix version: `phenix.refine --version 2>&1 | head -1`"
