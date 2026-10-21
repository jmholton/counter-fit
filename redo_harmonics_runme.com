#! /bin/tcsh -f
#
#   redo_harmonics_runme.com  - James Holton
#
#   Update harmonic holdback restraints for the geometry-correction loop.
#   Compares bad atoms in the START model vs THIS iteration; atoms that became
#   bad only in this iteration get harmonic positional restraints (sigma=0.1 A)
#   to prevent X-ray from pulling them further.  Waters that moved >0.1 A also
#   get harmonic restraints.
#
#   Writes:
#     phenix_holdme.eff       - reference_coordinate_restraints for phenix.refine
#     refmac_opts_holdback.txt / refmac_opts_new.txt  - REFMAC harmonic restraints
#
#   Usage: redo_harmonics_runme.com [startpdb=model.pdb] [prefix=start]

set startpdb = start_heavyweight.pdb
set prefix   = start

foreach Arg ( $* )
  set key = `echo $Arg | awk -F= '{print $1}'`
  set val = `echo $Arg | awk -F= '{print $2}'`
  if("$key" == "startpdb") set startpdb = $val
  if("$key" == "prefix")   set prefix   = $val
end

set thispdb = `ls -1rt ${prefix}_???.pdb | tail -n 1`
if("$thispdb" == "") then
  echo "ERROR: no ${prefix}_???.pdb found"
  exit 1
endif
set num = `echo $thispdb | awk -F "_" '{print substr($NF,1,3)}'`

if(! -e refmac_opts.txt) touch refmac_opts.txt

if(! -e this/refmac_opts_debump.txt) then
  `dirname $0`/molprobify_runme.com $thispdb keepgeo=1 writefudge=1
  mkdir -p this
  cp refmac_opts_de*.txt refmac_opts_fix*.txt this/
endif

foreach reason ( `ls this | awk -F "[_.]" '{print $(NF-1)}'` )
foreach stage ( start this )
cat ${stage}/refmac_opts_${reason}.txt |\
awk '/extern dist/{\
  print $5,$7,$9,$11;\
  print $14,$16,$18,$20}\
  /extern angle/{\
  print $5,$7,$9,$11;\
  print $14,$16,$18,$20;\
  print $23,$25,$27,$29}\
  /extern tor/{\
  print $5,$7,$9,$11;\
  print $14,$16,$18,$20;\
  print $23,$25,$27,$29;\
  print $32,$34,$36,$38;}' |\
sort -u |\
sort >! badatoms_${stage}.txt
end

awk '{print "THEN",$0}' badatoms_start.txt |\
cat - badatoms_this.txt |\
awk '/^THEN/{++seen[substr($0,6)];next}\
  ! seen[$0] {print}' |\
sort -k2g |\
tee hold_these_${reason}.txt

end

rmsd -v debug=1 $startpdb $thispdb |\
 awk '$6>0.1 && /HOH/ && /moved/{print $3,$4,$1,substr($2,1,1)}' |\
tee -a hold_these_water.txt

sort -u hold_these_*.txt |\
awk '!($1=="B" && $2=="30" && ($3~/^C[BGDE]/ || $3=="NZ" || $3~/^H[GDEZ]/)) && !($1=="B" && $2=="12" && ($3~/^[NCO]$/ || $3=="CA" || $3=="CB" || $3=="HA")) && !($1~/^[AB]$/ && ($3=="N" || $3=="CA" || $3=="C"))' |\
sort -k1,1 -k2,2g |\
awk '{\
  print "extern harmonic atinfo chain",$1,"resi",$2,"atom",$3,"alt",$4,"sigma",0.1;\
  }' |\
tee refmac_opts_holdback.txt

egrep "^extern harmonic" refmac_opts.txt |\
awk '!($5=="B" && $7=="30" && ($9~/^C[BGDE]/ || $9=="NZ" || $9~/^H[EGDZ]/)) && !($5=="B" && $7=="12" && ($9~/^[NCO]$/ || $9=="CA" || $9=="CB" || $9=="HA")) && !($5~/^[AB]$/ && ($9=="N" || $9=="CA" || $9=="C"))' >! previous.txt
sort -u previous.txt refmac_opts_holdback.txt |\
sort -k5,5 -k7,7g |\
cat >! new_harmonic.txt

egrep -v "extern harm" refmac_opts.txt >! refmac_opts_new.txt
cat new_harmonic.txt >> refmac_opts_new.txt

cat << EOF >! phenix_holdme.eff
refinement {
  pdb_interpretation {
    c_beta_restraints = True
    chir_volume_esd = 0.1
    reference_coordinate_restraints {
      enabled = True
      exclude_outliers = False
EOF
echo last |\
cat new_harmonic.txt - |\
awk 'BEGIN{\
   printf("      selection = \"");}\
  /extern harmonic/{sigma=$NF;\
     c=$5;r=$7;a=$9;f=$11;\
     ++n;if(n>1)printf("\n        or ");\
     printf("( name %s and resseq %s and chain %s and altid %s) ",a,r,c,f);\
  }\
  $0=="last"{print "\""}\
  END{\
     print "      sigma =",sigma;\
     print "    }";\
     print "  }";\
     print "}";\
   }' |\
cat >> phenix_holdme.eff
