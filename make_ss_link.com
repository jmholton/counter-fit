#! /bin/tcsh -f
#
#   make_ss_link.com  - James Holton  6-27-26
#
#   Auto-generate ss_link.cif and ss_link.eff from a PDB file.
#   Detects CYS-CYS pairs with SG-SG distance < 2.6 A and writes
#   apply_cif_link blocks for each conformer.
#
#   Uses NEW_SS geometry: CB-SG-SG = 100.2 deg (vs Phenix default 104.2).
#   geometry_restraints.edits cannot override SS angles in Phenix 2.1;
#   apply_cif_link is the only path.
#
#   Usage:  make_ss_link.com model.pdb
#   Output: ss_link.cif  ss_link.eff

set pdbfile = ""
foreach arg ( $* )
  if("$arg" =~ *.pdb) set pdbfile = $arg
end
if("$pdbfile" == "") then
  echo "Usage: make_ss_link.com model.pdb"
  exit 1
endif

set tmpdir = /tmp/make_ss_link_$$
mkdir -p $tmpdir

# ---- awk script to find SS pairs ----
cat >! ${tmpdir}/find_ss.awk << EOF
/^ATOM|^HETAT/ {
  if (\$3 != "SG") next
  alt = substr(\$0,17,1)
  res = substr(\$0,18,3)
  ch  = substr(\$0,22,1)
  rn  = substr(\$0,23,4) + 0
  if (res !~ /CYS/) next
  key = ch "_" rn "_" alt
  x[key] = \$7+0; y[key] = \$8+0; z[key] = \$9+0
  chain[key] = ch; resnum[key] = rn; altid[key] = alt
}
END {
  for (i in x) { for (j in x) {
    if (i < j) {
      d = sqrt((x[i]-x[j])^2 + (y[i]-y[j])^2 + (z[i]-z[j])^2)
      if (d >= 1.8 && d <= 2.6) {
        ai = altid[i]; aj = altid[j]
        if (ai == "" || aj == "" || ai == aj) {
          printf "%s %s %s %s %s %s\n", chain[i], resnum[i], (ai==""?" ":ai), \
                                         chain[j], resnum[j], (aj==""?" ":aj)
        }
      }
    }
  }}
}
EOF

awk -f ${tmpdir}/find_ss.awk $pdbfile >! ${tmpdir}/ss_pairs.txt

if (! -s ${tmpdir}/ss_pairs.txt) then
  echo "No SS bonds found in $pdbfile (SG-SG 1.8-2.6 A) — no files written"
  rm -rf $tmpdir
  exit 0
endif

set npairs = `wc -l < ${tmpdir}/ss_pairs.txt`
echo "Found $npairs SS bond(s):"
cat ${tmpdir}/ss_pairs.txt

# ---- Write ss_link.cif ----
cat >! ss_link.cif << EOF
data_link_NEW_SS

loop_
  _chem_link_bond.link_id
  _chem_link_bond.atom_1_comp_id
  _chem_link_bond.atom_id_1
  _chem_link_bond.atom_2_comp_id
  _chem_link_bond.atom_id_2
  _chem_link_bond.type
  _chem_link_bond.value_dist
  _chem_link_bond.value_dist_esd
  NEW_SS  1  SG  2  SG  single  2.031  0.020

loop_
  _chem_link_angle.link_id
  _chem_link_angle.atom_1_comp_id
  _chem_link_angle.atom_id_1
  _chem_link_angle.atom_2_comp_id
  _chem_link_angle.atom_id_2
  _chem_link_angle.atom_3_comp_id
  _chem_link_angle.atom_id_3
  _chem_link_angle.value_angle
  _chem_link_angle.value_angle_esd
  NEW_SS  1  CB  1  SG  2  SG  100.200  2.100
  NEW_SS  1  SG  2  SG  2  CB  100.200  2.100

loop_
  _chem_link_tor.link_id
  _chem_link_tor.id
  _chem_link_tor.atom_1_comp_id
  _chem_link_tor.atom_id_1
  _chem_link_tor.atom_2_comp_id
  _chem_link_tor.atom_id_2
  _chem_link_tor.atom_3_comp_id
  _chem_link_tor.atom_id_3
  _chem_link_tor.atom_4_comp_id
  _chem_link_tor.atom_id_4
  _chem_link_tor.value_angle
  _chem_link_tor.alt_value_angle
  _chem_link_tor.value_angle_esd
  _chem_link_tor.period
  NEW_SS  ss      1  CB  1  SG  2  SG  2  CB   93.0  -86.0        10.0  1
  NEW_SS  chi2_1  1  CA  1  CB  1  SG  2  SG   79.0  183.0,-73.0  20.0  1
  NEW_SS  chi2_2  1  SG  2  SG  2  CB  2  CA  -79.0 -183.0,73.0   20.0  1
EOF

# ---- Write ss_link.eff ----
echo "refinement.pdb_interpretation {" >! ss_link.eff

# awk script to generate apply_cif_link blocks
cat >! ${tmpdir}/gen_eff.awk << EOF
{
  ch1=\$1; rn1=\$2; al1=\$3
  ch2=\$4; rn2=\$5; al2=\$6
  sel1 = "chain " ch1 " and resname CYS and resseq " rn1
  sel2 = "chain " ch2 " and resname CYS and resseq " rn2
  if (al1 != " " && al1 != ch1) sel1 = sel1 " and altid " al1
  if (al2 != " " && al2 != ch2) sel2 = sel2 " and altid " al2
  print "  apply_cif_link {"
  print "    data_link = NEW_SS"
  print "    residue_selection_1 = " sel1
  print "    residue_selection_2 = " sel2
  print "  }"
}
EOF

awk -f ${tmpdir}/gen_eff.awk ${tmpdir}/ss_pairs.txt >> ss_link.eff
echo "}" >> ss_link.eff

rm -rf $tmpdir
echo ""
echo "Written: ss_link.cif  ss_link.eff"
