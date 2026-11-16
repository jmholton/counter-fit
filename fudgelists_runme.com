#! /bin/tcsh -f
#
#   fudgelists_runme.com  - the outlier lists that molprobify_runme.com writes with writefudge=1
#   (refmac_opts_declash, _debump, _fixbond, _fixang, _fixCBdev, _fixtorsion .txt), rebuilt from the files
#   wE_v2.py writes, so that redo_harmonics_runme.com can choose the atoms to hold (phenix_holdme.eff).
#   The blocks below are copied verbatim from code_v2's molprobify_runme.com (md5 457d87ea).
#   refmac_opts_fixplanes.txt is written empty: its "extern plane" lines need the .geo file, and
#   redo_harmonics takes no atoms from them.  Not written: refmac_opts_holdme.txt and phenix_opts_*.txt
#   (only reweight_opts_runme.com reads them, and ground_factory skips it with scorer=wE_v2).
#
#   Usage: fudgelists_runme.com P [sigma_fudge=3]
#     reads P_fullgeo.txt, P_clashes.txt, P_cbetadev.log; temporary files fudgelists_tmp*

set outprefix   = ""
set sigma_fudge = 3
set modulo      = 9999
set cbsigma     = 0.032
set cbdof       = 3
foreach Arg ( $* )
  set key = `echo $Arg | awk -F= '{print $1}'`
  set val = `echo $Arg | awk -F= '{print $2}'`
  if("$Arg" !~ *=*) set outprefix = $Arg
  if("$key" == "sigma_fudge") set sigma_fudge = $val
end
if("$outprefix" == "") then
  echo "usage: fudgelists_runme.com P   (reads P_fullgeo.txt, P_clashes.txt, P_cbetadev.log)"
  exit 1
endif
foreach f ( ${outprefix}_fullgeo.txt ${outprefix}_clashes.txt ${outprefix}_cbetadev.log )
  if(! -e $f) then
    echo "ERROR: fudgelists_runme.com: no $f"
    exit 1
  endif
end

set t = ./fudgelists_tmp
ln -sf ${outprefix}_fullgeo.txt ${t}_fullgeo.txt

# this is in Angstrom
# v2: the deviation is a 3D distance: x = dev/cbsigma is a chi(cbdof) deviate; energy = square of the Gaussian deviate
#     with the same (two-sided) tail probability.  (v1: energy=(dev/0.05)^2, i.e. a signed 1D Gaussian deviate.)
cat ${outprefix}_cbetadev.log |\
 awk -v modulo=$modulo -v cbsigma=$cbsigma -v cbdof=$cbdof -F ":" '! /^SUMM|^pdb/ && NF>5{x=$6/cbsigma;\
   n=($5-1)%modulo+1;\
   if(cbdof==2){Q="exp(-"x"**2/2.0)"}else{Q="(erfc("x"/sqrt(2.0))+sqrt(2.0/pi)*"x"*exp(-"x"**2/2.0))"};\
   arg="(-invnorm("Q"/2.0))**2";\
   if(x>37)arg=x*x;\
   gsub("\"","",$0);\
   print "print "arg", \" "n" "$0"\""}' |\
gnuplot |&\
awk '{$1+=0;print "CBETADEV",$0}' |\
sort -k2gr >! ${t}_cbetadev.txt
# CBETADEV energy resnum%modulo otherstuff

# translate hydrogens to bonded parents
cat ${t}_fullgeo.txt |\
awk '! /^BOND/{next}\
   $8~/^H/{print "PARENT",$8,$10,"-",$14}\
   $14~/^H/{print "PARENT",$14,$16,"-",$8}' |\
sort -u >! ${t}_parents.txt

# make anti-clash restraints
egrep "^CLASH" ${outprefix}_clashes.txt |\
awk '{print;\
   nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"|")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1);\
   print "|"id1"|";\
   print "|"id2"|";\
   a1=substr(id1,13,5);a2=substr(id2,13,5);\
   f1=substr(id1,8,1);f2=substr(id2,8,1);\
   t1=substr(id1,9,4);t2=substr(id2,9,4);\
   c1=substr(id1,2,1);c2=substr(id2,2,1);\
   r1=substr(id1,3,5);r2=substr(id2,3,5);\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
   print "select",$3,a1,f1,t1,c1,r1,a2,f2,t2,c2,r2}' |\
cat ${t}_parents.txt - |\
awk '/^PARENT/{parent[$2,$3]=$5;next}\
  parent[$3,$5]{$3=parent[$3,$5]}\
  parent[$8,$10]{$8=parent[$8,$10]}\
  {print}' |\
egrep "^select" >! ${t}_selections.txt

cat ${t}_selections.txt ${t}_fullgeo.txt |\
awk '/^select/{del[$3,$4,$5,$6,$7]=del[$8,$9,$10,$11,$12]=$2;next}\
  ! /^NONBON/{next}\
  {nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"-")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1);\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   sym=w[6];sub(/^symop=/,"",sym);symm="";if(sym!="")symm=" symm y ! symop "sym;\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
  }\
   del[a1,f1,t1,c1,r1] && del[a2,f2,t2,c2,r2]{\
   print ;\
   ideal=$5+del[a1,f1,t1,c1,r1];sigma=0.2;\
   print "extern interval first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "dmin",ideal,"dmax 9999","smin",sigma,"smax 9999" symm;\
   print "#extern dist first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "value",ideal,"sigma",sigma,"type 0" symm}' |\
awk '{gsub("chain _","");print}' |\
egrep "^extern|^#extern" >! refmac_opts_declash.txt

# key ljenergy delta obs ideal 1 | atoms1 | atoms2
awk '/^NONBON/ && $2>1' ${t}_fullgeo.txt |\
cat ${t}_parents.txt - |\
awk '/^PARENT/{parent[$2,$3]=$5;next}\
  {print;\
   nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"-")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1);\
   print "|"id1"|";\
   print "|"id2"|";\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   sym=w[6];sub(/^symop=/,"",sym);if(sym=="")sym="x,y,z";\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
   energy=$2;delta=$3;ideal=$5;\
   if(parent[a1,t1]){a1=parent[a1,t1];delta+=0.1}\
   if(parent[a2,t2]){a2=parent[a2,t2];delta+=0.1}\
   print "select",energy,delta,ideal,a1,f1,t1,c1,r1,a2,f2,t2,c2,r2,sym}' |\
egrep "^select" >! ${t}_selections.txt
# parent atoms selected, symmetry operator last

# now look for distance between parent atoms in the non-bond list
# if its not there, make something up
cat ${t}_selections.txt ${t}_fullgeo.txt |\
awk '/^select/{\
     a1 =$5;f1 =$6;t1 =$7;c1 =$8;r1 =$9;\
     a2=$10;f2=$11;t2=$12;c2=$13;r2=$14;\
     key=a1" "f1" "c1" "r1" "a2" "f2" "c2" "r2" "$15;\
     oldE[key]=$2;\
     olddel[key]=$3;\
     oldideal[key]=$4;\
     ++selected[key];next}\
  ! /^NONBON/{next}\
  {nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"-")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1);\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   sym=w[6];sub(/^symop=/,"",sym);if(sym=="")sym="x,y,z";\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
   energy=$2;delta=$3;ideal=$5;\
   key=a1" "f1" "c1" "r1" "a2" "f2" "c2" "r2" "sym;\
   #print del[key]+0,selected[key],key;\
  }\
  selected[key]{\
   print ;\
   dd=del[key];\
   ideal+=dd;delta+=dd;\
   sigma=0.4/(delta+1.1);\
   if(sigma<0.001)sigma=0.001;\
   symm="";if(sym!="x,y,z")symm=" symm y ! symop "sym;\
   print "extern interval first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "dmin",ideal,"dmax 9999","smin",sigma,"smax 9999" symm;\
   print "#extern dist first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "value",ideal,"sigma",sigma,"type 0" symm;\
   ++printed[key]}\
  END{for(key in selected){\
   if(selected[key] && ! printed[key]){\
     split(key,w);\
     a1=w[1];f1=w[2];c1=w[3];r1=w[4];\
     a2=w[5];f2=w[6];c2=w[7];r2=w[8];\
     symm="";if(w[9]!="x,y,z")symm=" symm y ! symop "w[9];\
     sigma=0.4/(olddel[key]+1.1);\
     ideal=oldideal[key]+olddel[key]+0.7;\
     print selected[key],( ! selected[key] ),del[key],key;\
     print "extern interval first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
        "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
        "dmin",ideal,"dmax 9999","smin",sigma,"smax 9999" symm;\
     print "#extern dist first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
        "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
        "value",ideal,"sigma",sigma,"type 0" symm;\
   };\
  }}' |\
tee ${t}debump.txt |\
awk '{gsub("chain _","");print}' |\
egrep "^extern|^#extern" >! refmac_opts_debump.txt

# key energy delta obs ideal sigma | atoms
egrep "^BOND " ${t}_fullgeo.txt |\
awk -v s=$sigma_fudge '$2>s{print;\
   nxt=substr($0,index($0,"|")+1);\
   i=index(nxt,"-")-1;\
   id1=substr(nxt,1,i);\
   nxt=substr(nxt,i+2);\
   id2=substr(nxt,1,i);\
   print "|"id1"|";\
   print "|"id2"|";\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);\
   gsub(" ","_",c1);gsub(" ","_",c2);\
   ideal=$5;sigma=$6;\
   if(sigma<0.001)sigma=0.001;\
   if(a1~/^H/ || a2~/^H/){next}\
   print "extern dist first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "seco chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "value",ideal,"sigma",sigma,"type 0"}' |\
tee ${t}fixbond.txt |\
awk '{gsub("chain _","");print}' |\
egrep "^extern" >! refmac_opts_fixbond.txt

# key energy delta obs ideal sigma | atoms
egrep "^ANGLE " ${t}_fullgeo.txt |\
awk -v s=$sigma_fudge '$2>s{print;\
   i=index($0,"|");\
   nxt=substr($0,i+1);\
   i=index(nxt,"-");\
   id1=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id2=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id3=substr(nxt,1);\
   print "|"id1"|";\
   print "|"id2"|";\
   print "|"id3"|";\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   split(id3,w);a3=w[1];f3=w[2];t3=w[3];c3=w[4];r3=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);gsub(" ","_",f3);gsub(" ","_",f4);\
   gsub(" ","_",c1);gsub(" ","_",c2);gsub(" ","_",c3);gsub(" ","_",c4);\
   ideal=$5;sigma=$6;\
   print "extern angle first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "next chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "next chain",c3,"resi",r3,"atom",a3,"altecode",f3,\
      "value",ideal,"sigma",sigma}' |\
tee ${t}fixang.txt |\
awk '{gsub("chain _","");print}' |\
egrep "^extern" >! refmac_opts_fixang.txt

# add cbetadev - pseudo-dihedrals N-C-CA-CB and C-N-CA-CB (mirrors phenix c_beta_restraints)
# ideal values from phenix mmtbx/geometry_restraints/c_beta.py get_cb_target_angle_pair()
# CBETADEV energy rank pdbname :altcode:resname:chain:resnum:deviation:angle:...
egrep "^CBETA" ${t}_cbetadev.txt |\
awk -v s=$sigma_fudge '$2>s{\
   n=split($0,w,":");\
   f=w[2];res=tolower(w[3]);c=w[4];r=w[5];\
   gsub(" ","",f);gsub(" ","",c);gsub(" ","",r);gsub(" ","",res);\
   ncab=122.8;cnab=-122.6;\
   if(res=="ala"){ncab=122.9;cnab=-122.6}\
   if(res=="pro"){ncab=115.1;cnab=-120.7}\
   if(res=="val"||res=="thr"||res=="ile"){ncab=123.4;cnab=-122.0}\
   if(res=="gly"){ncab=121.6;cnab=-121.6}\
   print "extern torsion first chain",c,"resi",r,"atom N altecode",f,\
      "next chain",c,"resi",r,"atom C altecode",f,\
      "next chain",c,"resi",r,"atom CA altecode",f,\
      "next chain",c,"resi",r,"atom CB altecode",f,\
      "value",ncab,"sigma 2.5 period 0";\
   print "extern torsion first chain",c,"resi",r,"atom C altecode",f,\
      "next chain",c,"resi",r,"atom N altecode",f,\
      "next chain",c,"resi",r,"atom CA altecode",f,\
      "next chain",c,"resi",r,"atom CB altecode",f,\
      "value",cnab,"sigma 2.5 period 0";\
}' |\
tee ${t}fixCBdev.txt |\
awk '{gsub("chain _","");print}' |\
egrep "^extern" >! refmac_opts_fixCBdev.txt

# key energy delta obs ideal sigma | atoms
egrep "^TOR" ${t}_fullgeo.txt |\
awk -v s=$sigma_fudge '$2>s{print;\
   i=index($0,"|");\
   nxt=substr($0,i+1);\
   i=index(nxt,"-");\
   id1=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id2=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id3=substr(nxt,1,i-1);\
   nxt=substr(nxt,i+1);i=index(nxt,"-");\
   id4=substr(nxt,1);\
   print "|"id1"|";\
   print "|"id2"|";\
   print "|"id3"|";\
   print "|"id4"|";\
   split(id1,w);a1=w[1];f1=w[2];t1=w[3];c1=w[4];r1=w[5];\
   split(id2,w);a2=w[1];f2=w[2];t2=w[3];c2=w[4];r2=w[5];\
   split(id3,w);a3=w[1];f3=w[2];t3=w[3];c3=w[4];r3=w[5];\
   split(id4,w);a4=w[1];f4=w[2];t4=w[3];c4=w[4];r4=w[5];\
   gsub(" ","_",f1);gsub(" ","_",f2);gsub(" ","_",f3);gsub(" ","_",f4);\
   gsub(" ","_",c1);gsub(" ","_",c2);gsub(" ","_",c3);gsub(" ","_",c4);\
   ideal=$5;sigma=$6;\
   print "extern torsion first chain",c1,"resi",r1,"atom",a1,"altecode",f1,\
      "next chain",c2,"resi",r2,"atom",a2,"altecode",f2,\
      "next chain",c3,"resi",r3,"atom",a3,"altecode",f3,\
      "next chain",c4,"resi",r4,"atom",a4,"altecode",f4,\
      "value",ideal,"sigma",sigma,"period 0"}' |\
tee ${t}fixtorsion.txt |\
awk '{gsub("chain _","");print}' |\
egrep "^extern" >! refmac_opts_fixtorsion.txt

# planes need the .geo file: redo_harmonics takes no atoms from "extern plane" lines
cat /dev/null >! refmac_opts_fixplanes.txt
