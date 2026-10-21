#! /bin/tcsh -f
#
#
set maxadd = 1

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if( $assign ) then
      # re-set any existing variables
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
#      if("$key" == "outfile") set outfile = "$Val"
    else
      # no equal sign
#      if("$Arg" =~ *.pdb ) set pdbfile = $Arg
    endif
    if("$key" == "debug") set debug = "$Val"
end

if(! -e refmac_opts.txt) touch refmac_opts.txt


set file = refmac_opts_declash.txt
if( ! -s $file ) set file = refmac_opts_debump.txt


egrep "extern interv" $file |\
head -n $maxadd |\
awk '$26<0.5 && $9!~/^H/ && $18!~/^H/' |\
awk '/^extern/{$26="new";print}' |\
cat - refmac_opts.txt |\
awk '! /^extern/{next} $2!="interval"{next}\
 {c1=$5; r1=$7; a1=$9;f1=$11;\
  c2=$14;r2=$16;a2=$18;f2=$20;sc=26;\
  key1=c1" "r1" "a1" "f1;\
  key2=c2" "r2" "a2" "f2;\
  line[key1,key2]=$0;\
  sigma=$sc;}\
  sigma=="new"{\
    ++seen[key1,key2];\
    ++seen[key2,key1];\
    next;\
  }\
  seen[key1,key2]{\
     $sc/=1.5;\
     print $0;++printed[key1,key2];\
  }\
  END{for(keys in seen){\
    if(! printed[keys] && line[keys]!=""){print line[keys]; ++printed[keys]};\
  }}' |\
 awk 'BEGIN{sc=26} $sc=="new"{$sc="4"}\
   $sc<0.1{$sc=0.1}\
   NF>10{print}' |\
sort -k26g |\
 tee new_debump.txt 




cat refmac_opts_fixbond.txt |\
awk '$9!~/^H/ && $18!~/^H/ && !($9=="SG" && $18=="SG")' |\
head -n $maxadd |\
awk '/^extern/{$24="new";print}' |\
cat - refmac_opts.txt |\
awk '! /^extern/{next} $2!="dist"{next}\
 {c1=$5; r1=$7; a1=$9;f1=$11;\
  c2=$14;r2=$16;a2=$18;f2=$20;\
  key1=c1" "r1" "a1" "f1;\
  key2=c2" "r2" "a2" "f2;\
  line[key1,key2]=$0;\
  sigma=$24;}\
  sigma=="new"{\
    ++seen[key1,key2];\
    ++seen[key2,key1];\
    next;\
  }\
  seen[key1,key2]{\
     $24/=1.5;\
     print $0;++printed[key1,key2];\
  }\
  ! seen[key1,key2]{print $0;}\
  END{for(keys in seen){\
    if(seen[keys] && ! printed[keys] && line[keys]!=""){print line[keys]; ++printed[keys]};\
  }}' |\
 awk '$24=="new"{$24="0.1"}\
   $24<0.001{$24=0.001}\
   NF>10{print}' |\
sort -k24g |\
 tee new_bonds.txt 


cat refmac_opts_fixang.txt |\
awk '$9!~/^H/ && $18!~/^H/ && $27!~/^H/ && !($9=="CB" && $18=="SG" && $27=="SG") && !($9=="SG" && $18=="SG" && $27=="CB") && !($9=="CA" && $18=="CB" && $27=="SG")' |\
head -n $maxadd |\
awk '/extern angle/{$NF="new";print}' |\
cat - refmac_opts.txt |\
awk '! /^extern/{next} $2!="angle"{next}\
 {c1=$5; r1=$7; a1=$9;f1=$11;\
  c2=$14;r2=$16;a2=$18;f2=$20;\
  c3=$23;r3=$25;a3=$27;f3=$29;\
  key1=c1" "r1" "a1" "f1;\
  key2=c2" "r2" "a2" "f2;\
  key3=c3" "r3" "a3" "f3;\
  line[key1,key2,key3]=$0;\
  sigma=$NF;}\
  sigma=="new"{\
    ++seen[key1,key2,key3];\
    ++seen[key2,key1,key3];\
    next;\
  }\
  seen[key1,key2,key3]{\
     $NF/=1.5;\
     print $0;++printed[key1,key2,key3];\
  }\
  ! seen[key1,key2,key3]{print $0;}\
  END{for(keys in seen){\
    if(seen[keys] && ! printed[keys] && line[keys]!=""){print line[keys]; ++printed[keys]};\
  }}' |\
 awk '$NF=="new"{$NF="10"}\
   $NF<0.01{$NF=0.01}\
   NF>10{print}' |\
sort -k33g |\
 tee new_angles.txt 


cat refmac_opts_fixtorsion.txt |\
awk '$9!~/^H/ && $18!~/^H/ && $27!~/^H/ && $36!~/^H/ && !($9~/^[CN]$/ && $18~/^[CN]$/ && $27=="CA" && $36=="CB") && !($5=="B" && $7=="30" && $25=="30") && !($9=="CA" && $18=="CB" && $27=="SG" && $36=="SG") && !($5=="B" && $7=="12" && $9=="N" && $18=="CA" && $27=="CB" && $36=="SG")' |\
head -n $maxadd |\
awk '/extern tors/{$42=$42 "new";print}\' |\
cat - refmac_opts.txt |\
awk '! /^extern/{next} $2!="torsion"{next}\
 {c1=$5; r1=$7; a1=$9;f1=$11;\
  c2=$14;r2=$16;a2=$18;f2=$20;\
  c3=$23;r3=$25;a3=$27;f3=$29;\
  c4=$32;r4=$34;a4=$36;f4=$38;\
  key1=c1" "r1" "a1" "f1;\
  key2=c2" "r2" "a2" "f2;\
  key3=c3" "r3" "a3" "f3;\
  key4=c4" "r4" "a4" "f4;\
  line[key1,key2,key3,key4]=$0;\
  sigma=$42;}\
  sigma~/new/{\
    ++seen[key1,key2,key3,key4];\
    ++seen[key2,key1,key3,key4];\
    next;\
  }\
  seen[key1,key2,key3,key4]{\
     $42/=1.5;\
     print $0;++printed[key1,key2,key3,key4];\
  }\
  ! seen[key1,key2,key3,key4]{print $0;}\
  END{for(keys in seen){\
    if(seen[keys] && ! printed[keys] && line[keys]!=""){print line[keys]; ++printed[keys]};\
  }}' |\
 awk '$42~/new/{$42=$42+0}\
   $42<0.01{$42=0.01}\
   NF>10{print}' |\
sort -k42g |\
 tee new_tors.txt

awk '!($9~/^[CN]$/ && $18~/^[CN]$/ && $27=="CA" && $36=="CB")' new_tors.txt >! new_tors_noCB.txt
mv new_tors_noCB.txt new_tors.txt


head -n $maxadd refmac_opts_fixplanes.txt |\
awk '/extern plane/{$NF="new";print}' |\
cat - refmac_opts.txt |\
awk '! /^extern/{next} $2!="plane"{next}\
 {key="";for(i=5;i<NF;i+=9){\
    c=$i;r=$(i+2);a=$(i+4);f=$(i+6);\
    key=key"  "c" "r" "a" "f;\
  }\
  #print "DEBUG",key;\
  line[key]=$0;\
  sigma=$NF;}\
 sigma=="new"{\
    ++seen[key];\
    next;\
  }\
  seen[key]{\
     #print "SEEN" key;\
     $NF/=1.5;\
     print $0;++printed[key];\
  }\
  END{for(keys in seen){\
    if(! printed[keys] && line[keys]!=""){print line[keys]; ++printed[keys]};\
  }}' |\
 awk '$NF=="new"{$NF=0.02}\
   $NF<1e-7{$NF=1e-7}\
   NF>10{print}' |\
 tee new_planes.txt 


echo "updating refmac_opts.txt into refmac_opts_new.txt"
egrep -v "^extern interval" refmac_opts.txt >! refmac_opts_temp.txt
egrep "^extern interval" new_debump.txt | sort -k26g |\
cat >> refmac_opts_temp.txt
mv refmac_opts_temp.txt refmac_opts_new.txt

egrep -v "^extern dist" refmac_opts_new.txt >! refmac_opts_temp.txt
egrep "^extern dist" new_bonds.txt | sort -k24g |\
cat >> refmac_opts_temp.txt
mv refmac_opts_temp.txt refmac_opts_new.txt

egrep -v "^extern angl" refmac_opts_new.txt >! refmac_opts_temp.txt
egrep "^extern angle" new_angles.txt | sort -k33g |\
cat >> refmac_opts_temp.txt
mv refmac_opts_temp.txt refmac_opts_new.txt

egrep -v "^extern tors" refmac_opts_new.txt >! refmac_opts_temp.txt
egrep "^extern tors" new_tors.txt | sort -k42g |\
cat >> refmac_opts_temp.txt
mv refmac_opts_temp.txt refmac_opts_new.txt

egrep -v "^extern plane" refmac_opts_new.txt >! refmac_opts_temp.txt
egrep "^extern plane" new_planes.txt |\
cat >> refmac_opts_temp.txt
mv refmac_opts_temp.txt refmac_opts_new.txt




cat default.geo |\
awk '/nonbonded pdb=/{key="NONBOND";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;getline;\
     obs=$1;ideal=$2;sigma=1;energy=lj(obs,ideal);}\
   /bond pdb=/{key="BOND";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;getline;\
     ideal=$1;obs=$2;sigma=$4;energy=$6}\
   /angle pdb=/{key="ANGLE";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;split($0,w,"\"");id3=w[2];\
     getline;getline;\
     ideal=$1;obs=$2;sigma=$4;energy=$6}\
   /dihedral pdb=/{key="TORSION";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;split($0,w,"\"");id3=w[2];\
     getline;split($0,w,"\"");id4=w[2];\
     getline;getline;\
     ideal=$1;obs=$2;sigma=$5;energy=$7}\
   /chirality pdb=/{key="CHIR";split($0,w,"\"");id1=w[2];\
     getline;split($0,w,"\"");id2=w[2];\
     getline;split($0,w,"\"");id3=w[2];\
     getline;split($0,w,"\"");id4=w[2];\
     getline;getline;\
     ideal=$2;obs=$3;sigma=$5;energy=$7}\
   /plane pdb=/{key="PLANE";ideal=0;obs=$(NF-4);sigma=$(NF-3);nids="";maxE=0;\
     while(NF && ! /sigma/){split($0,w,"\"");id=w[2];\
       a=substr(id,1,4);f=substr(id,5,1);t=substr(id,6,4);c=substr(id,10,1);r=substr(id,11,5);\
       if(c==" ")c="_";if(f==" ")f="_";\
       nid=sprintf("%4s %s %4s %s %5s",a,f,t,c,r);\
       if(! /plane/){obs=$(NF-2);sigma=$(NF-1)};\
       energy=(obs/sigma)**2;\
       if(maxE<energy){maxE=energy;wobs=obs};\
       nids=nids" - "nid;\
       getline;\
     }\
     gsub("^ - ","",nids);\
     print key,maxE+0,sigma+0,"|",nids;\
   }\
   id1!=""{\
     a1=substr(id1,1,4);a2=substr(id2,1,4);a3=substr(id3,1,4);a4=substr(id4,1,4);\
     f1=substr(id1,5,1);f2=substr(id2,5,1);f3=substr(id3,5,1);f4=substr(id4,5,1);\
     t1=substr(id1,6,4);t2=substr(id2,6,4);t3=substr(id3,6,4);t4=substr(id4,6,4);\
     c1=substr(id1,10,1);c2=substr(id2,10,1);c3=substr(id3,10,1);c4=substr(id4,10,1);\
     r1=substr(id1,11,5);r2=substr(id2,11,5);r3=substr(id3,11,5);r4=substr(id4,11,5);\
     gsub(" ","_",f1);gsub(" ","_",f2);gsub(" ","_",f3);gsub(" ","_",f4);\
     gsub(" ","_",c1);gsub(" ","_",c2);gsub(" ","_",c3);gsub(" ","_",c4);\
     nid1=sprintf("%4s %s %4s %s %5s",a1,f1,t1,c1,r1);\
     nid2=sprintf("%4s %s %4s %s %5s",a2,f2,t2,c2,r2);\
     nid3=sprintf("%4s %s %4s %s %5s",a3,f3,t3,c3,r3);\
     nid4=sprintf("%4s %s %4s %s %5s",a4,f4,t4,c4,r4);\
     ids=nid1; if(nid2~/[0-9]/)ids=ids" - "nid2; if(nid3~/[0-9]/)ids=ids" - "nid3;if(nid4~/[0-9]/)ids=ids" - "nid4;\
     print key,ideal,sigma+0,"|",ids;\
     id1=id2=id3=id4=""} \
   function lj0(r,r0) {if(r==0)return 1e40;return 4*((r0*2**(-1./6)/r)**12-(r0*2**(-1./6)/r)**6)}\
 function lj(r,r0) {return lj0(r,r0)-lj0(6,r0)}' |\
cat >! default_geo.txt





 egrep "^extern" new_debump.txt |\
 awk -v deadband=1 '$2=="interval"{\
     c1=$5;r1=$7;a1=$9;f1=$11;\
     c2=$14;r2=$16;a2=$18;f2=$20;\
     v=$22;s=$26;\
#   print "    #",c1,c2,a1,a2,r1,r2,f1,f2,v,s;\
   print "    bond {"\
   print "      action = *add";\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   print "      distance_ideal =",v+deadband;\
   print "      sigma =",s;\
   print "      slack =",deadband;\
   print "    }";}' |\
 cat >! phenix_opts_debump.txt


 egrep "^extern" new_bonds.txt |\
 awk '$2=="dist"{\
     c1=$5;r1=$7;a1=$9;f1=$11;\
     c2=$14;r2=$16;a2=$18;f2=$20;\
     v=$22;s=$24;\
#   print "    #",c1,c2,a1,a2,r1,r2,f1,f2,v,s;\
   print "    bond {"\
   print "      action = *change";\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   print "      distance_ideal =",v;\
   print "      sigma =",s;\
   print "    }";}' |\
 cat >! phenix_opts_fixbonds.txt

egrep "^extern" new_angles.txt |\
awk '$2=="angle"{\
     c1=$5;r1=$7;a1=$9;f1=$11;\
     c2=$14;r2=$16;a2=$18;f2=$20;\
     c3=$23;r3=$25;a3=$27;f3=$29;\
     v=$31;s=$33;\
#   print "    #",c1,c2,c3,a1,a2,a3,r1,r2,r3,f1,f2,f3,v,s;\
   print "    angle {"\
   print "      action = *change";\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   print "      atom_selection_3 = \"name",a3,"and resseq",r3,"and chain",c3,"and altid",f3 "\"";\
   print "      angle_ideal =",v;\
   print "      sigma =",s;\
   print "    }";\
   }' |\
cat >! phenix_opts_fixangles.txt


cat default_geo.txt new_tors.txt |\
awk '/^TOR/{key="";for(i=5;i<=NF;i+=6){a=$i;f=$(i+1);c=$(i+3);r=$(i+4);\
   key=key" "a" "f" "c" "r;};\
   ++dflt[key]}\
  ! /^extern/{next}\
  {key="";for(i=5;i+9<NF;i+=9){c=$i;r=$(i+2);a=$(i+4);f=$(i+6);\
     key=key" "a" "f" "c" "r;}\
     action="add";\
     c1=$5;r1=$7;a1=$9;f1=$11;\
     c2=$14;r2=$16;a2=$18;f2=$20;\
     c3=$23;r3=$25;a3=$27;f3=$29;\
     c4=$32;r4=$34;a4=$36;f4=$38;\
     v=$40;s=$42;\
#   print "     #",dflt[key],key;\
     if(dflt[key]) action="change";\
#   print "    #",c1,c2,c3,c4,a1,a2,a3,a4,r1,r2,r3,r4,f1,f2,f3,f4,v,s;\
   print "    dihedral {"\
   print "      action = *" action;\
   print "      atom_selection_1 = \"name",a1,"and resseq",r1,"and chain",c1,"and altid",f1 "\"";\
   print "      atom_selection_2 = \"name",a2,"and resseq",r2,"and chain",c2,"and altid",f2 "\"";\
   print "      atom_selection_3 = \"name",a3,"and resseq",r3,"and chain",c3,"and altid",f3 "\"";\
   print "      atom_selection_4 = \"name",a4,"and resseq",r4,"and chain",c4,"and altid",f4 "\"";\
   print "      angle_ideal =",v;\
   print "      sigma =",s;\
   print "    }";\
   }' |\
cat >! phenix_opts_fixtorsion.txt


cat default_geo.txt new_planes.txt |\
awk '/^PLANE/{key="";for(i=5;i<=NF;i+=6){a=$i;f=$(i+1);c=$(i+3);r=$(i+4);\
   key=key" "a" "f" "c" "r;};\
   ++dflt[key]}\
   ! /^extern plane/{next}\
   {sigma=$NF;key="";\
     for(i=5;i<NF;i+=9){c=$i;r=$(i+2);a=$(i+4);f=$(i+6); key=key" "a" "f" "c" "r;}\
     if(dflt[key]) next;\
     print "    planarity {";\
     printf("      atom_selection = \"");\
     for(i=5;i<NF;i+=9){\
       c=$i;r=$(i+2);a=$(i+4);f=$(i+6);\
       printf("( name %s and resseq %s and chain %s and altid %s) ",a,r,c,f);\
       if(i<NF-9){printf("\n        or ")};\
     }\
     if(sigma+0 < 0.001) sigma = 0.001;\
     print "\"\n      sigma =",sigma;\
     print "    }";\
   }' |\
cat >! phenix_opts_fixplanes.txt

cat new_debump.txt new_bonds.txt new_angles.txt new_tors.txt new_planes.txt |\
cat >! refmac_opts_extern.txt

if(! -e refmac_opts_fixCBdev.txt) touch refmac_opts_fixCBdev.txt
set cbdev_src = refmac_opts_fixCBdev.txt
if(-e refmac_opts_allCBdev.txt) set cbdev_src = refmac_opts_allCBdev.txt
awk '/extern tors/{$42=$42 "new";print}' $cbdev_src |\
cat - refmac_opts.txt |\
awk '! /^extern/{next} $2!="torsion"{next}\
 {c1=$5;r1=$7;a1=$9;f1=$11;\
  c2=$14;r2=$16;a2=$18;f2=$20;\
  c3=$23;r3=$25;a3=$27;f3=$29;\
  c4=$32;r4=$34;a4=$36;f4=$38;\
  key=c1" "r1" "a1" "f1" "c2" "r2" "a2" "f2" "c3" "r3" "a3" "f3" "c4" "r4" "a4" "f4;\
  line[key]=$0;\
  sigma=$42;}\
  sigma~/new/{\
    ++seen[key];\
    next;\
  }\
  seen[key]{\
    $42/=1.5;\
    print $0;++printed[key];\
  }\
  END{for(k in seen){\
    if(printed[k]==0 && line[k]!=""){print line[k]}\
  }}' |\
awk '$42~/new/{$42=$42+0}\
  $42<2.0{$42=2.0}\
  NF>10 && $27=="CA" && $36=="CB"{print}' |\
tee new_cbdev_tors.txt |\
awk '$2=="torsion"{\
  printf "    dihedral {\n";\
  printf "      action = *add\n";\
  printf "      atom_selection_1 = \"chain %s and resseq %s and name %s and altid %s\"\n",$5,$7,$9,$11;\
  printf "      atom_selection_2 = \"chain %s and resseq %s and name %s and altid %s\"\n",$14,$16,$18,$20;\
  printf "      atom_selection_3 = \"chain %s and resseq %s and name %s and altid %s\"\n",$23,$25,$27,$29;\
  printf "      atom_selection_4 = \"chain %s and resseq %s and name %s and altid %s\"\n",$32,$34,$36,$38;\
  printf "      angle_ideal = %s\n",$40;\
  printf "      sigma = %s\n",$42;\
  printf "    }\n";\
}' |\
cat >! phenix_opts_fixCBdev.txt
egrep "^extern tors" new_cbdev_tors.txt >> refmac_opts_new.txt

cat phenix_opts_debump.txt phenix_opts_fixbonds.txt phenix_opts_fixangles.txt \
  phenix_opts_fixtorsion.txt phenix_opts_fixplanes.txt phenix_opts_fixCBdev.txt |\
awk 'BEGIN{print "  geometry_restraints.edits {"} {print} END{print "  }"}' |\
awk 'BEGIN{print "refinement {"} {print} END{print "}"}' |\
cat >! phenix_fudge.eff






