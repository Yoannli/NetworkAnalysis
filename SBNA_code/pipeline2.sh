#!/bin/bash

## $1 is the pdb file (X.pdb)

out=${1%%.pdb}

echo "Cleaning PDB file"
source /home/yliy0004/ym65_scratch/yliy0004/phenix-1.20.1-4487/phenix_env.sh
phenix.pdbtools $1 remove="water"
cp $1 "$out"_original.pdb
mv $out"_modified.pdb" $out"_nowaters.pdb"
out=$out"_nowaters"

echo "Writing out polypeptide file..."
python /home/yliy0004/ym65_scratch/yliy0004/NetworkAnalysis/pdb2polypeptide.py $out".pdb"
awk '{print $1"\t"NR; print$2"\t"NR+1}' "$out"_polypeptidePairs | uniq > "$out"_AAorder

echo "Adding hydrogens..."
phenix.reduce -NOFLIP -Quiet $out".pdb" > $out"_H.pdb"
mv $out"_H.pdb" $out".pdb"

echo "Protonating waters..."
if grep -q "HOH" $out".pdb"
then
sed s/"AHOH"/" HOH"/g $out".pdb" | sed /"BHOH"/d > $out"_tmp.pdb"
phenix.ready_set add_h_to_water=True output_file_name=$out"_H" $out"_tmp.pdb"
rm $out"_tmp.pdb"
mv $out"_H.pdb" $out".pdb"
fi

echo "Making .phi file"
/home/yliy0004/ym65_scratch/yliy0004/Stride/stride $out".pdb" > $out".phi"
#fi

echo "Running full bond calculation script"
python /home/yliy0004/ym65_scratch/yliy0004/NetworkAnalysis/pdb2edgeSC.py $out".pdb"
#python pdb2edgeUNIQUE.py $out".pdb"

echo "Removing duplicates from BFactor file"
python /home/yliy0004/ym65_scratch/yliy0004/NetworkAnalysis/removeDuplicates-BFactor.py $out"_Bfactor"
mv $out"_Bfactor_duplicatesRemoved" $out"_Bfactor"

#echo "Finding water bonds"
#python /Users/rossin/Dropbox/Lab/SBNA/code/makeWaterBonds.py $out"_hb"
#awk '{split($1,x,"-");split($2,y,"-");if(x[1]x[2]!=y[1]y[2]){print}}' $out"_waterBonds" > $out"tmp"
#mv $out"tmp" $out"_waterBonds"
    
echo "Making "$out"_net"
awk '{print $1"\t"$2"\tMCMC\t10\tPP\tPP1\tPP2"}' $out"_polypeptidePairs" > $out"_net"
awk '{split($1,x,"-");if(s[x[1]"\t"x[2]]==""&&s[x[2]"\t"x[1]]==""){s[x[1]"\t"x[2]]=$2;next;}if(s[x[1]"\t"x[2]]!=""&&s[x[2]"\t"x[1]]==""){s[x[1]"\t"x[2]]=s[x[1]"\t"x[2]]+$2;next;}if(s[x[1]"\t"x[2]]==""&&s[x[2]"\t"x[1]]!=""){s[x[2]"\t"x[1]]=s[x[2]"\t"x[1]]+$2}}END{for(i in s){if(s[i]!=""&&s[i]>0){print i"\t"s[i]}}}' $out"_vdw" > $out"tmp"
awk '{print $1"-"$2;print $2"-"$1;}' $out"tmp" > $out"tmp2"
grep -f $out"tmp2" $out"_vdw" > $out"_vdw_noRepulse"
awk '{split($1,x,"-");if((x[1]~/PRO/&&x[5]=="CD")||(x[2]~/PRO/&&x[6]=="CD")){print x[1]"\t"x[2]"\t"x[3]x[4]"\t"$2+0"\tVDW\t"$3"\t"$4;}else{print x[1]"\t"x[2]"\t"x[3]x[4]"\t"$2"\tVDW\t"$3"\t"$4;}}' $out"_vdw_noRepulse" >> $out"_net"
awk '{print $1"\t"$2"\tSCSC\t"$3"\tPIPI\t"$4"\t"$5}' $out"_pipi2" >> $out"_net"
awk '{print $1"\t"$2"\tSCSC\t"$3"\tPICAT\t"$4"\t"$5}' $out"_pication2" >> $out"_net"
awk '{print $1"\t"$2"\tSCSC\t"$3"\tSS\t"$4"\t"$5}' $out"_disulf" >> $out"_net"
awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\tSCSC\t"$3"\tSB\t"$4"\t"$5}' $out"_saltBridges_Barlow" >> $out"_net"
awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\t"$3"\t"$10"\tHB\t"$11"\t"$12}' $out"_hb" | sed /HOH/d >> $out"_net"
#awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' $out"_waterBonds" >> $out"_net"
awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\t"$3"\t"3"\tMETAL\t"x[3]"\t"y[3]}' $out"_metal" >> $out"_net"
awk '{split($1,x,"-");split($2,y,"-");print x[1]x[2]"\t"y[1]y[2]"\tMC"$3"\t10\tDNA\tNT\t"x[3]}' $out"_DNA" >> $out"_net"


awk '{split($1,x,"-");if(s[x[1]"\t"x[2]]==""&&s[x[2]"\t"x[1]]==""){s[x[1]"\t"x[2]]=$2;next;}if(s[x[1]"\t"x[2]]!=""&&s[x[2]"\t"x[1]]==""){s[x[1]"\t"x[2]]=s[x[1]"\t"x[2]]+$2;next;}if(s[x[1]"\t"x[2]]==""&&s[x[2]"\t"x[1]]!=""){s[x[2]"\t"x[1]]=s[x[2]"\t"x[1]]+$2}}END{for(i in s){if(s[i]!=""&&s[i]>0){print i"\t"s[i]}}}' $out"_vdw2" > $out"tmp"
awk '{print $1"-"$2;print $2"-"$1;}' $out"tmp" > $out"tmp2"
grep -f $out"tmp2" $out"_vdw2" > $out"_vdw2_noRepulse"
awk '{split($1,x,"-");if((x[1]~/PRO/&&x[5]=="CD")||(x[2]~/PRO/&&x[6]=="CD")){print x[1]"\t"x[2]"\t"x[3]x[4]"\t"$2+0"\tVDW2\t"$3"\t"$4;}else{print x[1]"\t"x[2]"\t"x[3]x[4]"\t"$2"\tVDW2\t"$3"\t"$4;}}' $out"_vdw2_noRepulse" >> $out"_net"

echo "Removing any negative value edges"
awk '$4>0{print}' $out"_net" > $out"tmp"
mv $out"tmp" $out"_net"

echo "Creating QC file"
wc -l $out"_"* > $out".QC"

echo "Creating Centroid directory and moving files"
mkdir Centroid
cp "$out"_centroid* Centroid
cp $out".rsa" Centroid/$out"_centroid.rsa"
cp $out".phi" Centroid/$out"_centroid.phi"
mv $out"_centroidNetSC" Centroid/$out"_centroid_net"
mv $out"_centroidNetLigand" Centroid/$out"_centroid_ligand"
cp $out"_secondaryStructure" Centroid/$out"_centroid_secondaryStructure"
cp $out"_secondaryStructure2" Centroid/$out"_centroid_secondaryStructure2"
cp $out"_Bfactor" Centroid/$out"_centroid_Bfactor"

echo "Running energetic R script"
R --vanilla --args $out < /home/yliy0004/ym65_scratch/yliy0004/NetworkAnalysis/modularity_analysis_energetics.R

echo "Running centroid R script"
cd Centroid
R --vanilla --args $out"_centroid" < /home/yliy0004/ym65_scratch/yliy0004/NetworkAnalysis/modularity_analysis_centroid.R

