#!/bin/bash
out=${1%%.pdb}
source /home/yliy0004/ym65_scratch/yliy0004/phenix-1.20.1-4487/phenix_env.sh
phenix.pdbtools $1 keep="$2"
mkdir $out"_monomer"
mv $out"_modified.pdb" $out"_monomer/"$out"_monomer.pdb"
mkdir $out"_multimer"
cp $1 $out"_multimer/"$out"_multimer.pdb"
cd $out"_monomer"
sh /home/yliy0004/ym65_scratch/yliy0004/NetworkAnalysis/pipeline2.sh $out"_monomer.pdb"
cd ../$out"_multimer"
sh /home/yliy0004/ym65_scratch/yliy0004/NetworkAnalysis/pipeline2.sh $out"_multimer.pdb"
cd ../$out"_monomer"
R --vanilla --args ../$out"_multimer" <  /home/yliy0004/ym65_scratch/yliy0004/NetworkAnalysis/calculateFinalScore.R
