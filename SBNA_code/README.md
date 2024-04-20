# INSTRUCTIONS FOR RUNNING THE NETWORK ANALYSIS PIPELINE

E. coli TEM1 beta-lactamase (1BTL.pdb) will be used for the purposes of this tutorial.


## I. Initialize directory structure.

	1. Create a directory in your chosen location called "NetworkAnalysis".

	> mkdir NetworkAnalysis
	> cd NetworkAnalysis

	2. Populate this directory with the dependencies provided in this Git repo. You will need:
		- 1BTL.pdb (a sample PDB file - replace this with your own)
		- AAarea
		- allatoms_radii
		- atomicWeights
		- buildEdges.sh
		- calculateFinalScore.R
		- makeWaterBonds.py
		- modularity_analysis_centroid.R
		- modularity_analysis_energetics.R
		- pdb2edgeSC.py
		- pdb2polypeptide.py
		- terminalAtoms
		- uniqueAtoms

	3. For each structure you'd like to analyze, create a new subdirectory named for the PBD code, and include the relevant PDB file within that directory (here we've provided 1BTL by way of example). 

	> mkdir 1BTL
	> mv 1BTL.pdb 1BTL
	
## II. Install dependencies.

	1. Download PHENIX (Python-based Hierarchical ENvironment for Integrated Xtallography) for protonation.
		- Download binaries from https://www.phenix-online.org/download/ (after requesting a password via email).
		- Follow installation instructions and note the final location; if you are using a Mac, it will likely be installed in your Applications directory. 
		- Make sure that line 18 in the file buildEdges.sh reflects this path. 

	2. Download STRIDE for secondary structure assignment.
		- Download binaries from http://webclu.bio.wzw.tum.de/stride/install.html
		- Create a directory called "Stride" one level above your NetworkAnalysis folder and copy any Stride-related files into this directory.

## III. Build the network.

	1. Run the buildEdges.sh script from within the NetworkAnalysis folder.

	> sh buildEdges.sh 1BTL/1BTL.pdb

	2. You will now see a number of files within the 1BTL/ directory, which correspond to the intermediate and final network files. The final file to be analyzed is 1BTL_net.

## IV. Analyze the network.

There are 3 steps to analyzing the network: the energetics analysis, the centroid analysis, and the final analysis. Note that most of the analysis (energetics and centroid) is already conducted by buildEdges.sh. The final analysis needs to be run separately using the calculateFinalScore.R script. Run this from within the 1BTL/ directory, with the following usage:

	> cd 1BTL
    > R --vanilla --args 1BTL/1BTL/ < ../calculateFinalScore.R 

The output is named "FinalSum" and is a list of nodes (amino acids) with their scores as follows:

> head FinalSum

- GLY156	-1.12139304234615
- HIS158	-1.64504045031281
- ALA150	-1.9452762424746
- GLU147	-1.52075794536626
- LEU139	1.15367020700574
- LEU162	2.19848342502133
- ALA280	-0.961502115446914
- ILE279	-0.00840942137073508
- ARG164	1.77030331167432
- GLU171	-0.735560986652963

These are the final network scores.

**_(NB: The placement of hydrogens, as determined by the phenix.reduce algorithm is not entirely deterministic: reasons for this are described in Word, et al.(1999) "Asparagine and glutamine: using hydrogen atom contacts in the choice of sidechain amide orientation" J. Mol. Biol. 285, 1735-1747. This accounts for minor (not functionally relevant) discrepancies between the network scores calculated by individual runs of this pipeline.)_**

For questions related to the network analysis code, contact Elizabeth_Rossin@MEEI.HARVARD.EDU.
