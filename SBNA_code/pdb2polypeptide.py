import sys
import string

file = open(sys.argv[1],'r').readlines()
out = open(sys.argv[1].replace(".pdb","")+"_polypeptidePairs",'w')
A=""
B=""
Achain = ""
Bchain = ""

def isAminoAcid(acid):
    if acid in ["ARG","HIS","LYS","ASP","GLU","SER","THR","ASN","GLN","CYS","PRO","GLY","ALA","VAL","ILE","LEU","MET","PHE","TYR","TRP"]:
        return(True)
    else:
        return(False)

for line in file:
    if line[0:4]=="ATOM" and A!="":
        if line[26] in string.ascii_uppercase:
            continue
        residuename = line[17:20];residuename=residuename.lstrip().strip()
        chainname = line[21];chainname=chainname.lstrip().strip()
        residueposition = line[22:26]; residueposition=residueposition.lstrip().strip()
        if A==residuename+residueposition:
            continue
        if isAminoAcid(residuename) and residueposition.isdigit():
            B = residuename+residueposition
            Bchain = chainname
            if Achain==Bchain and float(B[3:])==(float(A[3:])+1):
                out.write(A+Achain+"\t"+B+Bchain+"\tSCORE\n")
    if line[0:4]=="ATOM":
        residuename = line[17:20];residuename=residuename.lstrip().strip()
        chainname = line[21];chainname=chainname.lstrip().strip()
        residueposition = line[22:26]; residueposition=residueposition.lstrip().strip()
        if isAminoAcid(residuename) and residueposition.isdigit():
            A = residuename+residueposition
            Achain = chainname

out.close()
        
