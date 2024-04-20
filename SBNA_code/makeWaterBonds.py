import sys
import itertools
import os

hbfile = open(sys.argv[1],'r').readlines()
outfile = open(sys.argv[1].strip("_hb")+"_waterBonds",'w')

HOH = {}
HOH2 = {}

for line in hbfile:
    line = line.split("\t")
    if "HOH" in line[0] and not "HOH" in line[1]:
        water = line[0].split("-")[0]+"-"+line[0].split("-")[1]
        if not water in HOH:
            HOH[water] = []
        HOH[water].append(line[1])
        HOH[water] = list(set(HOH[water]))
    if "HOH" in line[1] and not "HOH" in line[0]:
        water = line[1].split("-")[0]+"-"+line[1].split("-")[1]
        if not water in HOH:
            HOH[water] = []
        HOH[water].append(line[0])
        HOH[water] = list(set(HOH[water]))
    if "HOH" in line[0] and "HOH" in line[1]:
        water1 = line[0].split("-")[0]+"-"+line[0].split("-")[1]
        water2 = line[1].split("-")[0]+"-"+line[1].split("-")[1]
        if not water1 in HOH2:
            HOH2[water1] = []
        if not water2 in HOH2:
            HOH2[water2] = []
        HOH2[water1].append(water2)
        HOH2[water2].append(water1)
        HOH2[water1] = list(set(HOH2[water1]))
        HOH2[water2] = list(set(HOH2[water2]))

#Write out first order water interactions
for water in HOH:
    if len(HOH[water])>1:
        x = list(itertools.combinations(HOH[water],2))
        for item in x:
            outfile.write(item[0]+"\t"+item[1]+"\t"+water+"\n")

#Find second order interactions
for water1 in HOH:
    for atom1 in HOH[water1]:
        if water1 in HOH2:
            for water2 in HOH2[water1]:
                if water2 in HOH:
                    for atom2 in HOH[water2]:
                        if atom1 != atom2:
                            outfile.write(atom1+"\t"+atom2+"\t"+water1+","+water2+"\n")

outfile.close()

#Edge strength
waterfile = open(sys.argv[1].strip("_hb")+"_waterBonds",'r')
outfile = open(sys.argv[1].strip("_hb")+"_waterBonds2",'w')
singles = []

for line in waterfile:
    line = line.strip("\n").split()
    atom1 = line[0]; atom2 = line[1]
    if len(line[2].split(",")) == 1:
        water = "-".join(line[2].split("-")[0:2])
        os.system("grep -w "+atom1+" "+sys.argv[1]+" | grep "+water+" | uniq > "+sys.argv[1]+"tmp1")
        os.system("grep -w "+atom2+" "+sys.argv[1]+" | grep "+water+" | uniq > "+sys.argv[1]+"tmp2")
        tmp1 = open(sys.argv[1]+"tmp1",'r').readline().strip("\n").split()
        tmp2 = open(sys.argv[1]+"tmp2",'r').readline().strip("\n").split()
        weight1 = float(tmp1[9])
        weight2 = float(tmp2[9])
        if "SC" in tmp1[2]:
            type1 = "SC"
        else:
            type1 = "MC"
        if "SC" in tmp2[2]:
            type2 = "SC"
        else:
            type2 = "MC"
        pos1 = 10; pos2 = 10
        if "HOH" in tmp1[0]:
            pos1 = 11
        if "HOH" in tmp2[0]:
            pos2 = 11
        outfile.write(line[0]+"\t"+line[1]+"\t"+type1+type2+"\t"+str(sum([weight1,weight2])/2.0)+"\tWATER\t"+tmp1[pos1]+"\t"+tmp2[pos2]+"\t"+water+"\n")
        tmp = [atom1,atom2]
        tmp.sort()
        singles.append(tmp)
    else:
        water = line[2].split(",")
        water = ["-".join(water[0].split("-")[0:2]),"-".join(water[1].split("-")[0:2])]
        os.system("grep -w "+atom1+" "+sys.argv[1]+" | grep "+water[0]+" | uniq > "+sys.argv[1]+"tmp1")
        os.system("grep -w "+atom2+" "+sys.argv[1]+" | grep "+water[1]+" | uniq > "+sys.argv[1]+"tmp2")
        tmp1 = open(sys.argv[1]+"tmp1",'r').readline().strip("\n").split()
        tmp2 = open(sys.argv[1]+"tmp2",'r').readline().strip("\n").split()
        weight1 = float(tmp1[9])
        weight2 = float(tmp2[9])
        if "SC" in tmp1[2]:
            type1 = "SC"
        else:
            type1 = "MC"
        if "SC" in tmp2[2]:
            type2 = "SC"
        else:
            type2 = "MC"
        pos1 = 10; pos2 = 10
        if "HOH" in tmp1[0]:
            pos1 = 11
        if "HOH" in tmp2[0]:
            pos2 = 11
        tmp = [atom1,atom2]
        tmp.sort()
        if not tmp in singles:
            outfile.write(line[0]+"\t"+line[1]+"\t"+type1+type2+"\t"+str(sum([weight1,weight2])/2.0)+"\tWATER\t"+tmp1[pos1]+"\t"+tmp2[pos2]+"\t"+line[2]+"\n")



outfile.close()

os.system("mv "+sys.argv[1].strip("_hb")+"_waterBonds2 "+sys.argv[1].strip("_hb")+"_waterBonds")

