#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 17, 2019

#This script takes a .bed file (in command line input) and makes a list
#of all annotation classes within that .bed file.  Annotation classes should
#be in the 4th column and separated with ":" from the particular annotation.

import sys
arguments = sys.argv
file = sys.argv[1]
   
fin = open(file, "r")
anno_list = []

for line in fin:
    name = line.split("\t")[3].split(":")[0]
    if name not in anno_list: anno_list.append(name)

print(anno_list)
print(len(anno_list))
fin.close()
