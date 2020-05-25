#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 17, 2019

#This script will simply take the no_annotation.txt file generated after filtering
#and reformat the "no_annotation" entries so that they match the formatting of those
#with annotations

import sys
arguments = sys.argv
if len(sys.argv) == 1:
   input("By default, this script will run samples [\"S1\", \"S2\", \"S3\", \"S4\", \"S5\", \"S6\"].\nTo run a single file, include it in the argument.\n\nPress Enter to continue or Ctrl+C to break")
   prefixes = ["S1", "S2", "S3", "S4", "S5", "S6"]
else:
   prefixes = sys.argv[1:]

for prefix in prefixes:
    fin = open(prefix+"_no_annotation.txt", "r")
    fout = open(prefix+"_no_annotation_formatted.txt", "w+")

    for line in fin:
        values = line.split("\t")
        fout.write(line.rstrip()+"\t"+values[0]+"\t"+values[1]+"\t"+values[2]+"\tno_annotation\t"+values[4]+"\t"+values[5])

    fin.close()
    fout.close()
    
