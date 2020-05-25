#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 17, 2019

#This script will take a .bed file that has overlapping aligments from the same read
#(i.e. different numbers of mismatches or insertions) and output only the smallest alignment range.
#This is necessary to avoid "multi-aligners" within a gene that are not true "multi-mappers"

#Example:
#chr1	169270417	169270455	3544811-1	0	-
#chr1	169270417	169270533	3544811-1	0	-
#chr1	169270417	169282638	3544811-1	0	-
#chr13	9833632	        9833670	        3544811-1	0	-
#chr13	9833632	        9833749	        3544811-1	0	-
#chr16	11144239	11144356	3544811-1	0	+
#chr16	11144318	11144356	3544811-1	0	+
#chr9	56071439	56071477	3544811-1	0	-
#chr9	56071439	56071556	3544811-1	0	-

#Result:
#chr1	169270417	169270455	3544811-1	0	-
#chr13	9833632	        9833670	        3544811-1	0	-
#chr16	11144318	11144356	3544811-1	0	+
#chr9	56071439	56071477	3544811-1	0	-

import os
import sys
arguments = sys.argv
if len(sys.argv) == 1:
   input("By default, this script will run samples [\"S1\", \"S2\", \"S3\", \"S4\", \"S5\", \"S6\"].\nTo run a single file, include it in the argument.\n\nPress Enter to continue or Ctrl+C to break")
   prefixes = ["S1", "S2", "S3", "S4", "S5", "S6"]
else:
   prefixes = sys.argv[1:]

for prefix in prefixes:
    edit_left = edit_right = edit_overlap = 0

    #Left-bound alignment
    #Make a dictionary with chr_start_readid as keys and [chr,start,stop,read_id,score,strand] as values
    bed = {}
    fin = open(prefix+"_unique.bed", "r")

    for line in fin:
        values = line.split("\t")
        key = values[0]+"_"+values[1]+"_"+values[3]
        if (key not in bed): bed[key] = [values[0], values[1], values[2], values[3], values[4], values[5]]
        else:
            bed[key][2] = str(min(int(bed[key][2]),int(values[2])))
            edit_left += 1

    fin.close()
    print("Number of left-bound tidying actions:", edit_left)

    #Write the dictionary to a temporary file
    fout = open(prefix+"_unique_left-aligned_temp.txt", "w+")
    for key in bed:
        fout.write(bed[key][0]+"\t"+bed[key][1]+"\t"+bed[key][2]+"\t"+bed[key][3]+"\t"+bed[key][4]+"\t"+bed[key][5])
    fout.close()


    #Right-bound alignment
    #Make a dictionary with chr_end_readid as keys and [chr,start,stop,read_id,score,strand] as values
    bed = {}
    fin = open(prefix+"_unique_left-aligned_temp.txt", "r")

    for line in fin:
        values = line.split("\t")
        key = values[0]+"_"+values[2]+"_"+values[3]
        if (key not in bed): bed[key] = [values[0], values[1], values[2], values[3], values[4], values[5]]
        else:
            bed[key][1] = str(max(int(bed[key][1]),int(values[1])))
            edit_right += 1
            
    fin.close()
    print("Number of right-bound tidying actions:", edit_right)

    #Write the dictionary to a temporary file
    fout = open(prefix+"_unique_right-aligned_temp.txt", "w+")
    for key in bed:
        fout.write(bed[key][0]+"\t"+bed[key][1]+"\t"+bed[key][2]+"\t"+bed[key][3]+"\t"+bed[key][4]+"\t"+bed[key][5])
    fout.close()


    #Remove overlapping alignments with the same read_id
    #Make a dictionary with chr_readid as keys and a list of [[chr,start,stop,read_id,score,strand]] as values
    bed = {}
    fin = open(prefix+"_unique_right-aligned_temp.txt", "r")

    for line in fin:
        values = line.split("\t")
        key = values[0]+"_"+values[3]

        if (key not in bed):
            bed[key] = []
            bed[key].append([values[0], values[1], values[2], values[3], values[4], values[5]])
        #Check to see if the incoming alignment overlaps an existing one with the same read_id
        #If so, keep the left-most coordinates.  If not, add the alignment to the dictionary.
        else:
            distinct = True
            for i in range(len(bed[key])):
                if (bed[key][i][1] <= values[1]) and (bed[key][i][2] > values[1]):
                    edit_overlap += 1
                    distinct = False
                elif (bed[key][i][1] >= values[1]) and (bed[key][i][1] < values[2]):
                    bed[key][i][1] = values[1]
                    bed[key][i][2] = values[2]
                    edit_overlap += 1
                    distinct = False
            if (distinct): bed[key].append([values[0], values[1], values[2], values[3], values[4], values[5]])

    fin.close()
    print("Number of overlap removal tidying actions:", edit_overlap)


    #Write out the final dictionary
    fout = open(prefix+"_unique_tidy.bed", "w+")
    for key in bed:
        for i in range(len(bed[key])):
            fout.write(bed[key][i][0]+"\t"+bed[key][i][1]+"\t"+bed[key][i][2]+"\t"+bed[key][i][3]+"\t"+bed[key][i][4]+"\t"+bed[key][i][5])
    fout.close()

    print("Number of tidying actions total: ", edit_left+edit_right+edit_overlap)


    os.remove(prefix+"_unique_left-aligned_temp.txt")
    os.remove(prefix+"_unique_right-aligned_temp.txt")
