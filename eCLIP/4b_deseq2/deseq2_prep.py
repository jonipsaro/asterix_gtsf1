#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed November 08, 2019

#This script will take the list of all_logical_genes from the files listed in
#prefixes, and make a composite list of all_logical_genes and the corresponding
#readcounts.  The output from this file is the input for DESeq2.

prefixes = ["S1", "S2", "S3", "S4", "S5", "S6"]

all_annotations = {}

#Go through each file and build a composite dictionary with all annotations
for prefix in prefixes:
    fin = open(prefix+"/"+prefix+"-all-logical_genes.txt", "r")
    for line in fin:
        values = line.split("\t")
        if values[0] not in all_annotations:
            all_annotations[values[0]] = [0, 0, 0, 0, 0, 0]
    fin.close()


#Go through each file and insert the raw read count for each annotation into its dictionary entry
i = 0
for prefix in prefixes:
    fin = open(prefix+"/"+prefix+"-all-logical_genes.txt", "r")
    for line in fin:
        values = line.split("\t")
        if values[0] not in all_annotations: print("ERROR - MISSING KEY")
        else: all_annotations[values[0]][i] = int(values[1])
    fin.close()
    i = i+1


#Go through each entry of all_annotations and output
fout = open("S123456-deseq2_input.txt", "w+")
for annotation in all_annotations:
    fout.write(annotation+"\t")
    fout.write(str(all_annotations[annotation][0])+"\t")
    fout.write(str(all_annotations[annotation][1])+"\t")
    fout.write(str(all_annotations[annotation][2])+"\t")
    fout.write(str(all_annotations[annotation][3])+"\t")
    fout.write(str(all_annotations[annotation][4])+"\t")
    fout.write(str(all_annotations[annotation][5])+"\n")

fout.close()

#Print out summary statistics
print("Total number of annotations: ",len(all_annotations))
