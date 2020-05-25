#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 18, 2019

#This script will simply take multiple pileup files (annotation, read_count)
#and combine the entries into one file based on the annotation name.

#The script will process both the raw pileup files and the 
#normalized ones (annoation, read_count, normalized_count)

import math
file_list = ["S1", "S2", "S3", "S4", "S5", "S6"]		#Define the file list here


#Part I. Raw Pileup File Combination
i = 0
count_table = {}

for file in file_list:									#Open each file
	file_name = file+"-annotation_pileup_raw.txt"
	fin = open(file_name, "r")

	for line in fin:									#Make a dictionary for the annotations
		annotation = line.split("\t")[0]
		count = line.split("\t")[1].strip()

		#Add count value to each list
		if annotation not in count_table: count_table[annotation] = [0, 0, 0, 0, 0, 0]
		count_table[annotation][i] = int(count)

	fin.close()
	i += 1


fout = open("S123456-annotation_pileup_raw.txt", "w+")	#Write the output
for each in count_table:
	count_out = "\t".join(map(str, count_table[each]))
	fout.write(each+"\t"+count_out+"\n")
fout.close()

#Part II. Normalized & Rounded Pileup File Combination
i = 0
count_table_norm = {}

for file in file_list:									#Open each file
	file_name = file+"-annotation_pileup_normalized.txt"
	fin = open(file_name, "r")

	for line in fin:									#Make a dictionary for the annotations
		annotation = line.split("\t")[0]
		norm_count = line.split("\t")[2].strip()


		#Add count value to each list, round up to nearest whole number
		if annotation not in count_table_norm: count_table_norm[annotation] = [0, 0, 0, 0, 0, 0]
		count_table_norm[annotation][i] = int(math.ceil(float(norm_count)))

	fin.close()
	i += 1


fout = open("S123456-annotation_pileup_normalized_rounded.txt", "w+")	#Write the output
for each in count_table_norm:
	count_out = "\t".join(map(str, count_table_norm[each]))
	fout.write(each+"\t"+count_out+"\n")
fout.close()
