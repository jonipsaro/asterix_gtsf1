#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed November 08, 2019

#This script will take files of normalized tRNA reads & counts and output
#characteristics of those reads as they map onto a "model tRNA" of a specified length.

#The user must set the total_library_reads variable below (run get_total_library_reads.py to
#generate a list of ints).
#The length of the model tRNA can be set below in the model_tRNA_length variable.

#The following files are output:
	#A histogram of 5' read ends according to their position on the model tRNA (5_prime)
	#A histogram of the 3' read ends according to their positon on the model tRNA (3_prime)
	#A pileup of the span of positions covered on the model tRNA (span)
	#A histogram of read lengths (length)

#The following parameters must be updated to suit the particulars of the data:
	#The length of the model tRNA (model_tRNA_length)
	#Library sizes (run get_total_library_reads.py and copy the output into the total_library_reads list below)

#The following parameters may be updated:
	#One or more position in the model tRNA that you wish to extract bed coordinates for
	#	[get_sites] (for example, if you would like the bed coordinates for enriched 
	#	5' sites that correspond to position 22 in the tRNA)

#Setup
import math
import pandas as pd
import os

if (os.path.exists("tRNA_analysis") == False): os.mkdir("tRNA_analysis")

#Initialize variables
#Indicate the file list for the files to be processed
file_list = ["S1", "S2", "S3", "S4", "S5", "S6"]

#Indicate the size of the "model" tRNA to be used (i.e. all tRNA annotations will be scaled to this size)
model_tRNA_length = 73

#Open the bed file for optional output of particular sites
#To skip bed output, set get_sites = []
get_sites = [22]
if (len(get_sites) != 0): bed_out = open("tRNA_analysis/"+"-".join(map(str, get_sites))+"nt_coord.bed", "w+")

#These numbers should be generated using get_total_library_reads.py (list of ints)
total_library_reads = 

#Set up data frames
positions = list(range(1,model_tRNA_length+1))
all_5prime = pd.DataFrame(index = positions)
all_3prime = pd.DataFrame(index = positions)
all_span = pd.DataFrame(index = positions)
all_length = pd.DataFrame(index = positions)


#Open each file of tRNA reads and determine the sum of the 1/n-normalized rpm at each position
#This is done for the 5' end of each read, the 3' end of each read, and the span of each read
j = 0
for file_name in file_list:
	fin = open(file_name+"/"+file_name+"-tRNA_normalized.txt", "r")

	#Initialize or re-initialize totals for each file
	total_counts_5 = {}
	total_counts_3 = {}
	total_counts_span = {}
	total_counts_length = {}

	for each in positions: 
		total_counts_5[each] = 0
		total_counts_3[each] = 0
		total_counts_span[each] = 0
		total_counts_length[each] = 0

	#Extract values from the bed files
	for line in fin:
		values = line.split("\t")
		read_start = int(values[1])
		read_end = int(values[2])
		strand = values[5]
		anno_start = int(values[7])
		anno_end = int(values[8])
		norm_count = float(values[12])
		norm_count_rpm = (norm_count/total_library_reads[j])*1000000

		anno_length = anno_end-anno_start

		#Calculate 5' and 3' sites and the total span
		if strand == "+":
			_5_prime_dist = read_start-anno_start
			_3_prime_dist = anno_end-read_end
		elif strand == "-":
			_5_prime_dist = anno_end-read_end
			_3_prime_dist = read_start-anno_start

		_5_prime_site = math.ceil(float((_5_prime_dist+1)/anno_length)*model_tRNA_length)
		_3_prime_site = math.ceil(float(model_tRNA_length-((_3_prime_dist/anno_length)*model_tRNA_length)))
		span = list(range(_5_prime_site, _3_prime_site))
		length = _3_prime_site-_5_prime_site+1

		#Add the counts to the running totals
		total_counts_5[_5_prime_site] = total_counts_5[_5_prime_site]+norm_count_rpm
		total_counts_3[_3_prime_site] = total_counts_3[_3_prime_site]+norm_count_rpm
		total_counts_length[length] = total_counts_length[length]+norm_count_rpm

		for each in span:
			total_counts_span[each] = total_counts_span[each]+norm_count_rpm

		#Optional: Extract bed entries for particular 5' sites in the pulldowns
		#Note: Setting j >= 4 means that only files S5 and S6 will be output
		if (j >= 4) and (_5_prime_site in get_sites) and (len(get_sites) != 0):
			bed_out.write(line.strip()+"\t"+str(norm_count_rpm)+"\n")

	#Add entries to the appropriate dataframes
	_5_df = pd.DataFrame.from_dict(total_counts_5, orient='index')
	_3_df = pd.DataFrame.from_dict(total_counts_3, orient='index')
	_span_df = pd.DataFrame.from_dict(total_counts_span, orient='index')
	_length_df = pd.DataFrame.from_dict(total_counts_length, orient='index')

	_5_df.columns = [file_name]
	_3_df.columns = [file_name]
	_span_df.columns = [file_name]
	_length_df.columns = [file_name]

	all_5prime = pd.merge(all_5prime, _5_df, left_index=True, right_index=True)
	all_3prime = pd.merge(all_3prime, _3_df, left_index=True, right_index=True)
	all_span = pd.merge(all_span, _span_df, left_index=True, right_index=True)
	all_length = pd.merge(all_length, _length_df, left_index=True, right_index=True)

	j += 1
	fin.close()

#Add final formatting to dataframes
all_5prime = all_5prime.round(decimals=6)
all_3prime = all_3prime.round(decimals=6)
all_span = all_span.round(decimals=6)
all_length = all_length.round(decimals=6)

#Output
all_5prime.to_csv(r'tRNA_analysis/tRNA_5prime_positions.txt')
all_3prime.to_csv(r'tRNA_analysis/tRNA_3prime_positions.txt')
all_span.to_csv(r'tRNA_analysis/tRNA_read_span.txt')
all_length.to_csv(r'tRNA_analysis/tRNA_read_length.txt')

#Manage other files
if (len(get_sites) != 0): bed_out.close()