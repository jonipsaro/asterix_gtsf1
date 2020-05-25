#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 18, 2019

#This script will take output from normalize_multimappers.py
#And generate a table of reads in each annotation class in each sample

file_names = ["S1", "S2", "S3", "S4", "S5", "S6"]

samples = []
sample_sums = []
anno_types = []
total_counts = {}       #This dictionary has anno_type as the key and
                        #[S1, S2, ...] read counts as the entry

i = 0
for each in file_names:
    fin = open(each+"/"+each+"_normalize.log", "r")
    sample_sum = 0
    
    for line in fin:
        sample = line.split()[0]
        if sample not in samples: samples.append(sample)
        
        anno_type = line.split()[1].split(":")[0]
        if anno_type not in anno_types: anno_types.append(anno_type)

        anno_counts = int(line.split()[2].strip())

        if anno_type not in total_counts: total_counts[anno_type] = [None, None, None, None, None, None]
        total_counts[anno_type][i] = anno_counts

        sample_sum = sample_sum+anno_counts

    sample_sums.append(sample_sum)
    i += 1
    fin.close()


fout = open("total_annotations.log", "w+")
fout.write("annotation"+"\t"+"\t".join(samples)+"\n")

for each in total_counts:
    total_counts_out = map(str, total_counts[each])
    fout.write(each+"\t"+"\t".join(total_counts_out)+"\n")
    
sample_sums_out = map(str, sample_sums)
fout.write("total"+"\t"+"\t".join(sample_sums_out)+"\n")   #If you would like to include a total at the end, uncomment this line

fout.close() 
