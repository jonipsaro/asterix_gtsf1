#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed November 08, 2019

#This script will convert "logical gene" enrichment to "real gene" enrichment.  Two input files are required.
#1. The output from DESeq2 that has been parsed (polished) to be tab-delimited (logical_gene, fold_change, p-value)
#2. The input file for DESeq2 that has raw read counts for each logical_gene.

#I.     The DESeq2 output may be filtered by p-value & fold_enrichment
#       Note that the user can set the maxium p-value and minimum enrichment output in the pval and enrich variables
#II.    The filtered DESeq2 output will be used to generate dictionaries of logical_genes and real_genes.
#---------
#III.   Reads per million for each logical_gene is calculated from the raw data
#---------
#IV.    The average_rpm (of the two replicates) is 1/n normalized based on the number of annotations in the logical_gene
#       Each logical_gene is now given the following entry:
#       logical_gene:   [fold_change, p-value, reads_per_million/n]
#---------
#V.     Now, the tricky part. Each real_gene's (overall) fold_change is calculated from all contributing logical_genes weighted by rpm/n.
#Each logical_gene that contains the real_gene annotation contributes (rpm/n)*fold_change normalized to the total rpm/n for that real_gene.
#Effectively, this means that the more "unique" a contributing logical_gene is to the real_gene, the more it dominantes overall fold-change.
#---------
#VI.    Output

#I.     Filter DESeq2 output by p-value
pval = 1.0          #Change maximum p-value cutoff here (1.0 includes everything)
enrich = 0.0        #Change minimum enrichment cutoff here (0.0 includes everything)
deseq2_output = open("S123456-deseq2_polished_formatted.txt", "r")
deseq2_filtered = open("S123456-deseq2_polished-en"+str(enrich)+"-p"+str(pval)+".txt", "w+")

output_count = 0
for line in deseq2_output:
    values = line.split("\t")
    if (float(values[2]) > pval): pass
    elif (float(values[1]) < enrich): pass
    else:
        deseq2_filtered.write(line)
        output_count += 1
        
deseq2_output.close()
deseq2_filtered.close()
print("Part 1 - Filter DESeq2 output by pvalue, enrchiment ("+str(pval)+", "+str(enrich)+"): Done")


#IIa.   Generate dictionary with logical_gene annotations as the keys and a list of lg_classes
deseq2_filtered = open("S123456-deseq2_polished-en"+str(enrich)+"-p"+str(pval)+".txt", "r")
logical_genes = {}      #This dictionary will eventually hold values of [fold_change, p-value, reads_per_milion/n]
lg_classes = []

for line in deseq2_filtered:
    values = line.split("\t")
    logical_gene = values[0]
    if logical_gene in logical_genes: pass
    else:
        logical_genes[logical_gene] = [float(values[1]), float(values[2].strip()), None]
        if (logical_gene.split(":")[0] not in lg_classes): lg_classes.append(logical_gene.split(":")[0])       
        
deseq2_filtered.close()
print("Part 2a - Generate dictionaries of logical_gene annotations: Done")


#IIb.   Split the logical_gene annotations by class
lgbc = {}   #lgbc logical_genes_by_class has the format {lg_class: {lg_in_class:[fold change, p-value, None]; lg2_in_class:[...]; ...}; lg_class2: {...}}
for lg_class in lg_classes:
    lgbc[lg_class] = {}

for logical_gene in logical_genes:
    lg_class = logical_gene.split(":")[0]
    lgbc[lg_class][logical_gene] = logical_genes[logical_gene] 
    

#IIc.   Generate dictionary with unique real_gene annotations as the keys
deseq2_filtered = open("S123456-deseq2_polished-en"+str(enrich)+"-p"+str(pval)+".txt", "r")
real_genes = {}

for line in deseq2_filtered:
    logical_gene = line.split("\t")[0]
    for real_gene in logical_gene.split(";"):
        if real_gene in real_genes: pass
        else: real_genes[real_gene] = [None]
        
deseq2_filtered.close()
print("Parts 2b,2c - Generate dictionary of real_gene annotations: Done")



#III.   Calculate reads_per_million for each library
#Calculate total reads in each library
S1_total = S2_total = S3_total = S4_total = S5_total = S6_total = 0
deseq2_input = open("S123456-deseq2_input_formatted.txt", "r")
for line in deseq2_input:
    values = line.split("\t")
    S1_total += float(values[1])
    S2_total += float(values[2])
    S3_total += float(values[3])
    S4_total += float(values[4])
    S5_total += float(values[5])
    S6_total += float(values[6])

S1m = S1_total/1000000
S2m = S2_total/1000000
S3m = S3_total/1000000
S4m = S4_total/1000000
S5m = S5_total/1000000
S6m = S6_total/1000000

deseq2_input.close()

#Calculate reads per million and output to file
deseq2_input = open("S123456-deseq2_input_formatted.txt", "r")
deseq2_input_normalized = open("S123456-deseq2_input_normalized.txt", "w+")

for line in deseq2_input:
    values = line.split("\t")
    deseq2_input_normalized.write(values[0]+"\t")
    deseq2_input_normalized.write(str(float(values[1])/S1m)+"\t")
    deseq2_input_normalized.write(str(float(values[2])/S2m)+"\t")
    deseq2_input_normalized.write(str(float(values[3])/S3m)+"\t")
    deseq2_input_normalized.write(str(float(values[4])/S4m)+"\t")
    deseq2_input_normalized.write(str(float(values[5])/S5m)+"\t")
    deseq2_input_normalized.write(str(float(values[6])/S6m)+"\n")

deseq2_input.close()
deseq2_input_normalized.close()
print("Part 3 - Calculate reads_per_million for each libary: Done")


#IVa.    Calculate the average rpm
#Make a dictionary with logical_genes (from deseq2_input) as keys and average_rpm as the value
#This assumes that the data to be averaged are in columns 5 and 6
deseq2_input_normalized = open("S123456-deseq2_input_normalized.txt", "r")
logical_gene_avg_rpm = {}

for line in deseq2_input_normalized:
    values = line.split("\t")
    if values[0] in logical_gene_avg_rpm: print("ERROR - Duplicate logical_gene. Please check the input file")
    else: logical_gene_avg_rpm[values[0]] = (float(values[5])+float(values[6]))/2

deseq2_input_normalized.close()
print("Part 4a - Calculate the avg_rpm for the IP samples (S5,S6): Done")


#IVb.   For each logical_gene in logical_genes, retrieve the avg_rpm and 1/n normalize it
for logical_gene in logical_genes:
    avg_rpm = logical_gene_avg_rpm[logical_gene]
    n = len(logical_gene.split(";"))
    
    logical_genes[logical_gene][2] = avg_rpm/n      #logical_genes now has the structure logical_gene: [fold_change, p-value, reads_per_million/n]
    lg_class = logical_gene.split(":")[0]
    lgbc[lg_class][logical_gene][2] = avg_rpm/n     #as do lgbc entries
    
print("Part 4b - Retrieve avg_rpm and 1/n normalize: Done")


#V.     Calculate the weighted, normalized, fold_change for each real_gene
print("Part 5 - Calculate the weighted, normalized, fold_change for each real_gene")
print("Number of real_genes:",len(real_genes))
print("This may take a while depending on the p-value cutoff")
real_gene_count = 0
for real_gene in real_genes:
    total_reads = 0
    total_score = 0
    rg_class = real_gene.split(":")[0]
    
    for logical_gene in lgbc[rg_class]:
        if real_gene in logical_gene.split(";"):
            total_score = total_score + (logical_genes[logical_gene][0]*logical_genes[logical_gene][2])
            total_reads = total_reads + logical_genes[logical_gene][2]

    real_genes[real_gene] = total_score/total_reads

    real_gene_count += 1
    if (real_gene_count % 10000 == 0): print("Real_genes processed so far =",real_gene_count)


#VI.    Output
fout = open("S123456-real_genes-en"+str(enrich)+"-p"+str(pval)+".txt", "w+")
for real_gene in real_genes:
    fout.write(real_gene+"\t"+str(real_genes[real_gene])+"\n")
fout.close()
