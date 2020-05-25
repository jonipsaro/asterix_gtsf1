#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed November 08, 2019

#This script will take tab-delimited text files in the format
#chr start stop read_id-count score strand[+/-]  <-- [Columns 0-5] Information for the mapped read
#chr start stop annotation score strand[+/-] <-- [Columns 6-11] Information for the corresponding annotation

#And determine 'logical_genes' from multimappers (*-logical-genes.txt)
#The logical genes will be output as geneA;geneB;... followed by the read count (tab-delimited)

#Additionally, logical_gene will be associated with its corresponding reads (*-logical-gene-reads.txt)
#for later lookup.

#The user should specify the sample names in the prefixes variable (list of strings)
#The user should specify the gene classes in the suffixes variable (list of strings)

prefixes = ["S1", "S2", "S3", "S4", "S5", "S6"]

for prefix in prefixes:
    logfile = open(prefix+"/"+prefix+"-logical_genes.log", "w+")
    logfile.write("Summary for "+prefix+":\n")
    suffixes = ["rRNA", "snoRNA", "scaRNA", "snRNA", "ribozyme", "lncRNA", "tRNA", "pseudogene_tRNA", "mttRNA", "miRNA", "pri-miRNA", "exon", "TE_DNA", "TE_SINE", "TE_LINE", "TE_Satellite", "TE_LTR", "TE_RNA", "TE_RC", "TE_Other", "TE_Unknown", "intron", "pre-mRNA", "piRNA_cluster", "sncRNA", "miscRNA", "pseudogene", "asRNA", "no_annotation"]

    for suffix in suffixes:
        fin = open(prefix+"/"+prefix+"-"+suffix+".txt", "r")

        #Make a dictionary with read_ids as the key and list of annotations for the value
        read_annotations = {}
        count_duplicates = 0
        duplicate_annotations = []
        
        for line in fin:
            values = line.split("\t")
            if (values[3] not in read_annotations): read_annotations[values[3]] = []
            if (values[9] in read_annotations[values[3]]):
                count_duplicates += 1
                if (values[9] not in duplicate_annotations): duplicate_annotations.append(values[9])
                
            read_annotations[values[3]].append(values[9])

        logfile.write("\tUnique read ids - "+suffix+": "+str(len(read_annotations))+"\n")
        logfile.write("\tIdentical annotations that are multimappers - "+suffix+": "+str(count_duplicates)+"\n")
        logfile.write("\tList of annotations that map to more than one distinct locus: "+str(duplicate_annotations)+"\n")
        
        fin.close()

        #Go through the annotations dictionary, sort each annotation list, and make two new dictionaries:
        #1. With logical_genes as the key and read_count as the entry
        #2. With logical_genes as the key and a list of read_ids as the entry
        logical_gene_counts = {}
        logical_gene_read_ids = {}
        for read_id in read_annotations:
            read_annotations[read_id].sort()
            logical_gene = ";".join(read_annotations[read_id])
            if (logical_gene not in logical_gene_counts):
                logical_gene_counts[logical_gene] = 0
                logical_gene_read_ids[logical_gene] = []
            logical_gene_counts[logical_gene] += int(read_id.split("-")[1])
            logical_gene_read_ids[logical_gene].append(prefix+"_"+read_id)

        #Go through the logical_gene_count dictionary and output
        fout = open(prefix+"/"+prefix+"-"+suffix+"-logical_genes.txt", "w+")
        fout2 = open(prefix+"/"+prefix+"-"+suffix+"-logical_genes_to_readIDs.txt", "w+")
        ordered_logical_genes = sorted(logical_gene_counts.keys())
        for logical_gene in ordered_logical_genes:
            fout.write(logical_gene+"\t"+str(logical_gene_counts[logical_gene])+"\n")
            fout2.write(logical_gene+"\t"+str(logical_gene_read_ids[logical_gene])+"\n")
        logfile.write("\tNumber of logical genes - "+suffix+": "+str(len(logical_gene_counts))+"\n\n")
        
        fout.close()
        fout2.close()
        
    logfile.close()
