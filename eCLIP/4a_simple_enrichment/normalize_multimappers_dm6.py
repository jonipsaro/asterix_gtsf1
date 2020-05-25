#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 18, 2019

#This script will normalize multimapping read counts within each class
#The input files are the class-filtered files (bed+bed files).
#The output files will have one additional column which is the normalized read_count

#Note:  The suffixes list must correspond to the list of gene classes you wish to normalize
#This can/should be the same list from the get_annotation_classes script

import sys
arguments = sys.argv
if len(sys.argv) == 1:
   input("By default, this script will run samples [\"S1\", \"S2\", \"S3\", \"S4\", \"S5\", \"S6\"].\nTo run a single file, include it in the argument.\n\nPress Enter to continue or Ctrl+C to break")
   prefixes = ["S1", "S2", "S3", "S4", "S5", "S6"]
else:
   prefixes = sys.argv[1:]
   
for prefix in prefixes:
   logfile = open(prefix+"/"+prefix+"_normalize.log", "w+")

   suffixes = ["rRNA", "intron_28SrRNA-Psi", "pseudogene_5SrRNA-Psi", "pseudogene_5.8SrRNA-Psi", "pseudogene_28SrRNA-Psi", "pseudogene_2SrRNA-Psi", "pseudogene_18SrRNA-Psi", "snoRNA", "snRNA", "structRNA_7SLRNA", "structRNA_RNaseP", "structRNA_RNaseMRP", "structRNA_sbRNA", "structRNA_scaRNA", "intron_lncRNA", "lncRNA", "intron_tRNA", "tRNA", "pseudogene_tRNA", "mttRNA", "intron_hpRNA", "miRNA", "hairpin", "pri-miRNA", "exon", "TE_DNA", "TE_LINE", "TE_Satellite", "TE_LTR", "TE_RNA", "TE_Other", "TE_Unknown", "TE_RC", "TE_ARTEFACT", "intron_SteXh", "intron_Ste12DOR-RB", "intron_Ste", "intron_sisRNA", "intron", "mRNA", "piRNA_cluster", "intron_snmRNA", "snmRNA", "intron_His-Psi", "pseudogene_His-Psi", "pseudogene", "intron_asRNA", "asRNA", "no_annotation"]

   for suffix in suffixes:
     fin = open(prefix+"/"+prefix+"-"+suffix+".txt", "r")
     read_ids = {}

     #Make a dictionary with read_ids as the key and their count for the value
     for line in fin:
         values = line.split("\t")
         if (values[3] not in read_ids):
             read_ids[values[3]] = 1
         else:
             read_ids[values[3]] = read_ids[values[3]]+1
     fin.close()

     #Check each read_id in the input file and output its normalized read count
     fin = open(prefix+"/"+prefix+"-"+suffix+".txt", "r")
     fout = open(prefix+"/"+prefix+"-"+suffix+"_normalized.txt", "w+")
     for line in fin:
         values = line.split("\t")
         num_reads = int(values[3].split("-")[1])
         normalization = float(read_ids[values[3]])
         normalized = num_reads/normalization

         fout.write(line.rstrip()+"\t"+str(normalized)+"\n")
     fin.close()
     fout.close()

     #Calculate the total number of reads from each class (done two ways, these numbers should be the same)
     total1 = 0
     for read_id in read_ids:
         total1 = total1+int(read_id.split("-")[1])

     total2 = 0
     fin = open(prefix+"/"+prefix+"-"+suffix+"_normalized.txt", "r")
     for line in fin:
         total2 = total2+float(line.split("\t")[12])
     fin.close()

     if (abs(total1-total2) <= 1):
         print(prefix, suffix, ":", total1)
         logfile.write(prefix+" "+suffix+": "+str(total1)+"\n")

     else:
         print("WARNING!", prefix, suffix, "- Total from unique read ids:", total1)
         print("WARNING!", prefix, suffix, "- Total from summed normalized reads:", total2)
         logfile.write("\nWARNING! Unexpected difference in read count calculations:\n")
         logfile.write(prefix+" "+suffix+" from unique: "+str(total1)+"\n")
         logfile.write(prefix+" "+suffix+" from unique: "+str(total2)+"\n\n")

   logfile.close()
