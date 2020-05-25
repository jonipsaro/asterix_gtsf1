#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 17, 2019

#This script will take a bed file and filter it sequentially through individual annotation classes.
#Once a read is assigned to a class, it is blacklisted from lower priorty classes.

#Important:  The ordered_features list must reflect the classes present in the annotated file.

import sys
arguments = sys.argv
if len(sys.argv) == 1:
   input("By default, this script will run samples [\"S1\", \"S2\", \"S3\", \"S4\", \"S5\", \"S6\"].\nTo run a single file, include it in the argument.\n\nPress Enter to continue or Ctrl+C to break")
   samples = ["S1", "S2", "S3", "S4", "S5", "S6"]
else:
   samples = sys.argv[1:]


for sample_name in samples:
   folder = "abundance"

   if folder == "abundance":        #For filtering based on expected abundance, use this list:
      ordered_features = ["rRNA", "intron_28SrRNA-Psi", "pseudogene_5SrRNA-Psi", "pseudogene_5.8SrRNA-Psi", "pseudogene_28SrRNA-Psi", "pseudogene_2SrRNA-Psi", "pseudogene_18SrRNA-Psi", "snoRNA", "snRNA", "structRNA_7SLRNA", "structRNA_RNaseP", "structRNA_RNaseMRP", "structRNA_sbRNA", "structRNA_scaRNA", "intron_lncRNA", "lncRNA", "intron_tRNA", "tRNA", "pseudogene_tRNA", "mttRNA", "intron_hpRNA", "miRNA", "hairpin", "pri-miRNA", "exon", "TE_DNA", "TE_LINE", "TE_Satellite", "TE_LTR", "TE_RNA", "TE_Other", "TE_Unknown", "TE_RC", "TE_ARTEFACT", "intron_SteXh", "intron_Ste12DOR-RB", "intron_Ste", "intron_sisRNA", "intron", "mRNA", "piRNA_cluster", "intron_snmRNA", "snmRNA", "intron_His-Psi", "pseudogene_His-Psi", "pseudogene", "intron_asRNA", "asRNA", "no_annotation"]
   elif folder == "stringent":      #For higher stringency filtering that de-prioritizes tRNas, use this list:
      ordered_features = ["rRNA", "intron_28SrRNA-Psi", "pseudogene_5SrRNA-Psi", "pseudogene_5.8SrRNA-Psi", "pseudogene_28SrRNA-Psi", "pseudogene_2SrRNA-Psi", "pseudogene_18SrRNA-Psi", "snoRNA", "snRNA", "structRNA_7SLRNA", "structRNA_RNaseP", "structRNA_RNaseMRP", "structRNA_sbRNA", "structRNA_scaRNA", "intron_lncRNA", "lncRNA", "intron_hpRNA", "miRNA", "hairpin", "pri-miRNA", "exon", "intron_SteXh", "intron_Ste12DOR-RB", "intron_Ste", "intron_sisRNA", "intron", "mRNA", "intron_snmRNA", "snmRNA", "intron_His-Psi", "pseudogene_His-Psi", "pseudogene", "intron_asRNA", "asRNA", "piRNA_cluster", "TE_DNA", "TE_RC", "TE_LINE", "TE_Satellite", "TE_LTR", "TE_RNA", "TE_Other", "TE_Unknown", "TE_ARTEFACT", "intron_tRNA", "tRNA", "pseudogene_tRNA", "mttRNA", "no_annotation"]
   else: print("There is a problem with your folder designation. Please check this script (around line 16).")

   blacklist = set()
   logfile = open("../10_filtered/"+folder+"/"+sample_name+"/"+sample_name+".log", "w+")
   logfile.write("Filtering will be performed in the following order:\n")
   logfile.write(str(ordered_features)+"\n\n")
   logfile.write("The following values indicate the number of elements that were ignored because the read was already sorted to a higher priority class and blacklisted:\n")

   for feature in ordered_features:
      fin = open(sample_name+"_all.txt", "r")
      fout = open("../10_filtered/"+folder+"/"+sample_name+"/"+sample_name+"-"+feature+".txt", "w+")
      blacklist_add = set()
      blacklist_count = 0

      for line in fin:
         values = line.split("\t")
         if ((values[9].split(":")[0] == feature) and (values[3] not in blacklist)):
            fout.write(line)
            blacklist_add.add(values[3])

         elif (values[9].split(":")[0] == feature):
            blacklist_count = blacklist_count+1

      fin.close()
      fout.close()
      blacklist = blacklist|blacklist_add
      logfile.write("Number of "+feature+"s ignored: "+str(blacklist_count)+"\n")


   missing_features = 0
   fin = open(sample_name+"_all.txt", "r")
   for line in fin:
      values = line.split("\t")
      if (values[9].split(":")[0] not in ordered_features):
         logfile.write("WARNING - Missing feature: "+values[9].split(":")[0]+"\n")
         missing_features = missing_features+1
         ordered_features.append(values[9].split(":")[0])

   print("Number of missing features:",missing_features)
   fin.close()
   logfile.close()
