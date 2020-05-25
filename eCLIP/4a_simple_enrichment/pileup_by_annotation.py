#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 18, 2019

#This script will take a tab-delimited text file with the following format (all in one row):
#chr start stop read_id-count score strand[+/-]  <-- [Columns 0-5] Information for the mapped read
#chr start stop annotation score strand[+/-] <-- [Columns 6-11] Information for the corresponding annotation
#normalized_count <--- [Column 12]

#Importantly, the read_id-count column indicates the id of the collapsed read and the number of times that read was observed.

#Two files are generated.  Both files 1 and 2 are used for straightforward processing.  File 2 can be used to prepare inputs for DESeq2
#1.  -annotation_pileup_normalized.txt will have the following tab-delimited columns
#annotation     norm_total_pileup (based on in-class, multimap-normalized read counts)   rpm_norm_total_pileup (as before, but rpm)

#2.  -annotation_pileup_raw.txt will have the following tab-delimited columns
#annotation     raw_total_pileup (does NOT use multimap normalization)

from datetime import datetime
import sys
arguments = sys.argv
if len(sys.argv) == 1:
   input("By default, this script will run samples [\"S1\", \"S2\", \"S3\", \"S4\", \"S5\", \"S6\"].\nTo run a single file, include it in the argument.\n\nPress Enter to continue or Ctrl+C to break")
   prefixes = ["S1", "S2", "S3", "S4", "S5", "S6"]
else:
   prefixes = sys.argv[1:]
   
for prefix in prefixes:
    fin = open(prefix+"/"+prefix+"-all_normalized.txt", "r")
    fout = open(prefix+"/"+prefix+"-annotation_pileup_normalized.txt", "w+")
    fout2 = open(prefix+"/"+prefix+"-annotation_pileup_raw.txt", "w+")
    logfile = open(prefix+"/"+prefix+"-annotation_pileup.log", "w+")
    print("Run started at: ", datetime.now(), file=logfile)

    annotations_normalized = {}
    annotations_raw = {}
    total_norm_reads = 0
    total_raw_reads = 0

    for line in fin:
        values = line.split("\t")
        if values[9] not in annotations_normalized:
            annotations_normalized[values[9]] = float(values[12])
        else:
            annotations_normalized[values[9]] = annotations_normalized[values[9]]+float(values[12])
        total_norm_reads = total_norm_reads+float(values[12])

        if values[9] not in annotations_raw:
            annotations_raw[values[9]] = int(values[3].split("-")[1])
        else:
            annotations_raw[values[9]] = annotations_raw[values[9]]+int(values[3].split("-")[1])
        total_raw_reads = total_raw_reads + int(values[3].split("-")[1])
        
    
    for annotation in annotations_normalized:
        fout.write(annotation+"\t"+str(annotations_normalized[annotation])+"\t"+str((annotations_normalized[annotation]/total_norm_reads)*1000000)+"\n")

    for annotation in annotations_raw:
        fout2.write(annotation+"\t"+str(annotations_raw[annotation])+"\n")

    print("Number of annotations tabulated in "+prefix+":",len(annotations_normalized))
    print("Number of normalized reads in "+prefix+":",int(round(total_norm_reads, 0)))
    print("Number of raw reads (multi-counts multi-mappers) in "+prefix+":",total_raw_reads)

    print("Number of annotations tabulated in "+prefix+":",len(annotations_normalized), file=logfile)
    print("Number of normalized reads in "+prefix+":",int(round(total_norm_reads, 0)), file=logfile)
    print("Number of raw reads (multi-counts multi-mappers) in "+prefix+":",total_raw_reads, file=logfile)
    print("Run ended at: ", datetime.now(), file=logfile)

    fin.close()
    fout.close()
    fout2.close()
    logfile.close()
