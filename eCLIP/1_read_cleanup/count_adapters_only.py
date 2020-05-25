#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 17, 2019

#Notes:  There are some magic numbers/values in this script
#They include the adapter_dimer_length list (corresponding sizes for each library)
#and the adapter sequences searched for.

#This script is made for processing next generation sequencing data
#in Illumina format.
#
#Sequencing reads are read in, and scanned for size to determine how many reads
#are from adapter dimers.

from datetime import datetime

#initialize variables
#input the prefix for each file to be processed
#the remainder of the file name can be adjusted below in seqfile_name
files = ["S1", "S2", "S3", "S4", "S5", "S6"]

#indicate the corresponding size cutoff to be called an adapter dimer
#for eCLIP experiments, there are at least 10 nt for the inputs and 22 for the IPs
#I typically allow an additional 2 nt
adapter_dimer_length = [12, 12, 12, 24, 24, 24]

i = 0


#get file name
for each in files:
    total_in = 0            #total number of reads
    adapter_dimers = 0      #number of adapter_dimers
    true_adapter = 0        #number of verified adapters (only works for IP reads)
    
    seqfile_name = each+".extendedFrags_trimmed.fastq"

    logfile = "adapter_dimers.log"


#open files
    seqfile = open(seqfile_name, "r")
    logfile = open(logfile, "a+")

    print(str(datetime.now()))

    seqid = seqfile.readline()      #read the first sequence id

#check each read for its length
    while seqid != "":     
        seq = seqfile.readline().rstrip()       #read the remaining 3 lines corresponding to the first read
        plus = seqfile.readline()
        quality = seqfile.readline()

        if (len(seq) <= adapter_dimer_length[i]):                    #if the sequence is within the length of the adapters + 2nt
            adapter_dimers = adapter_dimers+1
            if (seq.find("CCTATAT") != -1) or (seq.find("TGCTATT") != -1): true_adapter += 1

        total_in = total_in+1
        seqid = seqfile.readline()
            
    print("total number of reads: ", total_in)
    print("number of reads that are adapter dimers", adapter_dimers)
    print("number of reads that are verified adapters", true_adapter)

    logfile.write(seqfile_name+"\ntotal reads: "+str(total_in)+"\nadapterdimers: "+str(adapter_dimers)+"\n")

    i += 1

#close files
    seqfile.close()
    logfile.close()
