#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 17, 2019

#!!!Note: For eCLIP, input and CLIP samples have slightly different barcoding
#strategies.  Make sure you are running this script only for CLIP samples!!!

#This script is made for processing next generation sequencing data
#in Illumina format.

#Sequencing reads are read in, and scanned for user defined barcodes.
#The PCR duplicates are removed and the barcodes are trimmed off.

#This script was designed to be used on eCLIP data that has had linker adapters
#removed.

#This script will produce four output files.
#The main output file will contain all unique reads that have been cropped
#to remove barcodes.
#The second output file will contain a list of PCR duplicates (uncropped).
#The third output file will contain a list of any CLIP reads that did not have
#a proper tag (at positions 6-12).
#The fourth output file will contain a list of reads that have problems (either
#unassigned bases or that are too short (<38 nt including the 22 nt of bardcodes).

from datetime import datetime

#initialize variables
total_in = 0            #total number of reads
unique = 0              #number of unique umis
duplicates = 0          #number of reads with double 3' adapter
bad_tag = 0             #number of CLIP reads that do not have a correct tag (positions 6-12)
bad_read = 0            #number of reads with undetermined nucleotide or of too short length


#the goal of this dictionary is to speedup lookup by dividing the hash table into 64 smaller ones
umi_dict = {
    "AAA": set(),"AAC": set(),"AAG": set(),"AAT": set(),"ACA": set(),"ACC": set(),"ACG": set(),"ACT": set(),
    "AGA": set(),"AGC": set(),"AGG": set(),"AGT": set(),"ATA": set(),"ATC": set(),"ATG": set(),"ATT": set(),
    "CAA": set(),"CAC": set(),"CAG": set(),"CAT": set(),"CCA": set(),"CCC": set(),"CCG": set(),"CCT": set(),
    "CGA": set(),"CGC": set(),"CGG": set(),"CGT": set(),"CTA": set(),"CTC": set(),"CTG": set(),"CTT": set(),
    "GAA": set(),"GAC": set(),"GAG": set(),"GAT": set(),"GCA": set(),"GCC": set(),"GCG": set(),"GCT": set(),
    "GGA": set(),"GGC": set(),"GGG": set(),"GGT": set(),"GTA": set(),"GTC": set(),"GTG": set(),"GTT": set(),
    "TAA": set(),"TAC": set(),"TAG": set(),"TAT": set(),"TCA": set(),"TCC": set(),"TCG": set(),"TCT": set(),
    "TGA": set(),"TGC": set(),"TGG": set(),"TGT": set(),"TTA": set(),"TTC": set(),"TTG": set(),"TTT": set()
    }


#get file name
seqfile = input("Please enter the original file name: ")
seqfile_name = seqfile[0:seqfile.rfind(".")]
seqfile_ext = seqfile[seqfile.rfind("."):len(seqfile)]

seqfile_unique = seqfile_name+"_unique"+seqfile_ext             #file of unique reads
seqfile_duplicates = seqfile_name+"_duplicates"+seqfile_ext     #file of pcr duplicates
seqfile_bad_tag = seqfile_name+"_bad_tags"+seqfile_ext          #file of experimental reads lacking a proper tag
seqfile_bad_read = seqfile_name+"_bad_reads"+seqfile_ext        #file of experimental reads with unassigned bases

logfile = seqfile_name+"_umi_collapse.log"


#open files
seqfile = open(seqfile, "r")
seqfile_unique = open(seqfile_unique, "w+")
seqfile_duplicates = open(seqfile_duplicates, "w+")
seqfile_bad_tag = open(seqfile_bad_tag, "w+")
seqfile_bad_read = open(seqfile_bad_read, "w+")
logfile = open(logfile, "w+")

print(str(datetime.now()))

seqid = seqfile.readline()      #read the first sequence id

while seqid!="":     
    seq = seqfile.readline().rstrip()    #read the remaining 3 lines corresponding to the first read
    plus = seqfile.readline()
    quality = seqfile.readline()

    umi = seq[0:12]+"_"+seq[len(seq)-10:len(seq)]   #makes a composite UMI based on standard eCLIP barcodes
                                                    #CLIP reads should have five random nt followed by a
                                                    #designed 7-mer.
    umi_key = umi[0:3]                              #makes a 3 nucleotide key for the UMI dictionary    

    if ("N" in seq) or (len(seq) < 38):             #if the sequence contains an undeterminened
        total_in = total_in+1                       #nucleotide, or is shorter than 38 nt (barcodes
        bad_read = bad_read+1                           #alone are 22), put it in the bad_reads file

        seqfile_bad_read.write(seqid)
        seqfile_bad_read.write(seq+"\n")
        seqfile_bad_read.write(plus)
        seqfile_bad_read.write(quality)        

    elif umi in umi_dict[umi_key]:                  #if the UMI is already known, write all the read info
        total_in = total_in+1                       #to the file of duplicates
        duplicates = duplicates+1
        
        seqfile_duplicates.write(seqid)
        seqfile_duplicates.write(seq+"\n")
        seqfile_duplicates.write(plus)
        seqfile_duplicates.write(quality)

    elif (umi[5:12] == "CCTATAT") or (umi[5:12] == "TGCTATT"):  #if the UMI is new, check to make sure the
        total_in = total_in+1                                   #linker is correct.  if so, output the
        unique = unique+1                                       #cropped sequence.
        
        seqfile_unique.write(seqid)
        seqfile_unique.write(seq[12:len(seq)-10]+"\n")
        seqfile_unique.write(plus)
        seqfile_unique.write(quality[12:len(seq)-10]+"\n")

        umi_dict[umi_key].add(umi)                              #add the new UMI to the set in the UMI dictionary

    else:                                                       #if the UMI is new but the linker is not
        total_in = total_in+1                                   #correct, output the read info to the
        bad_tag = bad_tag+1                                     #bad_tag file.

        seqfile_bad_tag.write(seqid)
        seqfile_bad_tag.write(seq+"\n")
        seqfile_bad_tag.write(plus)
        seqfile_bad_tag.write(quality)

        #do not add the new UMI to the set. since it has a bad tag, it will be excluded.
        

    seqid = seqfile.readline()  #read the next sequence id. this is done at the end of the loop
                                #to avoid an extra \n at the end of the document

#print run summary
print("Total number of reads input: ", total_in)
print("Total number of duplicates: ", duplicates)
print("Total number with bad tags: ", bad_tag)
print("Total number of bad reads: ", bad_read)
print("Total number of cropped reads kept: ", unique)
print("Percent duplicates: ", duplicates/total_in*100)
print("Percent bad tags: ", bad_tag/total_in*100)
print("Percent bad reads: ", bad_read/total_in*100)
print("Percent retained: ", unique/total_in*100)
print(str(datetime.now()))

print("Total number of reads input: ", total_in, file=logfile)
print("Total number of duplicates: ", duplicates, file=logfile)
print("Total number with bad tags: ", bad_tag, file=logfile)
print("Total number of bad reads: ", bad_read, file=logfile)
print("Total number of cropped reads kept: ", unique, file=logfile)
print("Percent duplicates: ", duplicates/total_in*100, file=logfile)
print("Percent bad tags: ", bad_tag/total_in*100, file=logfile)
print("Percent bad reads: ", bad_read/total_in*100, file=logfile)
print("Percent retained: ", unique/total_in*100, file=logfile)
print("Run ended at: ", datetime.now(), file=logfile)

#close files
seqfile.close()
seqfile_unique.close()
seqfile_duplicates.close()
seqfile_bad_tag.close()
seqfile_bad_read.close()
logfile.close()
