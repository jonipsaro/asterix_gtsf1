#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 17, 2019

#This script is made for processing next generation sequencing data
#in Illumina format.

#This script will generate the reverse complement of each read.
#Quality scores will also be reversed.
#Sequence identifiers will not be changed.

from datetime import datetime

def rev_comp(forward):
    reverse = forward[::-1]
    rc = reverse.lower().replace("a", "T").replace("t", "A").replace("c", "G").replace("g", "C")
    return rc

#initialize variables
total_in = 0    #total number of reads

#get file names
seqfile = input("Please enter the file name: ")
seqfile_name = seqfile[0:seqfile.rfind(".")]
seqfile_ext = seqfile[seqfile.rfind("."):len(seqfile)]

seqfile_rc= seqfile_name+"_rc"+seqfile_ext


#open files
seqfile = open(seqfile, "r")
seqfile_rc = open(seqfile_rc, "w+")

print(str(datetime.now()))


seq_id = seqfile.readline()

while seq_id !="":     
    seq = seqfile.readline()    #read the remaining 3 lines corresponding to the first read
    plus = seqfile.readline()
    quality = seqfile.readline()

    seqfile_rc.write(seq_id)
    seqfile_rc.write(rev_comp(seq.strip())+"\n")
    seqfile_rc.write(plus)
    seqfile_rc.write(quality.strip()[::-1]+"\n")

    total_in = total_in+1

    seq_id = seqfile.readline()  #read the next sequence id. this is done at the end of the loop
                                 #to avoid an extra \n at the end of the document

#print run summary
print("Total number of reads processed: ", total_in)
print(str(datetime.now()))

#close files
seqfile.close()
seqfile_rc.close()
