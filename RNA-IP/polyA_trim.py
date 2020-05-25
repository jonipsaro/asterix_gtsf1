#This script is made for processing next generation sequencing data
#in Illumina format.
#
#Sequencing reads are read in, and scanned for polyA stretches.
#Sequences are truncated at long PolyA stretches.
#[The length of these stretches can be adjusted using the variable polyA.]
#
#Also, if there is only a single nucleotide before a long strech of polyA
#that nucleotide will be ignored provided there are 4 or more A's preceeding it.
#
#Trimmed sequences are output to a new file.

#Example sequences and output:
#1 input.  AGTCAGTGTTGCAGTAGCGTCGGACAAAAAAAAAAAAAAAAA
#1 output. AGTCAGTGTTGCAGTAGCGTCGGAC
#2 input.  GTCTCGAATGTGGTTGGGACTGACGTGTTTAAAAAAAAAAAAAAATGTAAAAAAAAA
#2 output. GTCTCGAATGTGGTTGGGACTGACGTGTTT
#3 input.  GTCTCGAATGTGGTTGGGACTGACGTGTTTAAAATAAAAAAAAA
#3 output. GTCTCGAATGTGGTTGGGACTGACGTGTTT


#functions
def get_polyA(seq,polyA):
    first_trim = seq[0:polyA]      #remove any sequence that is downstream of a long stretch of polyA

    if first_trim[len(first_trim)-5:len(first_trim)-1] == "AAAA":   #determine if only a single nucleotide separates the
        second_trim = first_trim[0:len(first_trim)-1]               #...main A stretch from a shorter one upstream
        lastC = second_trim.rfind("C")
        lastG = second_trim.rfind("G")
        lastT = second_trim.rfind("T")

        lastN = max(lastC, lastG, lastT)    #finds the position of the last non-A nucleotide

        third_trim= second_trim[0:lastN+1]

        global second
        second = second + 1

        return third_trim
        
    else:
        global first
        first = first + 1
        return first_trim



#initialize variables
first = 0   #number of reads that are trimmed simply with polyA
second = 0  #number of reads that are trimmed a second time by skipping the last nucleotide

total = 0       #total number of reads read
excluded = 0    #number of reads without polyA

#get file name
seqfile = input("Please enter the original file name: ")
seqfile_name = seqfile[0:seqfile.rfind(".")]
seqfile_ext = seqfile[seqfile.rfind("."):len(seqfile)]

seqfile_Atrim = seqfile_name+"_polyA_trimmed"+seqfile_ext

#open files
seqfile = open(seqfile, "r")
seqfile_Atrim = open(seqfile_Atrim, "w+")

seqid = seqfile.readline()      #read the first sequence id

while seqid!="":     
    seq = seqfile.readline()    #read the remaining 3 lines corresponding to the first read
    plus = seqfile.readline()
    quality = seqfile.readline()

    polyA = seq.find("AAAAAAAA")  #determine if the sequence has 8 nt polyA and find the first instance
    
    if polyA == -1:             #count the total number of reads and those that are exculded (since no polyA)
        excluded = excluded+1   #do not write out the reads that do not have polyA
        total = total+1

    else:
        seqfile_Atrim.write(seqid)

        trimmed = get_polyA(seq, polyA)
        seqfile_Atrim.write(trimmed+"\n")
        seqfile_Atrim.write(plus)
        seqfile_Atrim.write(quality[0:len(trimmed)])

        total = total+1


    seqid = seqfile.readline()  #read the next sequence id. this is done at the end of the loop
                                #to avoid an extra \n at the end of the document
    if seqid!="" and polyA!=-1:
        seqfile_Atrim.write("\n")


#close files
seqfile.close()
seqfile_Atrim.close()

#print total and number of excluded
print("Total number of reads: ", total)
print("Total number excluded (no polyA): ", excluded)
print("Percent excluded: ", excluded/total*100)
print("Simply trimmed by polyA: ", first)
print("Trimmed with a second round of polyA: ", second)
