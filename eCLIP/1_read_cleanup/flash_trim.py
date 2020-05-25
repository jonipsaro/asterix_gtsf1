#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 17, 2019

#Notes:  There are some magic numbers/values in this script
#They include the adapter sequences searched for.

#This script is made for processing next generation sequencing data
#in Illumina format.

#This script was designed to be used on eCLIP data from paired-end libraries
#that were joined with the program FLASH using the following command:

#This script will then trim the FLASH output to remove all calls that are in the adapters
#by comparing the FLASH output to the paired end input.

#Additionally, this script will check for a single duplication of the 5'
#universal Illumina adapter and remove it if it is present.

from datetime import datetime

def rev_comp(forward):
    reverse = forward[::-1]
    rc = reverse.lower().replace("a", "T").replace("t", "A").replace("c", "G").replace("g", "C")
    return rc
 
#initialize variables
total_in = 0    #total number of reads
excluded = 0    #number of reads that could not be processed
trimmed = 0     #number of 'full length' reads with both adapters
double_adapter = 0  #number of reads with a double 5' adapter that is removed

#get file names
read1file = input("Please enter the read 1 file name: ")
read1file_name = read1file[0:read1file.rfind(".")]
read1file_ext = read1file[read1file.rfind("."):len(read1file)]

read2file = input("Please enter the read 2 file name: ")
read2file_name = read2file[0:read2file.rfind(".")]
read2file_ext = read2file[read2file.rfind("."):len(read2file)]

flashfile = input("Please enter the FLASH output file name: ")
flashfile_name = flashfile[0:flashfile.rfind(".")]
flashfile_ext = flashfile[flashfile.rfind("."):len(flashfile)]

flashfile_trimmed = flashfile_name+"_trimmed"+flashfile_ext  #file of FLASH trimmed reads

logfile = flashfile_name+"_trimmed.log"

#open files
read1file = open(read1file, "r")
read2file = open(read2file, "r")
flashfile = open(flashfile, "r")
flashfile_trimmed = open(flashfile_trimmed, "w+")
logfile = open(logfile, "w+")

print(str(datetime.now()))

print("Run started at: ",str(datetime.now()), file=logfile)
print("Read 1 file name: ",read1file, file=logfile)
print("Read 2 file name: ",read2file, file=logfile)
print("FLASH file name: ",flashfile, file=logfile)


ends = {}       #make a dictionary that will have the seq_id as the key and the ends as a list

read1_id = read1file.readline()      #read the first sequence id
read2_id = read2file.readline()


##Build a dictionary that lists each sequence id with the corresponding sequence trimming ends
while read1_id !="":
    seq1 = read1file.readline()     #read in lines from the sequence files
    end1 = seq1[0:20]               #figure out what the end of the forward read is
    plus1 = read1file.readline()
    quality1 = read1file.readline()

    seq2 = read2file.readline()
    end2 = seq2[0:20]
    end2_rc = rev_comp(end2)        #figure out the reverse complement of the reverse read
    plus2 = read2file.readline()
    quality2 = read2file.readline()

    ends[read1_id] = [end1, end2_rc]    #add the ends to the dictionary

    read1_id = read1file.readline()      #read the next sequence id
    read2_id = read2file.readline()



#Go through the FLASH output file and use the timming ends above to trim the output
flash_id = flashfile.readline()

while flash_id !="":     
    seq = flashfile.readline()    #read the remaining 3 lines corresponding to the first read
    plus = flashfile.readline()
    quality = flashfile.readline()

    trim1 = seq.find(ends[flash_id][0])     #trim1 and trim2 indicate the "ends" of the read as
    trim2 = seq.find(ends[flash_id][1])     #determined by the beginning of read1 and beginning of read2

    trim3 = seq.find("ACACGACGCTCTTCCGATCT")    #trim3 and trim4 indicate the "ends" of the read as determined
    trim4 = seq.find("AGATCGGAAGAGCACACGTC")    #by where the adapters stop.  This recovers some reads that
                                                #have been joined and have an incorrect call

    if (trim1 != -1) and (trim2 != -1):
        flashfile_trimmed.write(flash_id)
        seq_check = seq[trim1:trim2+20]
        qual_check = quality[trim1:trim2+20]
        
        if (seq_check.find("TCTTCCGATCT") == -1):       #checks to see if a duplicated adapter is present
            flashfile_trimmed.write(seq_check+"\n")     #if not, write the trimmed output
            flashfile_trimmed.write(plus)
            flashfile_trimmed.write(qual_check+"\n")

        else:
            flashfile_trimmed.write(seq_check[seq_check.find("TCTTCCGATCT")+11:]+"\n")  #if a duplicated adapter is present, remove it
            flashfile_trimmed.write(plus)                                               #before writing the output
            flashfile_trimmed.write(qual_check[seq_check.find("TCTTCCGATCT")+11:]+"\n")
            double_adapter = double_adapter+1
            
        trimmed = trimmed+1
        total_in = total_in+1
        
    elif (trim1 == -1) and (trim2 == -1) and (trim3 != -1) and (trim4 !=-1):
        flashfile_trimmed.write(flash_id)
        seq_check = seq[trim3+20:trim4]
        qual_check = quality[trim3+20:trim4]
        
        if (seq_check.find("TCTTCCGATCT") == -1):
            flashfile_trimmed.write(seq_check+"\n")
            flashfile_trimmed.write(plus)
            flashfile_trimmed.write(qual_check+"\n")

        else:
            flashfile_trimmed.write(seq_check[seq_check.find("TCTTCCGATCT")+11:]+"\n")
            flashfile_trimmed.write(plus)
            flashfile_trimmed.write(qual_check[seq_check.find("TCTTCCGATCT")+11:]+"\n")
            double_adapter = double_adapter+1

        trimmed = trimmed+1
        total_in = total_in+1

    elif (trim1 == -1) and (trim2 != -1) and (trim3 != -1):
        flashfile_trimmed.write(flash_id)
        seq_check = seq[trim3+20:trim2+20]
        qual_check = quality[trim3+20:trim2+20]
        
        if (seq_check.find("TCTTCCGATCT") == -1):
            flashfile_trimmed.write(seq_check+"\n")
            flashfile_trimmed.write(plus)
            flashfile_trimmed.write(qual_check+"\n")

        else:
            flashfile_trimmed.write(seq_check[seq_check.find("TCTTCCGATCT")+11:]+"\n")
            flashfile_trimmed.write(plus)
            flashfile_trimmed.write(qual_check[seq_check.find("TCTTCCGATCT")+11:]+"\n")
            double_adapter = double_adapter+1

        trimmed = trimmed+1
        total_in = total_in+1

    elif (trim2 == -1) and (trim1 != -1) and (trim4 != -1):
        flashfile_trimmed.write(flash_id)
        seq_check = seq[trim1:trim4]
        qual_check = quality[trim1:trim4]

        if (seq_check.find("TCTTCCGATCT") == -1):
            flashfile_trimmed.write(seq_check+"\n")
            flashfile_trimmed.write(plus)
            flashfile_trimmed.write(qual_check+"\n")

        else:
            flashfile_trimmed.write(seq_check[seq_check.find("TCTTCCGATCT")+11:]+"\n")
            flashfile_trimmed.write(plus)
            flashfile_trimmed.write(qual_check[seq_check.find("TCTTCCGATCT")+11:]+"\n")
            double_adapter = double_adapter+1

        trimmed = trimmed+1
        total_in = total_in+1

    else:
        excluded = excluded+1
        total_in = total_in+1

    flash_id = flashfile.readline()  #read the next sequence id. this is done at the end of the loop
                                #to avoid an extra \n at the end of the document

#print run summary
print("Total number of reads input: ", total_in)
print("Total number excluded: ", excluded)
print("Total number trimmed: ", trimmed)
print("Numer of reads with double 5' adapter: ", double_adapter)
print("Percent excluded: ", excluded/total_in*100)
print("Percent trimmed:  ", trimmed/total_in*100)
print(str(datetime.now()))

print("Total number of reads input: ", total_in, file=logfile)
print("Total number excluded: ", excluded, file=logfile)
print("Total number trimmed: ", trimmed, file=logfile)
print("Numer of reads with double 5' adapter: ", double_adapter, file=logfile)
print("Percent excluded: ", excluded/total_in*100, file=logfile)
print("Percent trimmed:  ", trimmed/total_in*100, file=logfile)
print("Run ended at: ", datetime.now(), file=logfile)

#close files
read1file.close()
read2file.close()
flashfile.close()
flashfile_trimmed.close()
logfile.close()
