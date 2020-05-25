#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed November 08, 2019

#This script will take output from DESeq2 for three combinations of comparisons:
#1. ip versus input
#2. ip versus background
#3. background versus input

#The user must specifiy the yield and reads variables below

#The input files should have the same number of lines and the same annotations (column 1)
#But the annotations do not need to be identically ordered.

#The output file will have the following tab-delimited columns:
#gene (or logical gene) name
#fold change (ip-(background*bg_scale))/input
#p-value[adjusted] (less significant value of [ip versus input; ip versus background])

#Note: DESeq2 output includes log2-fold changes.  These are converted to fold changes (no log2).

ip_v_input = open("S12356-deseq2_output.tabular", "r")      #Check that each filename matches the corresponding DESeq2 output
ip_v_bg = open("S456-deseq2_output.tabular", "r")
bg_v_input = open("S1234-deseq2_output.tabular", "r")

#Set variables for global scale correction 
#yield (concentration * volume)
yield_bg = 
yield_ip1 = 
yield_ip2 = 

#adapter trimmed reads - adapter dimers
reads_bg = 
reads_ip1 = 
reads_ip2 = 

bg_globscale1 = (yield_bg/yield_ip1)*(reads_bg/reads_ip1)
bg_globscale2 = (yield_bg/yield_ip2)*(reads_bg/reads_ip2)
bg_scale = (bg_globscale1+bg_globscale2)/2


#Create a dictionary with annotations as keys and
#[foldchange_ip_v_input, foldchange_ip_v_bg, foldchange_bg_v_input, p-value-adj_ip_v_input, p-value-adj_ip_v_bg]
#Again, note that the DESeq2 output is log2(fold_change).  Those are converted to fold_change.
deseq2 = {}
for line in ip_v_input.readlines()[1:]:
    values = line.split("\t")
    annotation = values[0].replace("\"", "")
    if annotation not in deseq2: deseq2[annotation] = [None,None,None,None,None]
    
    if (values[2] == "NA"): deseq2[annotation][0] = values[2]
    else: deseq2[annotation][0] = 2**float(values[2])
    
    if (values[6].strip() == "NA"): deseq2[annotation][3] = values[6].strip()
    else: deseq2[annotation][3] = float(values[6].strip())
ip_v_input.close()

for line in ip_v_bg.readlines()[1:]:
    values = line.split("\t")
    annotation = values[0].replace("\"", "")
    if annotation not in deseq2: print("ERROR - Missing Key: Check input files to ensure the first column is consistent between files")

    if (values[2] == "NA"): deseq2[annotation][1] = values[2]
    else: deseq2[annotation][1] = 2**float(values[2])

    if (values[6].strip() == "NA"): deseq2[annotation][4] = values[6].strip()
    else: deseq2[annotation][4] = float(values[6].strip())
ip_v_bg.close()

for line in bg_v_input.readlines()[1:]:
    values = line.split("\t")
    annotation = values[0].replace("\"", "")
    if annotation not in deseq2: print("ERROR - Missing Key: Check input files to ensure the first column is consistent between files")

    if (values[2] == "NA"): deseq2[annotation][2] = values[2]
    else: deseq2[annotation][2] = 2**float(values[2])
bg_v_input.close()


#Go through each entry of deseq2, calculate background-subtracted fold-change and output
fout = open("S123456-deseq2_polished.txt", "w+")
NA_count = 0

for annotation in deseq2:
    if (None in deseq2[annotation]):    #Check to make sure all columns have been filled
        print("ERROR - Missing Value: Check input files for consistent annotation column")
    elif ("NA" in deseq2[annotation]):  #Check DESeq2 output for "NA" values.  Exclude those entries.     
        NA_count += 1
    else:
        fout.write(annotation+"\t")                                                                 #annotation
        fout.write(str(float(deseq2[annotation][0])-(float(deseq2[annotation][2])*bg_scale))+"\t")  #background-subtracted fold enrichment
        pvalue = max(float(deseq2[annotation][3]), float(deseq2[annotation][4]))                    #p-value-adj (worse of ip_v_bg, ip_v_input)
        fout.write(str(pvalue)+"\n")

fout.close()

#Print out summary statistics
print("Total number of annotations: ",len(deseq2))
print("Total number of NA entries: ",NA_count)
