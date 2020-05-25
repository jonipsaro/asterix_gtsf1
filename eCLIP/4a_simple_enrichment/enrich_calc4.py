#Jonathan Ipsaro
#Cold Spring Harbor Laboratory
#Last reviewed October 18, 2019

#NOTE:  Users must change the yield_* and reads_* variables according to thier libraries.

#This script will determine fold_enrichment between the IP and input samples.
#First, a composite dictionary is built that includes each annotation and the
#corresponding normalized read counts.
#Additionally, a composite dictionary with raw_read_counts is built.
#Then read-per-million thresholds are applied to each annotation.
#For annotations with sufficient background and input reads, the fold_enrichment
#is cacluated as:
#(avg_normalized_ip - normalized_background*bk_scale) / avg_normalized_input

prefixes = ["S1", "S2", "S3", "S4", "S5", "S6"]
logfile = open("enrich_calc.log", "w+")

#These are the global scale factors for the background calculation
#They are determined by comparing the concentration of the background_IP to the IP concentrations to account
#for the background contribution in the IP.
#A value of 1 signifies that no background normalization is needed
#A value of 0 signifies that no background will be subtracted

#To calculate the global scale factors, change the measured concentrations and trimmed reads
#Read counts for trimmed_reads variables should be before removal of PCR duplicates

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


print("Background scale factor for IP1: ", bg_globscale1)
print("Background scale factor for IP2: ", bg_globscale2)
print("Background scale factor for IP1: ", bg_globscale1, file=logfile)
print("Background scale factor for IP2: ", bg_globscale2, file=logfile)

annotations_normalized = {}
annotations_raw = {}


#Go through each file and build composite dictionaries with all annotations
#(one dictionary for normalized values and one for raw values)
print("Part 1/6 - Building composite dictionary will all annotations")
for prefix in prefixes:
    fin = open(prefix+"-annotation_pileup_normalized.txt", "r")
    for line in fin:
        values = line.split("\t")
        if values[0] not in annotations_normalized:
            annotations_normalized[values[0]] = [None, None, None, None, None, None]
            annotations_raw[values[0]] = [None, None, None, None, None, None]           
    fin.close()
print("... Done")


#Go through each file and insert the 1/n-normalized read count into its dictionary entry
print("Part 2/6 - Adding 1/n-normalized read counts to each annotation")
i = 0
for prefix in prefixes:
    fin = open(prefix+"-annotation_pileup_normalized.txt", "r")
    for line in fin:
        values = line.split("\t")
        if values[0] not in annotations_normalized: print("ERROR - MISSING KEY")
        else: annotations_normalized[values[0]][i] = float(values[1])
    fin.close()
    i = i+1
print("... Done")


#Go through each file and insert the raw read count into its dictionary entry
print("Part 3/6 - Adding raw read counts to each annotation")
i = 0
total_read_list = []
for prefix in prefixes:                                     #First, calculate the number of reads in each library
    fin = open(prefix+"-annotation_pileup_normalized.txt", "r")
    total_reads = 0
    for line in fin:
        values = line.split("\t")
        total_reads += float(values[1])
    fin.close()

    mio_reads = total_reads/1000000
    total_read_list.append(total_reads)

    fin = open(prefix+"-annotation_pileup_raw.txt", "r")     #Then, calculate the rpm for each annotation
    for line in fin:
        values = line.split("\t")
        if values[0] not in annotations_raw: print("ERROR - MISSING KEY")
        else: annotations_raw[values[0]][i] = float(values[1])/mio_reads
    fin.close()
    i = i+1
print("... Done")

#Check each annotation's raw read counts to make sure they are above rpm_threshold
print("Part 4/6 - Read count thresholding")
rpm_threshold_input = 0.05      #Change thresholds here
rpm_threshold_ip = 0.5          #Note:  For libraries with ~20 mio reads each, these values correspond
threshold_pass = {}             #to 1 read in the input and 10 reads in the IP
total_excluded = 0

for annotation in annotations_raw:
    raw_rc = annotations_raw[annotation]    #rc stands for read_counts
    threshold_pass[annotation] = raw_rc     #by default, add the annotation to the threshold_pass dictionary
    
    if (raw_rc[0] == None): threshold_pass[annotation][0] = None                #change read counts below thresholds to None
    elif (raw_rc[0] < rpm_threshold_input): threshold_pass[annotation][0] = None

    if (raw_rc[1] == None): threshold_pass[annotation][1] = None
    elif (raw_rc[1] < rpm_threshold_input): threshold_pass[annotation][1] = None

    if (raw_rc[2] == None): threshold_pass[annotation][2] = None
    elif (raw_rc[2] < rpm_threshold_input): threshold_pass[annotation][2] = None

    if (raw_rc[3] == None): threshold_pass[annotation][2] = None
    elif (raw_rc[3] < rpm_threshold_ip): threshold_pass[annotation][3] = None

    if (raw_rc[4] == None): threshold_pass[annotation][2] = None
    elif (raw_rc[4] < rpm_threshold_ip): threshold_pass[annotation][4] = None

    if (raw_rc[5] == None): threshold_pass[annotation][2] = None
    elif (raw_rc[5] < rpm_threshold_ip): threshold_pass[annotation][5] = None


    if ((threshold_pass[annotation][1] == None) and (threshold_pass[annotation][2] == None)):           #no inputs
        del threshold_pass[annotation]
        continue
    if ((threshold_pass[annotation][1] == None) and (threshold_pass[annotation][4] != None)) and ((threshold_pass[annotation][2] == None) or (threshold_pass[annotation][5] != None)):  #no matched input/IP pairs
        del threshold_pass[annotation]
        continue
    if (threshold_pass[annotation][0] == None):         #no background input
        del threshold_pass[annotation]
print("... Done")


#Go through each entry of annotations_normalized
print("Part 5/6 - Subtracting Background")
#If the annotation is present in enough libraries, subtract the background
ip1_sub_bg_total = 0
ip2_sub_bg_total = 0

for annotation in annotations_normalized:
    rc = annotations_normalized[annotation]     #rc stands for read_counts
    bg_input = rc[0]                            #get each entry
    ip1_input = rc[1]
    ip2_input = rc[2]
    bg_ip = rc[3]
    ip1 = rc[4]
    ip2 = rc[5]

    #Subtract background for those annotations that passed the threshold
    if annotation in threshold_pass:
        threshold_ip1_input = threshold_pass[annotation][1]
        threshold_ip2_input = threshold_pass[annotation][2]
        threshold_ip1 = threshold_pass[annotation][4]
        threshold_ip2 = threshold_pass[annotation][5]

        #Scale the background IP to the experimental IPs based on the number of reads in their inputs
        #This allows the bg_ip to be subtracted linearly from each IP
        if bg_ip == None:           #set bg_ip_scaled to 0 if there were no reads
            bg_ip_scaled1 = 0
            bg_ip_scaled2 = 0
        else:
            if ip1_input != None: bg_ip_scaled1 = bg_ip*(ip1_input/bg_input)
            if ip2_input != None: bg_ip_scaled2 = bg_ip*(ip2_input/bg_input)

        #Subtract the background and calculate the size of the background-subtracted library
        if (threshold_ip1_input == None) or (threshold_ip1 == None): ip1_sub_bg = None
        else:
            ip1_sub_bg = ip1 - (bg_ip_scaled1*bg_globscale1)
            ip1_sub_bg_total += ip1_sub_bg
        annotations_normalized[annotation].append(ip1_sub_bg)

        if (threshold_ip2_input == None) or (threshold_ip2 == None): ip2_sub_bg = None
        else:
            ip2_sub_bg = ip2 - (bg_ip_scaled1*bg_globscale2)
            ip2_sub_bg_total += ip2_sub_bg
        annotations_normalized[annotation].append(ip2_sub_bg)

    else:
        annotations_normalized[annotation].append(None)
        annotations_normalized[annotation].append(None)

total_read_list.append(ip1_sub_bg_total)
total_read_list.append(ip2_sub_bg_total)


#Go through each entry of annotations_normalized
print("Part 6/6 - Calculating enrichment")
#If the annotation is present in enough libraries, calculate the enrichment
fout = open("S123456-enrich_per_annotation.txt", "w+")
fout2 = open("S123456-excluded.txt", "w+")

ip1_input_scale = total_read_list[1] / total_read_list[6]       #calculate the scale factor to account library size
ip2_input_scale = total_read_list[2] / total_read_list[7]

print("Scale factors for mapped library sizes (after background subtration):")
print("IP1:",ip1_input_scale)
print("IP2:",ip2_input_scale)
print("Scale factors for mapped library sizes (after background subtration):", file=logfile)
print("IP1:",ip1_input_scale, file=logfile)
print("IP2:",ip2_input_scale, file=logfile)

for annotation in annotations_normalized:
    rc = annotations_normalized[annotation]     #rc stands for read_counts
    ip1_input = rc[1]                           #get each entry
    ip2_input = rc[2]
    ip1_sub_bg = rc[6]
    ip2_sub_bg = rc[7]
    
    num_ip = 0
    enrich = -1

    #Subtract background for those annotations that passed the threshold
    if annotation in threshold_pass:
        threshold_ip1_input = threshold_pass[annotation][1]
        threshold_ip2_input = threshold_pass[annotation][2]
        threshold_ip1 = threshold_pass[annotation][4]
        threshold_ip2 = threshold_pass[annotation][5]

        if (threshold_ip1_input == None) or (threshold_ip1 == None): ip1_enrich = 0
        else:
            ip1_enrich = ip1_sub_bg * ip1_input_scale / ip1_input
            num_ip += 1

        if (threshold_ip2_input == None) or (threshold_ip2 == None): ip2_enrich = 0
        else:
            ip2_enrich = ip2_sub_bg * ip2_input_scale / ip2_input
            num_ip += 1

        #Calculate the average enrichment score
        if num_ip == 0: enrich = 0
        else: enrich = (ip1_enrich+ip2_enrich)/num_ip

    #Write the output
    #Note: enrich < 0 means background is higher than the average ip --> discarded
    #Note: enrich = 0 means that no ip reads above threshold were found; this could mean complete depletion
    #but this is hard to treat mathematically since we are working with 1/n normalized reads (rather than full cpm)
    if (enrich > 0):
        fout.write(annotation+"\t"+str(rc[0])+"\t"+str(rc[1])+"\t"+str(rc[2])+"\t"+str(rc[3])+"\t"+str(rc[4])+"\t"+str(rc[5])+"\t"+str(enrich)+"\n")
    else:
        total_excluded = total_excluded+1
        fout2.write(annotation+"\t"+str(rc[0])+"\t"+str(rc[1])+"\t"+str(rc[2])+"\t"+str(rc[3])+"\t"+str(rc[4])+"\t"+str(rc[5])+"\t"+str(enrich)+"\n")

fout.close()
fout2.close()
print("... Done")


#Print out summary statistics
print("Total number of annotations: ",len(annotations_normalized))
print("Number of annotations excluded (background was higher than IP): ", total_excluded)
logfile.close()
