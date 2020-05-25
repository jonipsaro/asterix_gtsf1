# RIP Data Processing Pipeline

## Introduction
RNA Immunoprecipitation allows for characterization of RNAs that co-purify with an isolated protein of interest.
This workflow begins data generated using the SMARTer smRNA-Seq Kit for Illumina sequencing (Takara). After sequencing, the fastq data are trimmed, size-selected, and quality-selected to generate the final dataset.

## Overview
The workflow is divided into the following steps:
1. Prepare, sequence, and demultiplex all sequencing libraries such that each fastq file corresponds to a single sample.
2. Trim the reads to remove the first 3 nucleotides and polyA tails (these are generated during library construction).
3. Size-select
4. Quality-select<BR>
**Most steps of this workflow can be completed with commonly available scripts on platforms such as Galaxy. One script for removal of the polyA-tail can be run in Python3.*

## Step 1 - Prepare, sequence, and demultiplex libraries

## Step 2 - Trim reads
1. The first 3 nucleotides can be removed in many ways. In our workflow, we used trimmomatic to simultaneously perform initial quality filtering<BR>
	```trimmomatic -threads 1 -phred33 fastq_in.fastqsanger fastq_out.fastqsanger HEADCROP:3 SLIDINGWINDOW:4:20```
2. PolyA tail removal was performed using *polyA_trim.py*. 
The script will look for the first instance of an 8 nucleotide polyA and crop the sequence (this includes some tolerance for small deviations if they are near to a large polyA span).
The length of these spans can be modified by replacing the variable $polyA and modifying the get_polyA function if desired.
The user will be prompted for the file name to operate on.<BR>
	```python polyA_trim.py```

## Step 3 - Size selection
This step was performed via Galaxy using the FASTX manipulation tools. In our workflow, we selected for sequences 30-85 nucleotides in length (inclusive).

## Step 4 - Quality-cutoff
This step was also performed via Galaxy using the FASTX manipulation tools.  In our workflow, we imposed a quality cutoff with 95% of bases in the sequence being obligated to pass a threshold score of 24.
