# eCLIP Pipeline - Parameters list and summary

## Read Cleanup - 1_read_cleanup directory
### Joining of paired-end reads with FLASH
```
    inputs		fastq files for the paired-end reads (2)
    output		fastq file with read pairs merged into a single sequence
    parameters
        -m		minimum overlap = 25 nt
        -M		maximum overlap = 175 nt
        -O		allow "outie" orientation
        -d		output directory
        -o		output prefix
    example
	    FLASH-1.2.11/flash -m 25 -M 175 -O -d 2_flash/S1 -o S1 1_raw_data/S1_R1.fastq 1_raw_data/S1_R2.fastq 2>&1 | tee 2_flash/S1/S1.log
```

### Adapter trimming – flash_trim.py
This script trims the FLASH output to remove all calls that are in the adapters by comparing the FLASH output to the paired end input. Additionally, this script will check for a single duplication of the 5' universal Illumina adapter and remove it.
```
    inputs		fastq files for the paired-end reads
                fastq file from FLASH
    output		single fastq file with read pairs merged into a single sequence
    parameters	(specified within flash_trim.py)
    adapters    ACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTC, TCTTCCGATCT
```

### PCR duplicate removal – umi_collapse.py & umi_collapse_inputs.py
This script will remove all reads that have arisen due to PCR duplication by comparing the random barcodes incorporated in each read. Additionally, reads that are very short (&lt; 16 nucleotides after removing barcodes) are excluded. Note that due to differences in the structure of the input and IP barcodes, there are two versions of this script (umi_collapse.py & and umi_collapse_inputs.py).
```
    input	    fastq file [adapter-trimmed, FLASH-merged reads]
    outputs	    fastq file of reads that are not from PCR duplication (*_unique) (main output)
			    fastq file of PCR duplicate reads (*_duplicates)
			    fastq file of reads that did not have a proper tag (*_bad_tags)
			    fastq file of reads that are too short or have unassigned bases (*_bad_reads)
			    logfile with summary counts
```

### Reverse complement – rev_comp.py
This is a very simple script to take the reverse complement of each fastq file. Scores are also reversed. There are equivalent scripts available on other platforms, including Galaxy.
```
    input		fastq file [PCR duplicate-removed, adapter-trimmed, FLASH-merged reads]
    output		fastq file of reverse-complemented reads and reversed quality scores
```

### Quality control check – FastQC

### Collapse – Part of FASTX toolkit, available on Galaxy
Collapse merges identical sequences in a fastq/fasta file into a single entry while maintaining read counts. This will allow for faster performance in later steps of the pipeline.
```
    input		fastq or fasta file [reverse complemented, trimmed, FLASH-merged reads]
    output		fasta file with identical reads collapsed
    example output
        >882-1176
        GCGCGACCCGCTCCGGGGACAGTGCCAGGTGGGGAGTTTGA
```

## Mapping - STAR, bedtools, and 2_mapping_cleanup directory

### Mapping with STAR
Collapsed reads are mapped to the genome using STAR. Importantly, multi-mapping reads were preserved.
```
    inputs		fastq file [read collapsed, reverse complemented, trimmed, FLASH-merged reads]
		        genome files
    outputs	    STAR generates several output files. For more information, see the STAR manual.
                The main alignment file is Aligned.sortedByCoord.out.bam
    parameters
        --genomeDir		                   directory in which the genome files are stored
        --runThreadN	                    number of threads to run
        --readFilesIn	                    input file
        --outSAMtype				        output file type = BAM (option SortedByCoordinate)
        --quantMode				            count number of reads while mapping
						                    (option GeneCounts for reads per gene;
					            	        option TranscriptomeSAM for transcriptome alignment)
        --outSAMAttributes			        output alignment information (see STAR manual)
        --outFilterMultimapNmax		        maximum number of multi-mappers allowed = 50
        --outFilterMismatchNoverReadLmax	maximum ratio of mismatches to read length = 0.05
        --outFilterMismatchNmax		        maximum number of mismatches = 999 = off
    example
        STAR_2.5.2b
        --genomeDir mm10
        --runThreadN 12  
        --readFilesIn  ../S1.extendedFrags_trimmed_unique_rc_collapsed.fasta
        --outSAMtype BAM SortedByCoordinate
        --quantMode GeneCounts TranscriptomeSAM
        --outSAMattributes NH HI nM AS NM MD
        --outFilterMultimapNmax 50
        --outFilterMismatchNoverReadLmax 0.05 
        --outFilterMismatchNmax 999
```

### Convert alignment file to .bed format using bedtools
```
    input		bam file output from STAR
    output		bed file
    example
        bedtools bamtobed -i S1/Aligned.sortedByCoord.out.bam > S1.bed
```

### Remove "multi-aligning" multi-mappers – unique_bed_tidy.py
Output from STAR can produce nearly identical alignments for a given read due to similar alignment scores for indels. These "multi-aligners" are collapsed into a single alignment to avoid normalization artefacts downstream. First, alignments with identical bounds are collapsed to a single alignment using the command line. Then unique_bed_tidy.py collapses overlapping alignments.
```
    input		bed file [generated by bamtobed from STAR .bam file]
    outputs	    bed file of alignments with identical bounds (*_unique.bed)
		        bed file of unique alignments (overlapping aligners collapsed; *_unique_tidy.bed)
    example
        sort S1.bed | uniq > S1_unique.bed
        then run unique_bed_tidy.py using S1_unique.bed as the input file.
```

## Annotation and Filtering - bedtools and 3_annotation_and_filtering directory

### Annotate reads using bedtools
```
    inputs		alignment bed file [converted STAR output, filtered as described above]
		        annotation bed file
    outputs	    plaintext file of mapped locations and the corresponding annotations (*_annotations.txt)
	            plaintext file of mapped locations that do not have annotations (*_no_annotation.txt)
    parameters
    -f		    fraction of overlap required between read and annotation = 1
    -wa	        output the complete read alignment .bed data
    -wb	        output the complete annotation .bed data
    -s		    enforce strandedness
    -a		    file a (alignment file)
    -b		    file b (annotation file)
    example
        bedtools intersect -f 1 -wa -wb -s -a S1_unique_tidy.bed -b mm10_annotations.bed > S1_annotated.txt; \
        bedtools intersect -f 1 -v -wa -s -a 8_bamtobed/S1_unique_tidy.bed -b mm10_annotations.bed > 9_annotated/S1_no_annotation.txt
```

### Format file of reads with no annotation so that it is consistent – format_no_annotation.py
```
    inputs		plaintext file of mapped locations that do not have annotations
    outputs	    reformatted file (*_no_annotation_formatted.txt)
```

### Create a composite file with all reads labeled with either their annotation or "no_annotation"
```
    example
        cat S1_annotated.txt S1_no_annotation_formatted.txt | sort > S1_all.txt
```

### Sequential filtering and blacklisting of reads by annotation class – filter_by_class.py
As this pipeline includes multi-mapping reads, implementing an appropriate treatment of the multi-mappers took significant consideration. We adopted a strategy inspired by TEsmall (O'Neill *et al.*, 2018) that serves as a pipeline to profile small RNAs derived from transposable elements. Ultimately, this approach sequentially assigns reads to gene classes (rRNAs, tRNAs, miRNAs, exons, etc.) then blacklists those reads from other gene classes. The priority order of gene classes can be modified; we performed the analysis using two such orders. First, we classified reads based on priority of the natural abundance of each class (e.g. rRNAs are known to be very prevalent and as such, take high priority). Second, as a very stringent filter, we adjusted the priority order to place the classes of interest (tRNAs, transposable elements, and piRNAs) lowest on the priority list.
```
    example abundance priority order:
        rRNA, snRNA, tRNA, miRNA/hairpin, exon, transposable element, intron, piRNA, no annotation
    example stringent priority order:
        rRNA, snRNA, miRNA/hairpin, exon, intron, piRNA cluster, transposable element, tRNA, no annotation
```
```
    inputs		plaintext file of mapped locations and annotations (*_all.txt)
    outputs	    plaintext file of mapped locations and annotations, filtered by priority (1 file per class)
    parameters	(specified within filter_by_class.py)
        ordered_features		priority list of gene classes (from highest to lowest priority)
    example	
       ordered_features = ["rRNA", "snRNA", "tRNA", "miRNA", "hairpin", "exon", "TE_DNA","TE_SINE", "TE_LINE", "TE_LTR", "TE_Other", "TE_Unknown", "TE_RC", "intron", "piRNA_cluster", "no_annotation"]
```

## Enrichment Calculation - Main Workflow - 4a_simple_enrichment directory

### Normalize multi-mappers (1/n) within each annotation class – normalize_multimappers.py
Classified reads are then 1/n-normalized within each class based on their total read count and the number of loci that they map to.
```
    inputs		plaintext files of mapped locations and their annotations (1 file per class)
    outputs	    plaintext files of mapped locations, annotations, and normalized read count (1 file per class)
```

### Make a composite, normalized annotation file with all post-filtering annotations
```
    inputs		plaintext files of mapped locations, annotations, and normalized read counts (1 file per class)
    output		composite file of mapped locations, annotations, and normalized read counts
    example	
        cat abundance/S1/*.txt > abundance/S1/S1-all_normalized.txt
```

### Calculate read-count pileup by annotation – pileup_by_annotation.py
Each annotation then receives a score corresponding to the total number of in-class normalized reads it received.
```
    input		composite file of mapped locations, annotations, and normalized read counts
    output		list of annotations, normalized read counts, and normalized read counts in rpm
```

### Calculate fold-enrichment – enrich_calc4.py
Libraries are then compared to determine fold enrichment. In our case, six libraries were processed as described above (three input libraries, one background library, and two replicate IP libraries). The fold-enrichment is calculated from annotations that have sufficient input and IP rpm read counts. Thresholds can be set individually (see below). Background rpm read counts are scaled compared to IP rpm read counts based on the amount of background library included in the sequencing submission. Fold enrichment is calculated as (average IP rpm count – background rpm count) / input rpm count. Inclusion of the background sample allows subtraction of RNAs that were pulled-down in the absence of the protein of interest.
```
    inputs	        	    list of annotations, normalized read counts, and normalized read counts in rpm
                            (1 file per library; 6 files in this study)
    output		            list of annotations above read thresholds and corresponding fold-enrichments
    parameters	            (specified within enrich_calc.py)
        bk_scale			global scale factor for background library
                            (based on amount included in library submission)
        rpm_threshold_input minimum average reads per million for the input annotations
                            (default = 0.05; corresponds to ~1 read per input library)
        rpm_threshold_ip	minimum average reads per million for the IP annotations
    				        (default = 0.5; corresponds to ~10 reads per IP library)
```

## Enrichment Calculation - Secondary Workflow - 4b_deseq2 directory

### Build logical gene model – make_logical_gene_list.py
In order to corroborate our main workflow, we implemented another pipeline that employs an alternative gene model that would simultaneously track multi-mapping reads while not immediately normalizing their read counts across the loci that they map to. To accomplish this, a "logical gene" model was created such that multimapping reads were assigned composite annotations. As an example, a multi-mapping read with two mappings (geneA and geneB) would be given a new annotation "geneA;geneB". This new gene model then maintains raw read counts (as required by packages like DESeq2), but could be deconvoluted downstream into the contributing annotations.
```
    inputs		plaintext files of mapped locations and their annotations (1 file per class)
    outputs	    list of logical genes with corresponding raw read counts
    		    list of logical genes with corresponding contributing reads
```

### Create composite files with all logical genes and the read counts from each library – command line and deseq2_prep.py
```
    inputs		lists of logical genes with corresponding raw read counts (1 file per class)
    outputs	    list of all logical genes (per library) and raw read counts (1 file per library)
    	        composite list of all logical genes (across all libraries)
    	        and each library's corresponding raw read count (1 file)	
    example
        cat S1/*.txt | sort > S1/S1-all-logical_genes.txt
        Then run deseq2_prep.py
```

### Determine fold-enrichment with DESeq2 – Run using Galaxy
```
    input		composite list of all logical genes and corresponding read counts
    outputs	    DESeq2 generates several output files. For more information, see the DESeq2
    manual. 
```
The main alignment file is *_output.tabular and contains each annotation, log2 fold changes, and p-value statistics.

### Calculate logical gene enrichment after background subtraction – deseq2_polish.py
DESeq2 output is then analyzed to determine the background-subtracted fold-enrichment for each logical gene. Fold enrichment is calculated as (average IP fold-change [compared to input] – background fold-change [compared to input]). Inclusion of the background sample allows subtraction of RNAs that were pulled-down in the absence of the protein of interest.
```
    inputs		    DESeq2 output for three combinations of files:
                        ip versus input, ip versus background, and background versus input
    output		    annotations, background-subtracted fold-changes [raw, not log2], and p-values
                    (less significant value of [ip versus input; ip versus background])
    parameters	    (specified within deseq2_polish.py)
        bk_scale	global scale factor for background library
                    (based on amount included in library submission; in our case 0.62)
```

### Add no_annotation class – command line
```
example
    sed 's/^no_annotation/no_annotation:no_annotation/' S123456-deseq2_input.txt > S123456-deseq2_input_formatted.txt ; \
    sed 's/^no_annotation/no_annotation:no_annotation/' S123456-deseq2_polished.txt > S123456-deseq2_polished_formatted.txt
```

### Calculate gene enrichment (convert logical genes to canonical annotations) – logical_to_real.py
This script will take the fold enrichment values from the background-subtracted DESeq2 output and apply them to the canonical annotations. Briefly, logical gene annotations that are more unique (i.e. logical genes that have fewer contributing canonical annotations) contribute more strongly to the fold-enrichment calculation for those that are less unique. Ultimately, this results in the reported fold-enrichment values to weight more unique reads to a greater extent than multi-mappers without excluding multi-mappers entirely. Fold-change and p-value filters can be specified.
```
    inputs		DESeq2 output that has been background-subtracted and formatted
    			DESeq2 input file that has raw read counts for each logical gene annotation
    output		plaintext file of annotations and their corresponding fold-change.
    parameters	(specified within logical_to_real.py)
    pval	    adjusted p-value cutoff (default = 1; include all entries)
    enrich	    minimum fold-enrichment (default = 0; include all entries)
```