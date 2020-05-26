# eCLIP Data Processing Pipeline

## Summary & Key Points
This folder has scripts for processing eCLIP next-generation sequencing libraries.
- These scripts accommodate multi-mapping reads and employ a sequential annotation class filtering scheme to assist with normalization. 
- Background is subtracted based on physical measurements made during library construction in addition to read counts upon sequencing.
- The workflow is outlined below.
- Each script is also annotated.
- Several scripts require either user input or modification based on the user's needs.
- **Please also see the eCLIP_parameters_summary.md file for cataloged explanations of inputs, outputs and parameters.**

## Introduction
Enhanced crosslinking and immunoprecipitation (eCLIP) is a molecular biological protocol that allows for identification and characterization of RNA species that are directly bound to a protein of interest. Briefly, this method uses UV crosslinking as a means to covalently link proteins to RNAs within cells, extract these complexes using affinity or immunological pull-down against a protein of interest, then subject the bound nucleic acids to next-generation sequencing.

In the processing of next-generation sequencing reads, many workflows discard sequences that have ambiguous genomic mapping (multi-mappers) in addition to reads that are common contaminants (including ribosomal RNAs and tRNAs). For our particular use case, we were interested in RNA species from many classes that are highly repetitive and abundant, including tRNAs, piRNAs, and retrotransposons.

To accommodate these needs, we implemented the pipeline below to allow for multi-mapping reads and implement additional background subtraction to better characterize the complexes being isolated.

## Notes
- Several of the scripts require input from the user. These instances are annotated in the markdown below.
- Almost all of the scripts provided here are written using Python 3 with the following dependencies: datetime, os, sys, numpy, scipy, math, pandas
- If combining paired-end reads, we recommend using FLASH, which can be downloaded here: https://ccb.jhu.edu/software/FLASH/
- Mapping is performed using STAR (https://github.com/alexdobin/STAR)
- Additional tools that are used in this workflow are bedtools (https://bedtools.readthedocs.io/) and collapse (from the FASTX utilities in Galaxy; https://usegalaxy.org/)
- DESeq2 (https://bioconductor.org/packages/release/bioc/html/DESeq2.html), which is only used for an alternative workflow is run in R. An example script is provided.
- Suggested file architecture and directory names are listed, but these can be adjusted.  Please note the if you want to customize these, changes may need to be made within the scripts as well.

## Overview
The pipeline is broken into several main sections:

0. Setup
1. Preliminary read processing
2. Mapping and cleanup
3. Annotation and filtering
4. Calculation of enrichment scores<BR>
**Note: Two methods are provided for enrichment score calculation, with the first being the main one used in our workflow.*

## Setup
1. Prepare, sequence, and demultiplex all sequencing libraries such that each fastq file corresponds to a single sample.  For paired-end libraries, each read should be in its own fastq file.
2. For ease in this workflow, the following sample designations are recommended:
	
	| Sample abbreviation | Description |
	| --- | --- |
	| S1 | Non-crosslinked input |
	| S2 | Crosslinked input (replicate 1) |
	| S3 | Crosslinked input (replicate 2) |
	| S4 | Non-crosslinked IP (background) |
	| S5 | Crosslinked IP (replicate 1) |
	| S6 | Crosslinked IP (replicate 2) |

3. Read annotations should be in tab-delimited bed format:
	```
	chromosome	start	stop	annotation_class:annotation	score (optional)	strand
	``` 

## Part 1 - Preliminary read processing
**Note: The beginning of this pipeline works off of a single fastq file. This can come directly from the sequencing of a single read, or, if paired-end reads are overlapping (our use case), they can be combined using FLASH.*

0. Set up directory architecture<BR>
	```mkdir 1_raw_data 2_flash 3_adapter_trimmed 4_remove_pcr_duplicates 5_reverse_complement 6_collapsed```
1. Place raw .fastq files in the 1_raw_data directory

2. If needed, join paired-ends using FLASH.  Example:<BR>
	```FLASH-1.2.11/flash -m 25 -M 175 -O -d 2_flash/S1 -o S1 1_raw_data/S1_R1.fastq 1_raw_data/S1_R2.fastq 2>&1 | tee 2_flash/S1/S1.log```

3. Remove adapters and extraneous sequences by comparing the FLASH input and output. 
    The script will prompt for the input and output file names.  Full paths are recommended. Output can be directed to the 3_adapter_trimmed folder.<BR>
    ```python 1_read_cleanup/flash_trim.py```

    Count the total number of adapter-trimmed reads and the number of adapter dimers.  This will be needed for library normalization in the final enrichment step.<BR>
	```python 1_read_cleanup/count_adapters_only.py```

4. Unique Molecular Identifier (UMI) collapse
    Remove PCR duplicates by comparing UMIs within each file.  The script will prompt to specify the file names for the inputs. Output can be directed to the 4_remove_pcr_duplicates folder.<BR>
    **Note: The library construction is slightly different for the input samples versus the CLIP samples.  Make sure to run the correct version of umi_collapse.py.*<BR>
	```python 1_read_cleanup/umi_collapse.py```<BR>
	```python 1_read_cleanup/umi_collapse_inputs.py```

5. Reverse complement
    The script will prompt to specify the file names for the inputs. Output can be directed to the 5_reverse_complement folder.<BR>
    **Note: eCLIP libraries are stranded with the reverse complement corresponding to sense strand of the input.*<BR>
	```python 1_read_cleanup/rev_comp.py```

6. Collapse reads
    This is performed to greatly speed up mapping. We used the "collapse" tool in Galaxy.  Ouput can be directed to the 6_collapsed folder.  The output will be a FASTA file with the format:
	```
	>readnumber-count
	SEQUENCE
	```

## Part 2 - Mapping
0. Set up directory architecture<BR>
	```mkdir 7_star 8_bamtobed```
    
    Generate genomes (mm9, mm10, and dm6 examples are below, but please be aware that you need to download and name the genome files accordingly)

    **mm9**
	```
	STAR_2.5.2b 
	--runThreadN 12
	--runMode genomeGenerate
	--genomeDir mm9
	--genomeFastaFiles NCBIM37.genome.fa
	--sjdbGTFfile gencode.vM1.annotation.gtf
	--sjdbOverhang 100
	--limitGenomeGenerateRAM 100000000000
	```

    **mm10**
	```
	STAR_2.5.2b 
	--runThreadN 12
	--runMode genomeGenerate
	--genomeDir genomeGRCm38_M15
	--genomeFastaFiles GRCm38.primary_assembly.genome.fa
	--sjdbGTFfile gencode.vM15.annotation.spikeins.gtf
	--sjdbOverhang 100
	--limitGenomeGenerateRAM 100000000000
	```

    **dm6**
	```
	STAR_2.5.2b 
	--runThreadN 12
	--runMode genomeGenerate
	--genomeDir dm6
	--genomeFastaFiles dmel-all-chromosome-r6.27.fasta
	--sjdbGTFfile dmel-all-r6.27.gtf
	--sjdbOverhang 100
	--limitGenomeGenerateRAM 100000000000
	```

1. Run STAR
    These parameters allow for up to 50 multimapping loci.
    Output should be place in the 7_star directory.
	```
	STAR_2.5.2b
	--genomeDir genomeGRCm38_M15
	--runThreadN 12
	--readFilesIn  S1.extendedFrags_trimmed_unique_rc_collapsed.fasta 
	--outSAMtype BAM SortedByCoordinate
	--quantMode GeneCounts TranscriptomeSAM
	--outSAMattributes NH HI nM AS NM MD
	--outFilterMultimapNmax 50
	--outFilterMismatchNoverReadLmax 0.05
	--outFilterMismatchNmax 999
	```

2. Generate .bed files from STAR output
    This will be need to be run for each sample.  Output should be placed in the 8_bamtobed folder.<BR></BR>
    <I>*Note: uniqued files seem to have about 0.1% fewer entries ~10K entries for 10M read library.<BR>
    *Note: It is a good idea to check the formatting of the .bed files.  Locus formatting should match that in the main annotation .bed that will be used.  Depending on the STAR input files, the STAR output may or may not have some common formatting features (such as "chr" beginning each entry).  Check for consistency and include an awk command below if minor formatting changes are required.<BR>
    *Note: For the mm9 and mm10 assembly annotations, the formatting is consistent so no addition awk command is needed. It is needed, however, for the dm6 annotations.</I><BR>
   
    **Example - mm10:**<BR>
	```bedtools bamtobed -i S1/Aligned.sortedByCoord.out.bam > S1.bed; sort -k1,1 -k2,2n S1.bed | uniq > S1_unique.bed```

    **Example - dm6:**<BR>
	```bedtools bamtobed -i S1/Aligned.sortedByCoord.out.bam > S1_temp.bed; awk '{print "chr" $0}' S1_temp.bed > S1.bed; rm S1_temp.bed; sort -k1,1 -k2,2n S1.bed | uniq > S1_unique.bed```

3. Clean up the mapping to remove "multi-aligners" that are not true "multi-mappers."<BR></BR>
    ```python 2_mapping_cleanup/unique_bed_tidy.py```

## Part 3a - Annotation
Annotation is performed with bedtools. The -f 1 option requires complete intersection "a" within "b"; -s enforces strandedness. The path to the annotation .bed file must be specified for file "b".

1.  Annotate<BR>
    ```bedtools intersect -f 1 -wa -wb -s -a 8_bamtobed/S1_unique_tidy.bed -b mm10_annotations/mm10_annotations.bed > 9_annotated/S1_annotated.txt```

2. Format the no_annotation file so that is consistent with the annotated file list<BR>
	```python 3_annotation_and_filtering/format_no_annotation.py```

3.  Create a composite file with all reads labeled with either their annotation or "no annotation"<BR>
	```cat S1_annotated.txt S1_no_annotation_formatted.txt | sort -k1,1 -k2,2n > S1_all.txt```

## Part 3b - Filtering
Filtering is performed to ascribe multimapping reads to a particular class.
Most usually, the filtering order is based on transcript abundance, but the order can be modified with the scripts.<BR>

0. Set up directories. Subdirectories can be useful if you want to try multiple filtering orders (see below and the filter_by_class.py script).<BR>
	```mkdir 10_filtered```

1. Determine ALL of the annotation classes in the main annotation .bed file<BR>
	```python 3_annotation_and_filtering/get_annotation_classes.py _path_to_file_name_```

2. **CRITICAL STEP!** - Determine the preferred filtering order.
    - Using all of the classes from get_annotation_classes, make an ordered list.
    - Highest priority classes go first.
    - Once a read assigns to a class, it is blacklisted from all lower priority classes. Thus, the order is VERY IMPORTANT to the analysis.
    - Make this a list and add it to filter_by_class.py file.
    **Example:**<BR>
	```ordered_features = ["rRNA", "snoRNA", "scaRNA", "snRNA", "ribozyme", "lncRNA", "tRNA", "pseudogene_tRNA", "mttRNA", "miRNA", "pri-miRNA", "exon", "TE_DNA", "TE_SINE", "TE_LINE", "TE_Satellite", "TE_LTR", "TE_RNA", "TE_RC", "TE_Other", "TE_Unknown", "intron", "pre-mRNA", "piRNA_cluster", "sncRNA", "miscRNA", "pseudogene", "asRNA", "no_annotation"]```

3. Perform filtering. This step is somewhat slow, so it may be best to run it on a server. Output will go to the 10_filtered folder<BR>
	```python 3_annotation_and_filtering/filter_by_class.py```

## Part 4 - Calculate Enrichment
### Option 1 - Simple Enrichment
This option will normalize multimappers within each class, calculate annotation pileups, subtract background, and calculate fold-changes.<BR>
	
0. Set up directories<BR>
	```mkdir 11a_normalized 12a_pileup 13a_enrichment```

1. Normalize multi-mappers within each class (1/n-normalization).<BR>
	```python 4a_simple_enrichment/normalize_multimappers.py```

2. It is recommend to sort the output of this script into individual subdirectories (S1/ S2/ ... S6/)

3. Make a composite file for each of the libraries:<BR>
    ```cat S1/*.txt > S1/S1-all_normalized.txt```

4. Generate the summary:<BR>
    ```python 4a_simple_enrichment/reads_per_annotation.py```

5. Calculate the annotation pileups for each library.<BR>
	```python 4a_simple_enrichment/pileup_by_annotation.py```

6. Move the output files to the 12a_pileup directory.

7. Combine pileups into a single file. This is useful for Option 2 (below), but is convenient to do here.<BR>
    ```python 4a_simple_enrichment/all_annotation_pileup.py```

8. Calculate the enrichment for each library
    - Be sure to set the bg_globscale factors. To set these correctly, the yields of each IP library (concentration * volume) and their corresponding read counts are required.*
	- Concentrations should come from Bioanalyzer quantification of
		non-adapter dimers.
	- Yields should come from elution volumes.
	- Read counts should be those from after adapter-trimming (subtracting the number of adapter dimers).
	- Both of these values can be obtained from the output of count_adapters_only.py (adapter_dimers.log); generated above.
	- Typical values for bg_globscale are between 0.5 and 0.9 but may be outside this range. Values closer to 0 indicate low amounts of DNA in the background IP compared to the experimental samples.  Values higher than 1 indicate that the background had more material than the experimental pull-down (suspicious of bad library prep).
	- Low read-count thresholds can also be adjusted in the script.  A value which gives 1 read per input library and 10 reads per pull-down library seem reasonable.<BR>
	```python 4a_simple_enrichment/enrich_calc4.py```

4. File clean-up
    - Move the enrichment file outputs to the 13a_enrichment directory.
    - Make the final enrichment file consistent and sorted by enrichment score<BR>
        ```sort -nk8 -r S123456-enrich_per_annotation.txt > S123456-enrich_per_annotation_sorted.txt; \```<BR>
	```sed -e 's/no_annotation/no_annotation:no_annotation/' -e 's/:/        /1' S123456-enrich_per_annotation_sorted.txt | awk '{print $2 "\t" $1 "\t" $9}' > S123456-enrich_simple.txt```

### Option 2 - DESeq2
This option relies on building an alternative gene model to track multimappers. The advantage of this pipeline is that it provides more statistics using the DESeq2 pipeline.  The disadvantage is that deconvolution of the alternative gene model likely introduces more approximations into the data.

Code is provided, however, due to the complexities of this option, it is recommended only for expert users with working familiarity with DESeq2 and as such, not as much documentation is provided.

There are many ways to customize this pipeline and we have tried several variations (including rounding of normalized reads using the alternative gene model, rounding of normalized reads using the conventional gene model, and different means of subtracting background.) The workflow below includes all inputs and the alternative gene model.

#### Explanation:
Since DESeq2 does not support fractional reads, a new gene model is built.
- Multimappers in a given class are given new annotation in the form of "logical_genes"
- Example: 
    - read1 maps to geneA and geneB; 		read2 maps to geneA; 		read3 maps to geneB
    - Three logical gene annotations will be made:
        -  read1 --> geneA;geneB
        -  read2 --> geneA
        -  read3 --> geneB
    - Ultimately, this allows multimappers to be annotated in multiple genes while preserving raw read counts (important for precision determination in DESeq2).
    - Additionally, lists of logical_genes with corresponding read_ids (e.g. S1_readid-count) are output.
    - After enrichment analysis, the multimappers can be 1/n normalized and distributed to the 
	appropriate annotation.

1. Build the alternative gene model (i.e. "logical genes")<BR>
<I>*Note: Make sure all classes are accounted for.  This should be the same list as used in filter_by_class.py (above).  Class order does not matter for this script.</I><BR>
	```python 4b_deseq2/make_logical_gene_list.py```
Move the outputs to a directory 11b_logical_genes.

2. Create composite files will all logical_genes and read_ids
**Example:**<BR>
	```cat S1/*.txt | sort > S1/S1-logical_abundance_lookup.txt```
	```cat S1/*.txt | sort > S1/S1-all-logical_genes.txt```

3. Run DESeq2 - Preparation
A prep script was written that takes *-all-logical_genes.txt files for each sample and builds a composite file: [logical_gene, S1_value, S2_value, ...]<BR>
	```python 4b_deseq2/deseq2_prep.py```

4. Run DESeq2
    1. Launch the script run_deseq2_eclip_v1.Rmd in R studio. <B>This script assumes the column order as follows: bg_input, ip1_input, ip2_input, bg_ip, ip1, ip2</B>
    2. Set the folder and working directory.<BR></BR>
	**DESeq2 will be run three times:**
    	1. ip versus all inputs
    	2. bg_ip versus all inputs
    	3. experimental ips versus background ip

5. Calculate fold-enrichments and extract statistics
Background-subtracted fold enrichments are calculated by comparing the various DESeq2 outputs.
The output file will have the following tab-delimited columns:
    ```
    gene name (or logical gene name)    fold enrichment     p-value[adjusted]
    ```
    <I>*Note: Make sure to adjust the background scaling factors at the top of the script.  These should be the same as those used in the simple enrichment pipeline.<BR>
    *Note: DESeq2 output includes log2-fold enrichments.  These are converted to fold enrichments (not log2).<BR>
    *Note: p-values will be the less significant value of [ip versus input; ip versus background])</I><BR>
    ```python 4b_deseq2/deseq2_polish.py```

6. File cleanup - Add no_annotation class<BR>
	```sed 's/^no_annotation/no_annotation:no_annotation/' S123456-deseq2_input.txt > S123456-deseq2_input_formatted.txt```<BR>
	```sed 's/^no_annotation/no_annotation:no_annotation/' S123456-deseq2_polished.txt > S123456-deseq2_polished_formatted.txt```

7. Convert from logical_genes to normal annotations
Note: p-value thresholds can be set in the logical_to_real.py script.<BR>
	```python 4b_deseq2/logical_to_real.py```
