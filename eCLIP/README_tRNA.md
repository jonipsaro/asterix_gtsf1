# tRNA characterization

## Citation
>Asterix/Gtsf1 links tRNAs and piRNA silencing of retrotransposons<BR>
>Jonathan J. Ipsaro, Paul A. O'Brien, Shibani Bhattacharya, Arthur G. Palmer 3rd, Leemor Joshua-Tor<BR>
>Cell Rep. 2021 Mar 30;34(13):108914. doi: 10.1016/j.celrep.2021.108914.<BR>
>PMID: 33789107

## Setup
Begin by running the eCLIP processing pipeline at least until the step of generating the 1/n-normalized files.

## Generate overall tRNA information
**Note: This script will specifically look at the tRNA_normalized.txt file. If you only have a composite file, first extract the tRNA annotations*
1. Adjust any of the script parameters, including the length of the model tRNA.
2. If there are particular sites of interest, they can be specified in the get_sites variable. If specified, this will generate a bed file listing the reads that have the specified 5prime site.
	```python 5_tRNA_analysis/tRNA_by_position.py```
3. Final processing such as background subtraction can be carried out as in the main eCLIP workflow and should perform the same global scale correction using yields and reads.
