# Maternal anti-viral antibodies in human newborns

Here are scripts to perform analysis of systematic viral epitope scanning (VirScan). 
The principles behind the codes were described in [Xu et al. (2015)] (http://doi.org/10.1126/science.aaa0698)

The scripts have been used to obtain the results published
“The repertoire of maternal anti-viral antibodies in human newborns”
Christian Pou, Dieudonné Nkulikiyimfura, Ewa Henckel, Axel Olin, Tadepally Lakshmikanth, Jaromir Mikes, Jun Wang, Yang Chen, Anna-Karin Bernhardsson, Anna Gustafsson, Kajsa Bohlin and Petter Brodin. (doi: 10.1038/s41591-019-0392-8)

# Dependencies
- Bowtie/1.2.0
- Samtools/1.1
- Python/3.5.4 (numpy/1.11.3, pandas/0.17.1, matplotlib/2.1.2)

# Repo description
`db`       <-- db directory with the input virus library  
`src/`     <-- src directory contains python scripts  
`index/`   <-- index directory with an index for the reference sequences  

# Documentation for each step of data analysis pipeline

The scripts available in this repo are all located in `src/`. The usage descriptions are provided in this readme

## Align sequencing data to reference

Use Bowtie to align the reads to a reference sequence.  
Before you can align, you need to build an index for the reference sequences.   
`bowtie-build REFERENCE.fasta REFERENCE_OUTPUT_NAME`  

This has been done for you and an index for the reference sequences is found in the index directory.  

`bzip2 -dc path/to/FASTQ_SEQUENCES.fastq.bz2 | bowtie -n 3 -l 30 -e 1000 --tryhard --nomaqround --norc --best --sam --quiet path/to/REFERENCE_OUTPUT_NAME - | samtools-1.1/bin/samtools view -u - | samtools-1.1/bin/samtools sort -T BAM_OUTPUT_NAME.bam - -o $@`  

where  
* `path/to/FASTQ_SEQUENCES.fastq.bz2` is the path to the FASTQ file containing the sequencing reads
*  `path/to/REFERENCE_OUTPUT_NAME` is the path to bowtie index (without the `.1.ebwt` suffix)
*  `BAM_OUTPUT_NAME.bam` is the name of the bam file to which you will save the alignment results

## Count reads

Use samtools `idxstats` command to count the number of reads for each sequence. 
Prior to that, it is necessary to index the bam file:  
 
`samtools index BAM_OUTPUT_NAME.bam`

Then, count the number of reads for each reference sequence and converts it into a csv file.  

`samtool idxstats BAM_OUTPUT_NAME.bam | cut -f 1,3 | sed -e '/^\*\t/d' -e '1 i id\tSAMPLE_ID' | tr "\\t" "," >COUNT_FILE.csv`  

Where `SAMPLE_ID` is the id of the sample and `COUNT_FILE.csv` is the name of the count file.  

To save space, compress the csv files with gzip:

`gzip COUNT_FILE.csv`

## Calculate zero-inflated p-values from counts

The python script `calc_zipval.py` will calculate zero-inflated p-values for a set of output counts based on a set of input counts.

`python calc_zipval.py OUTPUT.count.csv.gz INPUT.count.csv.gz log_directory >OUTPUT.zipval.csv`  

*  `OUTPUT.count.csv` is the output read count
*  `INTPUT.count.csv` is the input read count
*  `log_directory` is a directory where the script will save several plots showing model fits
*  `OUTPUT.zipval.csv` is the resulting zero-inflated p-values. 

## Call hits from replicate zero-inflated p-values

The python `script call_hits.py` will call hits based on replicate zero-inflated p-values.   

`python call_hits.py REPLICATE1.zipval.csv.gz REPLICATE2.zipval.csv.gz THRESHOLD log_directory >OUTPUT.zihit.csv`

*  `REPLICATE1.zipval.csv.gz` and `REPLICATE2.zipval.csv.gz` are the two replicate zero-inflated pvalue files
*  `THRESHOLD` is the threshold zero-inflated p-value for calling hits (2.3 for VirScan)
*  `log_directory` is a directory where the script will save plots showing correlation between the replicates
*  `OUTPUT.zihit.csv` is the resulting hits. First column is oligo id, second is True/False

## Calculate virus scores from hits

The python script calc_scores.py calculates virus scores using the maximum parsimony approach. 

`python calc_scores.py ZIHITS.zihit.csv.gz METADATA.csv.gz NHITS.BEADS.csv.gz NHITS.SAMPS.csv.gz GROUPING_LEVEL EPITOPE_LEN >OUTPUT.ziscore.spp.csv`

* ` ZIHITS.zihit.csv.gz` is the gzipped results of call_hits.py
*  `METADATA.csv.gz` is the file containing the metadata for the virus library
*  `NHITS.BEADS.csv.gz` is a two column gzipped csv file, column 1 is oligo id, column 2 is the number of beads samples in which that oligo was a hit
*  `NHITS.SAMPS.csv.gz` is a two column gzipped csv file, column 1 is oligo id, column 2 is the number of non-beads samples in which that oligo was a hit
*  `GROUPING_LEVEL` can be Species or Organism, depending on 
*  `EPITOPE_LEN` should be 7, the length of a linear epitope
*  `OUTPUT.ziscore.spp.csv` is a two column csv file.

## Combine multiple cvss into one table

The `concat_tables.py` script combines multiple csv files into one.

