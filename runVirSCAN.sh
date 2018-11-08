#!/bin/bash
#set -e
#set -u
#set -pipefail

# VirSCAN pipeline

#SBATCH -A b2016266
#SBATCH -t 01:00:00 #--qos=short #
#SBATCH -p core
#SBATCH -n 1
#SBATCH -o Result_VirSCAN_20170330_baby96.out
#SBATCH -J runVirSCAN_baby96

# First add the modules needed for the program(s) your're about to run:
# bioinfo-tools
# samtools
# python
# bowtie

# cd to the directory where the data is (of course, use your own user name):
cd /proj/b2016266/nobackup/Microbiome/VirScan_pipeline/

#specify the input samples file, where the path for each sample file is on new line
#sample_info=run5_sample_info.txt

#creates Bash array from the first column
#sample_NGI=($(cut -f 1 "$sample_info")) #get the first column, NGI id of sample
#sample_ID=($(cut -f 2 "$sample_info")) #get the first column, ID of sample 

#creating sample path/name
#initialize variable to count samples
#count=0
#for NGI in ${sample_NGI[@]}
#do
 #   count=$((count+1))
    #sample path name
 #   sample_file="../Data/run5/${NGI}/02-FASTQ/170421_D00410_0403_AHJHH2BCXY/${NGI}_S${count}_L001_R1_001.fastq.gz"
    
    #change format .gz to .bzip2
  #  sample_bz2="../Data/run5/${NGI}/02-FASTQ/170421_D00410_0403_AHJHH2BCXY/${NGI}_S${count}_L001_R1_001.fastq.bz2"
   # bzip2 ${sample_file} > $sample_bz2

    #sample_bam=(basename $sample_file _L001_R1_001.fastq.bz2)
    ##ALIGN: decompress the bzipped sequence data, use bowtie to align and output a sorted compressed BAM alignment file   
    #echo "Reads alignment to the reference genome for sample ... ${sample_bam}"
   # bzip2 -dc ${sample_bz2} | bowtie -n 3 -l 30 -e 1000 --tryhard --nomaqround --norc --best --sam --quiet /proj/b2016266/nobackup/Microbiome/SampleData/vir2 - | samtools view -u -| samtools sort -T ${sample_bam}.bam -o ${sample_bam}.bam

    ##COUNT READS    
   # echo "Counting reads for sample...${sample_bam}"
    #Index aligned bam file to create .bai file                                                                                         
   # samtools index ${sample_bam}.bam
    #count the number of reads for each reference sequence and convert into a csv file                                                                               
    #samtools idxstats ${sample_bam}.bam | cut -f 1,3 | sed -e '/^\*\t/d' -e '1 i id\tSAMPLE_ID' | tr "\\t" "," > RESULTS/run5/${sample_bam}.count.csv
    #compress csv to gzip                                                                                                                                                           
    #gzip RESULTS/run5/${sample_bam}.count.csv

    #Calculate zero-inflated p-values from counts                                                                                                                 
    #echo "calculate zero-inflated p-values from counts for sample .....${sample_bam}"
    #mkdir RESULTS/run5/log_${smaple_bam}                                                                                                                              
    #python calc_zipval.py RESULTS/run5/${smaple_bam}.count.csv.gz vir2.170307.count.csv.gz RESULTS/run5/log_${smaple_bam} > RESULTS/run5/${sample_bam}.zipval.csv                                           
    #gzip RESULTS/run5/${sample_bam}.zipval.csv                                                                                                                              

    #print results paths in file
    #touch results_zipval.txt
    #result_info=results_zipval.txt

    #echo RESULTS/run5/${sample_bam}.zipval.csv >> results_zipval.txt
#done

#file containing zipval files info
#result_info=results_zipval.txt 

    #Call hits from replicate zero-inflated p-values                                                                                                            
    #mkdir log_Control_hits                                                                                                                                                       
    #python call_hits.py P6412_control_L1.zipval.csv.gz P6412_control_L2.zipval.csv.gz log_Control_hits 2.3 > P6412_control.zihit.csv                      
    #gzip P6412_control.zihit.csv                                                                                                                          
    
    #calculate virus scores from hits                                                                                                                                               
    #python calc_scores.py ZIHITS.zihit.csv.gz METADATA 



##############################
#result_file="$(basename ${sample_files[0]} .S53_L001_R1_001.fastq.gz)"

#if [ ${#sample_files[@]} -eq 2 ]
#then
#    zcat ${sample_files[0]} ${sample_files[1]} | bzip2>${result_file}_input.fastq.bz2 
#    exit 1
#else
#    zcat ${sample_files[0]} ${sample_files[1]} | bzip2>${result_file}_repl_1.fastq.bz2
#    zcat ${sample_files[2]} ${sample_files[3]} | bzip2>${result_file}_repl_2.fastq.bz2
#fi

########## 
#concatenate fastq files from two lanes
zcat ../Data/run4/Baby/P7805_1096/P7805_1096_S96_L001_R1_001.fastq.gz ../Data/run4/Baby/P7805_1096/P7805_1096_S96_L002_R1_001.fastq.gz | bzip2 > ../Data/run4/Baby/P7805_1096/P7805_1096.fastq.bz2

# Randomly downsampling the number of reads                                                                                                                                  
#bunzip2 ../Data/run4/Baby/P7805_1069/P7805_1069_S69_L001_R1_001.fastq.gz
#seqtk sample -s100 ../Data/run4/Baby/P7805_1069_S69_L001_R1_001.fastq 500000 | bzip2 > ../Data/run4/Baby/P7805_1069/P7805_1069.fastq.bz2

##ALIGN: decompress the bzipped sequence data, use bowtie to align and output a sorted compressed BAM alignment file                                                        
bzip2 -dc ../Data/run4/Baby/P7805_1096/P7805_1096.fastq.bz2 | bowtie -n 3 -l 30 -e 1000 --tryhard --nomaqround --norc --best --sam --quiet /proj/b2016266/nobackup/Microbiome/SampleData/vir2 - | samtools view -u - | samtools sort -T P7805_1096.bam -o P7805_1096.bam

##COUNT READS                                                                                                                                                               
#Index aligned bam file to create .bai file                                                                                                                                
samtools index P7805_1096.bam

#count the number of reads for each reference sequence and convert into a csv file                                                                                         
samtools idxstats P7805_1096.bam | cut -f 1,3 | sed -e '/^\*\t/d' -e '1 i id\tSAMPLE_ID' | tr "\\t" "," > P7805_1096.count.csv
#compress csv to gzip                                                                                                                                                      
gzip P7805_1096.count.csv

    




# Randomly downsampling the number of reads
#bunzip2 /proj/b2016266/INBOX/P.Brodin_16_02/P6412_beads.fastq.bz2 

#seqtk sample -s100 /proj/b2016266/INBOX/P.Brodin_16_02/P6412_1002/161031_BH5GW3BCXY/P6412_control_L2.fastq 750000 > /proj/b2016266/INBOX/P.Brodin_16_02/P6412_1002/161031_BH5GW3BCXY/P6412_control_L2_750k.fastq

#bzip2 /proj/b2016266/INBOX/P.Brodin_16_02/P6412_1002/161031_BH5GW3BCXY/P6412_control_L2_750k.fastq


## PART A: Align sequencing data to reference
#Build Bowtie index


##ALIGN: decompress the bzipped sequence data, use bowtie to align and output a sorted compressed BAM alignment file
#bzip2 -dc ../Data/VZ/P7466_1009/P7466_1009_13_vz_repl_2.fastq.bz2 | bowtie -n 3 -l 30 -e 1000 --tryhard --nomaqround --norc --best --sam --quiet /proj/b2016266/nobackup/Microbiome/SampleData/vir2 - | samtools view -u - | samtools sort -T P7466_1009_13_vz_repl_2.bam -o P7466_1009_13_vz_repl_2.bam

##COUNT READS
#Index aligned bam file to create .bai file
#samtools index P7466_1009_13_vz_repl_2.bam 

#count the number of reads for each reference sequence and convert into a csv file
#samtools idxstats P7466_1009_13_vz_repl_2.bam | cut -f 1,3 | sed -e '/^\*\t/d' -e '1 i id\tSAMPLE_ID' | tr "\\t" "," > P7466_1009_13_vz_repl_2.count.csv
#compress csv to gzip
#gzip P7466_1009_13_vz_repl_2.count.csv

#Calculate zero-inflated p-values from counts
#mkdir ../SampleData/log_SABES_2B_A5
#python calc_zipval.py SABES_2B_A5.count.csv.gz vir2.160114.count.csv.gz ../log_SABES_2B_A5 > SABES_2B_A5.zipval.csv

#gzip SABES_2B_A5.zipval.csv 

#Call hits from replicate zero-inflated p-values
#mkdir log_Control_hits
#python call_hits.py P6412_control_L1.zipval.csv.gz P6412_control_L2.zipval.csv.gz log_Control_hits 2.3 > P6412_control.zihit.csv

#gzip P6412_control.zihit.csv

#calculate virus scores from hits
#python calc_scores.py ZIHITS.zihit.csv.gz METADATA 