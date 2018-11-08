#!/bin/bash
#set -e
#set -u
#set -pipefail

#VirSCAN pipeline

#SBATCH -A b2016266
#SBATCH -t 02:00:00 #--qos=short #
#SBATCH -p core
#SBATCH -n 1
#SBATCH -o Result_VirSCAN_20170530.out
#SBATCH -J runVirSCAN_run4

# First add the modules needed for the program(s) your're about to run:
# bioinfo-tools
# samtools
# python
# bowtie

# cd to the directory where the data is (of course, use your own user name):
cd /proj/b2016266/nobackup/Microbiome/VirScan_pipeline/

#specify the input samples file, where the path for each sample file is on new line
sample_info=run4_sample_seq_info.csv

#creates Bash array from the first column
sample_NGI=($(cut -d ',' -f 1 "$sample_info")) #get the first column, NGI id of sample
sample_ID=($(cut -d ',' -f 2 "$sample_info")) #get the first column, ID of sample 

#creating sample path/name
#initialize variable to count samples
#count=0
#for NGI in ${sample_NGI[@]}
#do
#    count=$((count+1))

    #sample path name                                                                                                                         
#    sample_file1="../Data/run4/${NGI}/02-FASTQ/170327_D00410_0396_AHFMJVBCXY/${NGI}_S${count}_L001_R1_001.fastq.gz"
#    sample_file2="../Data/run4/${NGI}/02-FASTQ/170327_D00410_0396_AHFMJVBCXY/${NGI}_S${count}_L002_R1_001.fastq.gz"
    #change format .gz to .bzip2                                                                                                                                                    
#    sample_bz2="../Data/run4/${NGI}/02-FASTQ/170327_D00410_0396_AHFMJVBCXY/${NGI}.fastq.bz2"
#    sample_bam=$(basename $sample_bz2 .fastq.bz2)
#    zcat  ${sample_file1} ${sample_file2} | bzip2 > $sample_bz2
    
    ##ALIGN: decompress the bzipped sequence data, use bowtie to align and output a sorted compressed BAM alignment file   
#    echo "Reads alignment to the reference genome for sample ... ${sample_bam}"
#    bzip2 -dc ${sample_bz2} | bowtie -n 3 -l 30 -e 1000 --tryhard --nomaqround --norc --best --sam --quiet /proj/b2016266/nobackup/Microbiome/SampleData/vir2 - | samtools view -u -| samtools sort -T ${sample_bam}.bam -o ${sample_bam}.bam

    ##COUNT READS    
#    echo "COUNTING READS FOR SAMPLE...${sample_bam}"
    #Index aligned bam file to create .bai file                                                                                         
#    samtools index ${sample_bam}.bam
    #count the number of reads for each reference sequence and convert into a csv file                                                                               
#    samtools idxstats ${sample_bam}.bam | cut -f 1,3 | sed -e '/^\*\t/d' -e '1 i id\tSAMPLE_ID' | tr "\\t" "," > ../Result/run4/count/${sample_bam}.count.csv
    #compress csv to gzip                                                                                                                                 
#    gzip ../Result/run4/count/${sample_bam}.count.csv
    
    #Calculate zero-inflated p-values from counts                                                                                                   
#    echo "calculate zero-inflated p-values from counts for sample .....${sample_bam}"
#    mkdir ../Result/run4/Figures/log_${sample_bam}                                                                                                                              
#    python calc_zipval.py ../Result/run4/count/${sample_bam}.count.csv.gz vir2.170307.count.csv.gz ../Result/run4/Figures/log_${sample_bam} > ../Result/run4/zipval/${sample_bam}.zipval.csv                                           
#    gzip ../Result/run4/zipval/${sample_bam}.zipval.csv                                                                                                   
#done

#initialize sample id count                                                                                                                                                         
#y=1

#mkdir RESULTS/run5/zipval
#for samp in ${sample_NGI[@]}
#do
   
    #Calculate zero-inflated p-values from counts                                                                                                                               
 #   echo "CALCULATING ZERO-INFLATED P-VALUES FOR SAMPLE .....${samp}"                                                                                             
 #   mkdir RESULTS/run5/Figures/log_${samp}_S${y}                                                                                                                                
   
 #   python calc_zipval.py RESULTS/run5/${samp}_S${y}.count.csv.gz vir2.170307.count.csv.gz RESULTS/run5/Figures/log_${samp}_S${y} > RESULTS/run5/zipval/${samp}_S${y}.zipval.csv 
  #  gzip RESULTS/run5/zipval/${samp}_S${y}.zipval.csv                                                                                                     

   # touch results_zipval.txt 
   # echo RESULTS/run5/zipval/${samp}_S${y}.zipval.csv >> results_zipval.txt
    
   # y=$((y+1))
#done

#####

#file containing zipval files info                                                                                                                                                  
result_info=run4_sample_info.txt
#subj=($(cut -f 1 "$result_info")) #subject e.g. 135, K510
samp=($(cut -f 3 "$result_info")) #sample ID
rep1=($(cut -f 1 "$result_info")) #replicate 1
rep2=($(cut -f 2 "$result_info")) #replicate 2

#initialize sample id count                                                                                     
cnt_id=0
x1=1
x2=2
#y=0

#mkdir RESULTS/run5/hits
#touch final_baby_files.txt
for samp in ${samp[@]}
do                                                                                                                                                                        
   # samp_id=$(basename ${sample_ID[${x1}]} _2)
    #Call hits from replicate zero-inflated p-values                                                                                                            
    echo "CALL HITS FROM REPLICATE ZERO-INFLATED P-VALUES FOR SAMPLE .....$samp"
    mkdir ../Result/run4/Figures/log_${rep1[${cnt_id}]}_${rep2[${cnt_id}]}_hits                                                                                                                            
    echo "Repl_1     ../Result/run4/zipval/${rep1[${cnt_id}]}.zipval.csv.gz"
    echo "Repl_2     ../Result/run4/zipval/${rep2[${cnt_id}]}.zipval.csv.gz"
                 
    python call_hits.py ../Result/run4/zipval/${rep1[${cnt_id}]}.zipval.csv.gz ../Result/run4/zipval/${rep2[${cnt_id}]}.zipval.csv.gz ../Result/run4/Figures/log_${rep1[${cnt_id}]}_${rep2[${cnt_id}]}_hits 2.3 > ../Result/run4/hits/${rep1[${cnt_id}]}_${rep2[${cnt_id}]}.zihit.csv                      
    
    gzip ../Result/run4/hits/${rep1[${cnt_id}]}_${rep2[${cnt_id}]}.zihit.csv                                                                                                                          
    #calculate virus scores from hits                                                                                                                                      
    echo "hit_sample        ../Result/run4/hits/${rep1[${cnt_id}]}_${rep2[${cnt_id}]}.zihit.csv  "
    python calc_scores.py ../Result/run4/hits/${rep1[${cnt_id}]}_${rep2[${cnt_id}]}.zihit.csv.gz VIR2.csv.gz vir2.nhits.beads.csv.gz vir2.nhits.samps.csv.gz 'Species' 7 > ../Result/run4/scores/${samp[${cnt_id}]}_ziscore.csv 

    #update counting variable
    cnt_id=$((cnt_id+1))
    x1=$((x1+2))
    x2=$((x2+2))
 #   y=$((y+1))
    
    #running status
    echo "FINISH CALCULATING VIRUS SCORE FOR SAMPLE...$samp"
 #   echo "Seq replicate number ++++++++++++ .........${y}"

    #echo "RESULTS/run5/scores/${sample_NGI[${cnt_id}]}_${sample_NGI[${x1}]}.ziscore.csv">>final_baby_files.txt    
done
######
#calculate the virus score
#file containing zipval files info                                                                                                                                                  
#result_info=run5_sample_info.txt
#subj=($(cut -f 1 "$result_info")) #subject e.g. 135, K510                                                                                                                          
#samp=($(cut -f 2 "$result_info")) #sample ID                                                                                                                                       
#rep1=($(cut -f 3 "$result_info")) #replicate 1                                                                                                                                     
#rep2=($(cut -f 4 "$result_info")) #replicate 2
#week=($(cut -f 5 "$result_info")) 
#days=($(cut -f 6 "$result_info")) #Gestational age 
#k=0
#for i in ${subj[@]}
#do 
 #   echo "CALCULATING SCORE FOR SAMPLE ... RESULTS/run5/hits/${rep1[${k}]}_${rep2[${k}]}.zihit.csv.gz "                                                                            
#    python calc_scores.py RESULTS/run5/hits/${rep1[${k}]}_${rep2[${k}]}.zihit.csv.gz VIR2.csv.gz vir2.nhits.beads.csv.gz vir2.nhits.samps.csv.gz 'Species' 7 > RESULTS/run5/New_scores/${subj[${k}]}_${week[${k}]}_${days[${k}]}_${samp[${k}]}_ziscore.csv                                                                                                          

 #   k=$((k+1))

#done 



####
#baby_info=final_baby.txt
#baby_id=($(cut -f 1 "$baby_info"))
#subj_id=($(cut -f 2 "$baby_info"))

#touch final_info.txt

#for score in ${baby_id[@]}
#do
#    echo "RESULTS/run5/scores/${score}.ziscore.csv" >> final_info.txt
#done
#baby_score_info=final_info.txt
#baby_score=($(cut -f 1 "$baby_score_info"))

#python concat_tables.py ${baby_score[@]} > COMBINED_baby_virus_score.csv




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
#zcat ../Data/run4/Baby/P7805_1096/P7805_1096_S96_L001_R1_001.fastq.gz ../Data/run4/Baby/P7805_1096/P7805_1096_S96_L002_R1_001.fastq.gz | bzip2 > ../Data/run4/Baby/P7805_1096/P7805_1096.fastq.bz2

# Randomly downsampling the number of reads                                                                                                                                  
#bunzip2 ../Data/run4/Baby/P7805_1069/P7805_1069_S69_L001_R1_001.fastq.gz
#seqtk sample -s100 ../Data/run4/Baby/P7805_1069_S69_L001_R1_001.fastq 500000 | bzip2 > ../Data/run4/Baby/P7805_1069/P7805_1069.fastq.bz2

##ALIGN: decompress the bzipped sequence data, use bowtie to align and output a sorted compressed BAM alignment file                                                        
#bzip2 -dc ../Data/run4/Baby/P7805_1096/P7805_1096.fastq.bz2 | bowtie -n 3 -l 30 -e 1000 --tryhard --nomaqround --norc --best --sam --quiet /proj/b2016266/nobackup/Microbiome/SampleData/vir2 - | samtools view -u - | samtools sort -T P7805_1096.bam -o P7805_1096.bam

##COUNT READS                                                                                                                                                               
#Index aligned bam file to create .bai file                                                                                                                                
#samtools index P7805_1096.bam

#count the number of reads for each reference sequence and convert into a csv file                                                                                         
#samtools idxstats P7805_1096.bam | cut -f 1,3 | sed -e '/^\*\t/d' -e '1 i id\tSAMPLE_ID' | tr "\\t" "," > P7805_1096.count.csv
#compress csv to gzip                                                                                                                                                      
#gzip P7805_1096.count.csv

    




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