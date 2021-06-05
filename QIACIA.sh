#!/bin/bash -l

#SBATCH -A sens2018120
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 7-00:00:00
#SBATCH -J QIACIA
#0.1: Inputs, Constants, and Variables:
#0.1.1: Color-constants:
RED='\e[31m'
NC='\e[0m'

printf "\n${RED}1. Write the sample without the suffix (.sra).\n2. Write the reference without the suffix (.fa/.fasta)\n3. Remember to be in the correct directory where your samples are${NC}\n\n" &
##0.1.2: Input
#wait; read -p "Sample: " SAMPLE 
#wait; read -p "Reference: " RFE
SAMPLE=ERS1066752 #Need to change this variable for each sample.
RFE=HS_concat

##0.1.3: Constants:
RF=$RFE.fa
#Rackham
#SRA=/scratch/gusalaba/ncbi/sra
#TMP=/scratch/gusalaba/ncbi/tmp
#FASTQ=/scratch/gusalaba/ncbi/fastq
#QC=/scratch/gusalaba/ncbi/qualitycheck
#REF=/scratch/gusalaba/ncbi/refseq
#BAM=/scratch/gusalaba/ncbi/bam
#ALN=/scratch/gusalaba/ncbi/alignment
#Bianca
SRA=/proj/nobackup/sens2018120/wharf/gusalaba/gusalaba-sens2018120/data/sra
TMP=/proj/nobackup/sens2018120/wharf/gusalaba/gusalaba-sens2018120/data/tmp
FASTQ=/proj/nobackup/sens2018120/wharf/gusalaba/gusalaba-sens2018120/data/fastq/WGS_Matched_NPC
QC=/proj/nobackup/sens2018120/wharf/gusalaba/gusalaba-sens2018120/data/qualitycheck
REF=/proj/nobackup/sens2018120/wharf/gusalaba/gusalaba-sens2018120/data/reference-genome
BAM=/proj/nobackup/sens2018120/wharf/gusalaba/gusalaba-sens2018120/data/bam
ALN=/proj/nobackup/sens2018120/wharf/gusalaba/gusalaba-sens2018120/data/alignment
REFCHECKBWA=$REF/$RF.sa
REFCHECKSAM=$REF/$RF.fai
REFCHECKDIC=$REF/$RFE.dict
wait

##0.1.4: Variables:
FASTQC_Flag=0
Index_Ref_BWA_Flag=0
Index_Ref_SAM_Flag=0
Create_Dict_Flag=0
Aligning_Flag=0
Conversion_Flag=0
Indexing_Alignment_Flag=0
RTC_Flag=0
Indel_Flag=0
Mark_Dup_Flag=0
Indexing_Alignment2_Flag=0
wait

printf "\n${RED}Initiate the QIACIA script: ${NC}\n" 
#1.1: Setting the module environment:
printf "\n${RED}1: Initiate setting the environment.${NC}\n 1.1: Initiate module load.\n"
module load bioinfo-tools FastQC bwa samtools picard GATK && module list #Fix how to load GATK as a module.
printf "Module load: Done\n" &
wait

#1.2: Setting the directories:
##1.2.1: QualityControl:
printf "\n1.2.1: Making a Quality Control Directory:"
cd $QC
if [ -d $SAMPLE ] ; then printf " Directory already exists.\n"
elif [ ! -d $SAMPLE ] ; then mkdir $SAMPLE 
printf " Done.\n"
fi
wait

##1.2.2: BAM:
printf "\n1.2.2: Making a BAM Directory:" 
cd $BAM
if [ -d $SAMPLE ] ; then printf " Directory already exists.\n"
elif [ ! -d $SAMPLE ] ; then mkdir $SAMPLE
printf " Done.\n"
fi
wait

#2: Checking the quality of the reads through FastQC:
printf "\n${RED}2: Checking Quality: ${NC}"
if [ -f $QC/$SAMPLE/$SAMPLE.sra_1_fastqc.html ] ; then printf "Quality Control already exists.\n"
elif [ ! -f $QC/$SAMPLE/$SAMPLE.sra* ] ; then cd $FASTQ && printf "\n\n${RED}Initiate FastQC${NC}\n" && fastqc $SAMPLE.sra_1.fastq $SAMPLE.sra_2.fastq -o $QC/$SAMPLE/ -d $TMP && printf "\nQuality Check: Done\n"
fi
wait

#3: Aligning:
printf "\n${RED}3: Initiate Aligning Stage\n${NC}"
##3.1: Indexing:
###3.1.1: Indexing with BWA index:
if [ -f $REFCHECKBWA ] ; then printf "\n3.1.1: Reference Already Indexed with BWA, proceeding to the next stage.\n" && Index_Ref_BWA_Flag=1
elif [ ! -f $REFCHECKBWA ] ; then printf "\n${RED}3.1: Initiate Indexing of Reference Sequence with BWA.${NC}\n" && cd $REF && bwa index -a bwtsw $RF && Index_Ref_BWA_Flag=1 || Index_Ref_BWA_Flag=2 && printf "\nIndexing of reference sequence with BWA: Done\n" 
fi
wait

###3.1.2: Indexing with samtools:
if [ -f $REFCHECKSAM ] ; then printf "\n3.1.2: Reference Already Indexed with Samtools, proceeding to the next stage.\n" && Index_Ref_SAM_Flag=1
elif [ ! -f $REFCHECKSAM ] && [ $Index_Ref_BWA_Flag == 1 ] ; then printf "\n${RED}3.2: Initiate Indexing of Reference Sequence with Samtools.${NC}\n" && cd $REF && samtools faidx $RF && Index_Ref_SAM_Flag=1 || Index_Ref_SAM_Flag=2 && printf "\nIndexing of reference sequence with Samtools: Done.\n"
elif [ $Index_Ref_BWA_Flag != 1 ] ; then printf "\nError in BWA-Indexing Stage, exiting the program.\n" && exit
fi
wait

##3.2: Creating Dictionary:
if [ -f $REFCHECKDIC ] ; then printf "\n3.2: Dictionary Already Created, proceeding to the next stage.\n" && Create_Dict_Flag=1
elif [ ! -f $REFCHECKDIC ] && [ $Index_Ref_BWA_Flag == 1 ] && [ $Index_Ref_SAM_Flag == 1 ]; then printf "\n${RED}3.3: Initiate Creation of Dictionary.${NC}\n" && cd $ALN && java -Xmx16G -jar $PICARD_HOME/picard.jar CreateSequenceDictionary R=$REF/$RF O=$REF/$RFE.dict && Create_Dict_Flag=1 || Create_Dict_Flag=2 && printf "\nCreation of Dictionary: Done.\n"
elif [ $Index_Ref_SAM_Flag != 1 ] ; then printf "\nError in SAM-Indexing Stage, exiting the program.\n" && exit
fi
wait

##3.3: Aligning the sequences to the reference:
#!!!Add a checkpoint for if one already has aligned the sequence to skip the alignment again. However, what happends if the resulting product is faulty?
if [ -f $BAM/$SAMPLE/$SAMPLE-ALN-PE.bam ] || [ -f $ALN/$SAMPLE-ALN-PE.sam ] ; then printf "\n3.3: Alignment already completed, skipping to the indexing of .bam file.\n" && Aligning_Flag=1 
elif [ ! -f $BAM/$SAMPLE/$SAMPLE-ALN-PE.bam ] && [ $Create_Dict_Flag == 1 ] ; then printf "\n${RED}3.3: Initiate Aligning (BWA MEM) to the Reference Sequence.${NC}\n" && cd $REF && bwa mem -t 10 $RF $FASTQ/$SAMPLE.sra_1.fastq $FASTQ/$SAMPLE.sra_2.fastq > $ALN/$SAMPLE-ALN-PE.sam && Aligning_Flag=1 || Aligning_Flag=2 && printf "\nAlignment with BWA MEM: Done.\n"
elif [ $Index_Ref_BWA_Flag != 1 ] || [ $Index_Ref_BWA_Flag != 1 ] || [ $Index_Ref_BWA_Flag != 1 ] ; then printf "\nError in Indexing Stage, exiting the program.\n" && exit
fi
wait

##3.4: Conversion into .bam file as well as adding ReadGroups:
#!!!Add a checkpoint for if one already has converted the sequence to skip the conversion again. However, what happends if the resulting product is faulty?
if [ -f $BAM/$SAMPLE/$SAMPLE-ALN-PE.bam ] ; then printf "\n3.4: Converison from .sam into .bam already completed, skipping to the indexing of .bam file.\n" && Conversion_Flag=1 
elif [ ! -f $BAM/$SAMPLE/$SAMPLE-ALN-PE.bam ] && [ $Aligning_Flag == 1 ] ; then printf "\n${RED}3.4: Initiate Conversion from .sam into .bam${NC}\n" && cd $ALN && java -Xmx16G -jar $PICARD_HOME/picard.jar AddOrReplaceReadGroups INPUT=$SAMPLE-ALN-PE.sam OUTPUT=$BAM/$SAMPLE/$SAMPLE-ALN-PE.bam SORT_ORDER=coordinate RGID=$SAMPLE-1 RGLB=$SAMPLE-lib RGPL=ILLUMINA RGPU=$SAMPLE-01 RGSM=$SAMPLE && Conversion_Flag=1 || Conversion_Flag=2 && printf "\nConversion from .sam into .bam file: Done\n" 
elif [ $Aligning_Flag != 1 ] ; then printf "\nError in Aligning Stage, exiting the program.\n" && exit
fi
wait

##3.5: Indexing the .bam file:
if [ -f $BAM/$SAMPLE/$SAMPLE-ALN-PE.bai ] && [ $Conversion_Flag == 1 ] ; then printf "\n3.5: Dictionary Already Created, proceeding to the next stage.\n" && Indexing_Alignment_Flag=1
elif [ ! -f $BAM/$SAMPLE/SAMPLE-ALN-PE.bai ] && [ $Conversion_Flag == 1 ] ; then printf "\n${RED}3.5: Initiate indexing the .bam file.${NC}\n" && cd $BAM/$SAMPLE && java -Xmx16G -jar $PICARD_HOME/picard.jar BuildBamIndex INPUT=$SAMPLE-ALN-PE.bam && Indexing_Alingment_Flag=1 || Indexing_Alingment_Flag=2 && printf "\nIndexing of the -bam file: Done.\n"
elif [ $Conversion_Flag != 1 ] ; then printf "\nError in Conversion Stage, Exiting the program.\n" && exit
fi
wait

##3.6: Removing the intermediary files:
#Make a condition that only removes the file if the file actually exists there.
if [ $Index_Ref_BWA_Flag == 1 ] && [ $Index_Ref_SAM_Flag == 1 ] && [ $Create_Dict_Flag == 1 ] && [ $Aligning_Flag == 1 ] && [ $Conversion_Flag == 1 ] ; then printf "\n3.6: Initiate the Removal of intermediary .sam file.\n" && rm $ALN/$SAMPLE-ALN-PE.sam && printf "\nRemoval of intermediary .sam file: Done.\n" #.sam file. 
elif [ $Index_Ref_BWA_Flag != 1 ] || [ $Index_Ref_SAM_Flag != 1 ] || [ $Create_Dict_Flag != 1 ] || [ $Aligning_Flag != 1 ] || [ $Conversion_Flag != 1 ] ; then printf "\nError in an intermediary stage, exiting the removal.\n" && exit
fi
wait

#4: Post-alignment processing:
##4.1: Local Realignment:
printf "\n${RED}4. Initiate Post-alignment processing.\n${NC}"
###4.1.1: Realignment Target Creator:
#!!!Add a checkpoint for if one already has converted the sequence to skip the conversion again. However, what happends if the resulting product is faulty?
if [ $Conversion_Flag == 1 ] ; then printf  "\n${RED}4.1.1: Initiate Creation of Realigning Targets.\n${NC}" && cd $BAM/$SAMPLE && java -Xmx16G -jar $GATK_HOME/GenomeAnalysisTK.jar -I $SAMPLE-ALN-PE.bam -R $REF/$RF -T RealignerTargetCreator -o $SAMPLE.intervals && RTC_Flag=1 || RTC_Flag=2 && printf "\nCreation of Indel targets: Done.\n"
elif [ $Conversion_Flag != 1 ] ; then printf "\nError in an intermediary stage, exiting the program.\n" && exit
fi
wait

###4.1.2: Indel Realinger:
#!!!Add a checkpoint for if one already has converted the sequence to skip the conversion again. However, what happends if the resulting product is faulty?
if [ $RTC_Flag == 1 ] ; then printf  "\n${RED}4.1.2: Initiate Local Realignment based in Indels.${NC}\n" && cd $BAM/$SAMPLE && java -Xmx16G -jar $GATK_HOME/GenomeAnalysisTK.jar -T IndelRealigner -I $SAMPLE-ALN-PE.bam -R $REF/$RF -o $SAMPLE-IDR.bam -targetIntervals $SAMPLE.intervals && Indel_Flag=1 || Indel_Flag=2 && printf "\nLocal Realignment: Done.\n"
elif [ $RTC_Flag != 1 ] ; then printf "\nError in the target creation stage, exiting the program.\n" && exit
fi
wait

##4.2: Duplicate Removal:
#!!!Add a checkpoint for if one already has converted the sequence to skip the conversion again. However, what happends if the resulting product is faulty?
if [ $Indel_Flag == 1 ] ; then printf  "\n${RED}4.2: Initate Duplicate Marking.${NC}\n" && cd $BAM/$SAMPLE && java -Xmx16G -jar $PICARD_HOME/picard.jar MarkDuplicates INPUT=$SAMPLE-IDR.bam OUTPUT=$SAMPLE-DR.bam METRICS_FILE=$SAMPLE-DR-Metrics.txt && Mark_Dup_Flag=1 || Mark_Dup_Flag=2 && printf "\nMarking of Duplicates: Done.\n"
elif [ $Indel_Flag != 1 ] ; then printf "\nError in the indel realignment stage, exiting the program.\n" && exit
fi

##4.3: Indexing:
if [ $Mark_Dup_Flag == 1 ] ; then printf "\n${RED}4.3: Initiate indexing the Realigned+Duplicate Marked file.${NC}\n" && cd $BAM/$SAMPLE && java -Xmx16G -jar $PICARD_HOME/picard.jar BuildBamIndex INPUT=$SAMPLE-IDR.bam && Indexing_Alingment2_Flag=1 || Indexing_Alingment2_Flag=2 && printf "\nIndexing of the Realigned+Duplicate Markerd file: Done.\n"
elif [ $Mark_Dup_Flag != 1 ] ; then printf "\nError in Marking of the duplicate Stage, exiting the program.\n" && exit
fi
wait

##Base Quality Score Recalibration
#Add the BQSR program? Also, one needs then known sites of SNPs? Or something else.
#Finished
printf "\nThere will be several files created from the script:\n1. Quality Control files using FastQC under $QC.\n2. Indexed reference Sequence under $REF.\n3. .bam file of the aligned sample which will be named '$SAMPLE-ALN-PE.bam' under $BAM/$SAMPLE.\n4. .bam file of the aligned sample which have been processed with Indel-Realigner, named '$SAMPLE-IDR.bam' under $BAM/$SAMPLE.\n5. .bam file of the aligned sample which have been processed with MarkDuplicates, named '$SAMPLE-DR.bam' under $BAM/$SAMPLE.\n6. The indexed version of the aligned sample, named '$SAMPLE-MD.bai' under $BAM/$SAMPLE.\n\nThe Program has finished, exiting the program.\n\nFollowing data is for diagnostics: \n"
printf "FASTQC_Flag: $FASTQC_Flag.\nIndex_Ref_BWA_Flag: $Index_Ref_BWA_Flag.\nIndex_Ref_SAM_Flag: $Index_Ref_SAM_Flag.\nCreate_Dict_Flag: $Create_Dict_Flag.\nAligning_Flag: $Aligning_Flag.\n Conversion_Flag: $Conversion_Flag.\nIndexing_Alignment_Flag: $Indexing_Alignment_Flag.\nRTC_Flag: $RTC_Flag.\nIndel_Flag: $Indel_Flag.\nMark_Dup_Flag: $Mark_Dup_Flag.\nIndexing_Alignment2_Flag: $Indexing_Alignment2_Flag.\n"
