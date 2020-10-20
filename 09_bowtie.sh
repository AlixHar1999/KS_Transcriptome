#!/bin/bash

# Define the key directories and files

PROJROOT=~/KS_transcriptome
readDIR=${PROJROOT}/02_trimmed_data/fastp_100/trimmed_data
FASTAINPUT=${PROJROOT}/03_assembly/All_Fasta_Files
BOWTIEINDEX=${PROJROOT}/04_assembly_assessment/bowtie/index
BOWTIEOUT=${PROJROOT}/04_assembly_assessment/bowtie/alignments

# Global parameters
CPU=20

##---------------------------##
    ##Checking Directories##
##---------------------------##

# Check the project root exists
if [[ -d ${PROJROOT} ]]; then
	echo -e "Found ${PROJROOT}\n"
     else
	echo -e "${PROJROOT} not found. \nExiting with error code 1" >&2
	exit 1
fi

# Check the fasta files exists
if [[ -d ${FASTAINPUT} ]]; then
	echo -e "Found ${FASTAINPUT}\n"
     else
	echo -e "${FASTAINPUT} not found. \nExiting with error code 2" >&2
	exit 2
fi

# Checking if bowtie index directory exists if not create it
if [[ -d ${BOWTIEINDEX} ]]; then
        echo -e "Found ${BOWTIEINDEX}\n"
else
        echo -e "${BOWTIEINDEX} not found. \nCreating directory ${BOWTIEINDEX}\n"
        mkdir -p ${BOWTIEINDEX}
fi

# Checking if bowtie output directory exists if not create it
if [[ -d ${BOWTIEOUT} ]]; then
        echo -e "Found ${BOWTIEOUT}\n"
else
        echo -e "${BOWTIEOUT} not found. \nCreating directory ${BOWTIEOUT}\n"
        mkdir -p ${BOWTIEOUT}
fi

##---------------------------##
      ##Running bowtie##
##---------------------------##

R1=${readDIR}/KS_LRS_fastp100_R1.fq.gz
R2=${readDIR}/KS_LRS_fastp100_R2.fq.gz

for FASTA in ${FASTAINPUT}/*fa ; 
  do 

  FASTABASE=$(basename ${FASTA})
  echo -e "Found ${FASTABASE}. Creating bowtie index then running bowtie."

  bowtie2-build --threads ${CPU} ${FASTAINPUT}/${FASTABASE} ${BOWTIEINDEX}/${FASTABASE}
  
  bowtie2 -p ${CPU} -q --no-unal -k 20 -x ${BOWTIEINDEX}/${FASTABASE} -1 ${R1} -2 ${R2}  \
     2>${BOWTIEOUT}/${FASTABASE%.fa}_align_stats.txt | samtools view -@10 -Sb -o ${BOWTIEOUT}/${FASTABASE%.fa}_bowtie2.bam 
  
done 