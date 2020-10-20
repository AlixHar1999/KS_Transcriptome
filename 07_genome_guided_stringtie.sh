#!/bin/bash

# Define the key directories

PROJROOT=~/KS_transcriptome
readOUT=${PROJROOT}/03_assembly/stringtie_all_tissues
BAMSTAR=${PROJROOT}/05_aligned_rnaseq_data/STAR/bam
BAMHISAT=${PROJROOT}/05_aligned_rnaseq_data/HISAT2/bam
GENOME=${PROJROOT}/KS_genome/KS_genome.canu_smartdenovo.fa
BUSCOOUT=${PROJROOT}/04_assembly_assessment/busco

# Global parameters
CPU=20
RAM=250G

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

# Check the parent Trinity output directory exists
if [[ -d ${readOUT} ]]; then
	echo -e "Found ${readOUT}\n"
     else
	echo -e "${readOUT} not found. \nCreating directory ${readOUT}\n"
	mkdir -p ${readOUT}
fi

# Check the STAR bam file directory exists
if [[ -d ${BAMSTAR} ]]; then
	echo -e "Found ${BAMSTAR}\n"
     else
	echo -e "${BAMSTAR} not found. \nExiting with error code 2" >&2
	exit 2
fi

# Check the HISAT2 bam file directory exists
if [[ -d ${BAMHISAT} ]]; then
	echo -e "Found ${BAMHISAT}\n"
     else
	echo -e "${BAMHISAT} not found. \nExiting with error code 3" >&2
	exit 3
fi  
##---------------------------##
    ##Running Stringtie##
##---------------------------##

echo -e "Running Stringtie using STAR alignment\n"

stringtie \
  -o ${readOUT}/KS_LRS100_ref_guided_stringtie_star/KS_LRS100_stringtie_star.gtf \
  -p ${CPU} \
  ${BAMSTAR}/*.bam

gffread \
  -w ${readOUT}/KS_LRS100_ref_guided_stringtie_star/KS_LRS100_stringtie_transcripts.fa \
  -g ${GENOME} \
  ${readOUT}/KS_LRS100_ref_guided_stringtie_star/KS_LRS100_stringtie_star.gtf

echo -e "Running Stringtie using HISAT2 alignment\n"

stringtie \
  -o ${readOUT}/KS_LRS100_ref_guided_stringtie_hisat/KS_LRS100_stringtie_hisat.gtf \
  -p ${CPU} \
  ${BAMHISAT}/*.bam

gffread \
  -w ${readOUT}/KS_LRS100_ref_guided_stringtie_hisat/KS_LRS100_stringtie_transcripts.fa \
  -g ${GENOME} \
  ${readOUT}/KS_LRS100_ref_guided_stringtie_hisat/KS_LRS100_stringtie_hisat.gtf
##---------------------------##
      ##Running BUSCO##
##---------------------------##

cd ${BUSCOOUT}
FASTA1=${readOUT}/KS_LRS100_ref_guided_stringtie_star/KS_LRS100_stringtie_transcripts.fa

echo -e "Running BUSCO for STAR assembly"
busco \
  -c ${CPU} \
  -i ${FASTA1} \
  -l fabales_odb10 \
  -o fabales_LRS100_ref_guided_stringtie_star \
  -m transcriptome

FASTA2=${readOUT}/KS_LRS100_ref_guided_stringtie_hisat/KS_LRS100_stringtie_transcripts.fa

echo -e "Running BUSCO for HISAT2 assembly"
busco \
  -c ${CPU} \
  -i ${FASTA2} \
  -l fabales_odb10 \
  -o fabales_LRS100_ref_guided_stringtie_hisat \
  -m transcriptome