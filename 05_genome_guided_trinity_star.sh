#!/bin/bash

# Define the key directories

PROJROOT=~/KS_transcriptome
readDIR=${PROJROOT}/02_trimmed_data/fastp_100/trimmed_data
readOUT=${PROJROOT}/03_assembly/trinity_all_tissues
GENOME=${PROJROOT}/KS_genome/KS_genome.canu_smartdenovo.fa
GENOMEINDEX=${PROJROOT}/KS_genome_index/STAR
BAMOUT=${PROJROOT}/05_aligned_rnaseq_data/STAR/bam
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

# Check the trimmed data directory exists
if [[ -d ${readDIR} ]]; then
	echo -e "Found ${readDIR}\n"
     else
	echo -e "${readDIR} not found. \nExiting with error code 2" >&2
	exit 2
fi

# Check the parent Trinity output directory exists
if [[ -d ${readOUT} ]]; then
	echo -e "Found ${readOUT}\n"
     else
	echo -e "${readOUT} not found. \nCreating directory ${readOUT}\n"
	mkdir -p ${readOUT}
fi

# Check the genome index output directory exists
if [[ -d ${GENOMEINDEX} ]]; then
	echo -e "Found ${GENOMEINDEX}\n"
     else
	echo -e "${GENOMEINDEX} not found. \nCreating directory ${GENOMEINDEX}\n"
	mkdir -p ${GENOMEINDEX}
fi

# Check the genome index output directory exists
if [[ -d ${BAMOUT} ]]; then
	echo -e "Found ${BAMOUT}\n"
     else
	echo -e "${BAMOUT} not found. \nCreating directory ${BAMOUT}\n"
	mkdir -p ${BAMOUT}
fi

##---------------------------##
      ##Genome Index##
##---------------------------##

echo -e "Running genome index\n"

STAR \
  --runMode genomeGenerate \
  --runThreadN ${CPU} \
  --genomeDir ${GENOMEINDEX} \
  --genomeFastaFiles ${GENOME}

##---------------------------##
      ##RNA Alignment##
##---------------------------##

echo -e "Running RNAseq alignment\n"
# STAR will be run using default settings apart from the option 
# --outSAMstrandField intronMotif which is needed for use in assembly and 
# --outSAMtype BAM SortedByCoordinate to produce a sorted .bam file without using samtools.
STAR \
  --runThreadN ${CPU} \
  --genomeDir ${GENOMEINDEX} \
  --readFilesIn ${readDIR}/KS_LRS_fastp100_R1.fq.gz ${readDIR}/KS_LRS_fastp100_R2.fq.gz \
  --readFilesCommand gunzip -c \
  --outFileNamePrefix ${BAMOUT}/KS_LRS100. \
  --outSAMstrandField intronMotif \
  --outSAMtype BAM SortedByCoordinate
  
##---------------------------##
      ##Running Trinity##
##---------------------------##

echo -e "Running genome guided Trinity using all tissue samples\n"

Trinity \
  --genome_guided_bam ${BAMOUT}/*.bam \
  --genome_guided_max_intron 10000 \
  --max_memory ${RAM} \
  --CPU ${CPU} \
  --output ${readOUT}/KS_LRS100_ref_guided_trinity_star

echo -e "Done"

##---------------------------##
      ##Running BUSCO##
##---------------------------##

cd ${BUSCOOUT}
FASTA=${PROJROOT}/03_assembly/trinity_all_tissues/KS_LRS100_ref_guided_trinity_star/Trinity-GG.fasta

echo -e "Running BUSCO for root assembly"
busco \
  -c ${CPU} \
  -i ${FASTA} \
  -l fabales_odb10 \
  -o fabales_LRS100_ref_guided_trinity_star \
  -m transcriptome