#!/bin/bash

# Define the key directories

PROJROOT=~/KS_transcriptome
BUSCOOUT=${PROJROOT}/04_assembly_assessment/busco

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

# Check the BUSCO output directory exists
if [[ -d ${BUSCOOUT} ]]; then
	echo -e "Found ${BUSCOOUT}\n"
     else
	echo -e "${BUSCOOUT} not found. \nCreating directory ${BUSCOOUT}\n"
	mkdir -p ${BUSCOOUT}
fi

##---------------------------##
          ##BUSCO##
##---------------------------##

cd ${BUSCOOUT}

# Running BUSCO for leaf assembly

LEAFFASTA=${PROJROOT}/03_assembly/trinity_tissue_specific/KS_leaf30_trinity/Trinity.fasta

echo -e "Running BUSCO for leaf assembly"

busco \
  -c ${CPU} \
  -i ${LEAFFASTA} \
  -l fabales_odb10 \
  -o fabales_leaf30 \
  -m transcriptome
  
# Running BUSCO for root assembly
  
ROOTFASTA=${PROJROOT}/03_assembly/trinity_tissue_specific/KS_root30_trinity/Trinity.fasta

echo -e "Running BUSCO for root assembly"

busco \
  -c ${CPU} \
  -i ${ROOTFASTA} \
  -l fabales_odb10 \
  -o fabales_root30 \
  -m transcriptome
  
# Running BUSCO for stem assembly
  
STEMFASTA=${PROJROOT}/03_assembly/trinity_tissue_specific/KS_stem30_trinity/Trinity.fasta

echo -e "Running BUSCO for root assembly"

busco \
  -c ${CPU} \
  -i ${STEMFASTA} \
  -l fabales_odb10 \
  -o fabales_stem30 \
  -m transcriptome
  
# Running BUSCO for leaf, root and stem combined assembly with min length 30
  
LRSFASTA=${PROJROOT}/03_assembly/trinity_all_tissues/KS_LRS30_trinity/Trinity.fasta

echo -e "Running BUSCO for root assembly"

busco \
  -c ${CPU} \
  -i ${LRSFASTA} \
  -l fabales_odb10 \
  -o fabales_LRS30 \
  -m transcriptome
  
# Running BUSCO for leaf, root and stem combined assembly with min length 100
  
LRS100FASTA=${PROJROOT}/03_assembly/trinity_all_tissues/KS_LRS100_trinity/Trinity.fasta

echo -e "Running BUSCO for root assembly"

busco \
  -c ${CPU} \
  -i ${LRS100FASTA} \
  -l fabales_odb10 \
  -o fabales_LRS100 \
  -m transcriptome