#!/bin/bash

# Define the key directories and files

PROJROOT=~/KS_transcriptome
FASTAINPUT=${PROJROOT}/03_assembly/All_Fasta_Files
BLAST1=${PROJROOT}/04_assembly_assessment/blast/swissprot/uniprot_sprot.fasta
BLASTOUT=${PROJROOT}/04_assembly_assessment/blast/blastout
TRINITY_HOME=~/anaconda3/opt/trinity-2.1.1

# Global parameters
CPU=22

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

# Check the uniprot_sprot.fasta files exists
if [[ -f ${BLAST1} ]]; then
	echo -e "Found ${BLAST1}\n"
     else
	echo -e "${BLAST1} not found. \nExiting with error code 3" >&2
	exit 3
fi

# Checking if blast output directory exists if not create it
if [[ -d ${BLASTOUT} ]]; then
        echo -e "Found ${BLASTOUT}\n"
else
        echo -e "${BLASTOUT} not found. \nCreating directory ${BLASTOUT}\n"
        mkdir -p ${BLASTOUT}
fi

##---------------------------##
      ##BLAST Database##
##---------------------------##

echo -e " Making blast database from swiss prot genes"

makeblastdb -in ${BLAST1} -dbtype prot

##---------------------------##
      ##Running BLAST##
##---------------------------##

cd ${BLASTOUT}

for FASTA in ${FASTAINPUT}/*fa ; 
  do 

  FASTABASE=$(basename ${FASTA})
  echo -e "Found ${FASTABASE}. Running blast search using SwissProt database"

  blastx \
    -query ${FASTAINPUT}/${FASTABASE} \
    -db ${BLAST1} \
    -out ${FASTABASE%.fa}_blastx.outfmt6 \
    -evalue 1e-20 -num_threads ${CPU} \
    -max_target_seqs 1 -outfmt 6
    
  perl ${TRINITY_HOME}/util/analyze_blastPlus_topHit_coverage.pl ${FASTABASE%.fa}_blastx.outfmt6 ${FASTAINPUT}/${FASTABASE} ${BLAST1}
  
done 