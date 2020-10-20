#!/bin/bash

# Define the key directories

PROJROOT=~/KS_transcriptome
readDIR=${PROJROOT}/02_trimmed_data/fastp_100/trimmed_data
readOUT1=${PROJROOT}/03_assembly/trinity_tissue_specific
readOUT2=${PROJROOT}/03_assembly/trinity_all_tissues

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
if [[ -d ${readOUT1} ]]; then
	echo -e "Found ${readOUT1}\n"
     else
	echo -e "${readOUT1} not found. \nCreating directory ${readOUT1}\n"
	mkdir -p ${readOUT1}
fi

# Check the triple data Trinity output directory exists
if [[ -d ${readOUT2} ]]; then
	echo -e "Found ${readOUT2}\n"
     else
	echo -e "${readOUT2} not found. \nCreating directory ${readOUT2}\n"
	mkdir -p ${readOUT2}
fi

##-----------------------##
    ##Running Trinity##
##-----------------------##
R1=${readDIR}/KS_LRS_fastp100_R1.fq.gz
R2=${readDIR}/KS_LRS_fastp100_R2.fq.gz

echo -e "Running Trinity using all tissues\n"

Trinity \
  --seqType fq \
  --left ${R1} --right ${R2} \
  --CPU ${CPU} \
  --max_memory ${RAM} \
  --output ${readOUT2}/KS_LRS100_trinity \
  --full_cleanup   

echo -e "Done"