#!/bin/bash

## Define the key directories

PROJROOT=~/KS_transcriptome
RAWFQ=${PROJROOT}/01_raw_data/fasta
RAWQC=${PROJROOT}/01_raw_data/fastqc

## Global parameters
THREADS=10

##----------------------------------##
    	##Checking Directories##
##----------------------------------##

## Check the project root exists
if [[ -d ${PROJROOT} ]]; then
	echo -e "Found ${PROJROOT}\n"
     else
	echo -e "${PROJROOT} not found. \nExiting with error code 1" >&2
	exit 1
fi

## Check raw data directory exists
if [[ -d ${RAWFQ} ]]; then
        echo -e "Found ${RAWFQ}\n"
     else
        echo -e "${RAWFQ} not found. \nExiting with error code 2" >&2
        exit 2
fi

## Checking if fastqc directory exists if not create it
if [[ -d ${RAWQC} ]]; then
	echo -e "Found ${RAWQC}\n"
else
  	echo -e "${RAWQC} not found. \nCreating directory ${RAWQC}"
	mkdir ${RAWQC}
fi

##----------------------------------##
	       ##Running FastQC##
##----------------------------------##

# Check for existing FastQC files and only run FastQC if required

RAWFQC=$(find "${RAWQC}" -name *zip | head -n1) #Define a variable called RAWFQC and assign it a command. Find will see if any zip files exist in the RAWQC directory.
if [ -z "${RAWFQC}" ]
  then #If the RAWFQC variable is empty (no files were in the RAWFQC directory) then FastQC will be run.
    echo -e "Running FastQC" #Print to Stdout that FastQC is running.
    fastqc -t ${THREADS} -o ${RAWQC} ${RAWFQ}/*gz #Run FastQC for all .gz files in the RAWFQ directory.
  else
    echo -e "FastQC output in ${RAWQC} detected from a previous run and will not be rerun.\n" #If the files already exist then say so and do not run FastQC.
fi

echo -e "Done"

