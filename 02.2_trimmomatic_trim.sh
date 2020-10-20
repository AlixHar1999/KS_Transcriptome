#!/bin/bash

# Define the key directories

PROJROOT=~/KS_transcriptome
RAWFQ=${PROJROOT}/01_raw_data/fasta
TRIMMED=${PROJROOT}/02_trimmed_data/trimmomatic
CONDA_PREFIX=~/anaconda3

# Global parameters
THREADS=10

##----------------------------------##
	     ##Checking Directories##
##----------------------------------##

# Check the project root exists
if [[ -d ${PROJROOT} ]]; then
        echo -e "Found ${PROJROOT}\n"
     else
        echo -e "${PROJROOT} not found. \nExiting with error code 1" >&2
        exit 1
fi

# Checking if 02_trimmed_data directory exists if not create it
if [[ -d ${TRIMMED} ]]; then
        echo -e "Found ${TRIMMED}\n"
else
        echo -e "${TRIMMED} not found. \nCreating directory ${TRIMMED}\n"
        mkdir -p ${TRIMMED}
fi

# Checking if trimmed_data directory exists if not create it
if [[ -d ${TRIMMED}/trimmed_data ]]; then
        echo -e "Found ${TRIMMED}/trimmed_data\n"
else
        echo -e "${TRIMMED}/trimmed_data not found. \nCreating directory ${TRIMMED}/trimmed_data\n"
        mkdir -p ${TRIMMED}/trimmed_data
fi

# Checking if fastqc directory exists if not create it
if [[ -d ${TRIMMED}/fastqc ]]; then
        echo -e "Found ${TRIMMED}/fastqc\n"
else
        echo -e "${TRIMMED}/fastqc not found. \nCreating directory ${TRIMMED}/fastqc\n"
        mkdir -p ${TRIMMED}/fastqc
fi

# Checking if discarded directory exists if not create it
if [[ -d ${TRIMMED}/discarded ]]; then
        echo -e "Found ${TRIMMED}/discarded\n"
else
        echo -e "${TRIMMED}/discarded not found. \nCreating directory ${TRIMMED}/discarded\n"
        mkdir -p ${TRIMMED}/discarded
fi

# Checking if log directory exists if not create it
if [[ -d ${TRIMMED}/log ]]; then
        echo -e "Found ${TRIMMED}/log\n"
else
        echo -e "${TRIMMED}/log not found. \nCreating directory ${TRIMMED}/log\n"
        mkdir -p ${TRIMMED}/log
fi

##----------------------------------##
	        ##Running Fastp##
##----------------------------------##

# Run trimmomatic on the raw data to trim the data if it has not already been done.

TRIMMOMATIC=$(find "${TRIMMED}/trimmed_data" -name *.gz | head -n1) #Define a variable called TRIMMOMATIC and assign it a command. #Find will see if any .gz files exist in the trimmed_data directory.
if [ -z "${TRIMMOMATIC}" ] #If the TRIMMOMATIC variable is empty (no files were in the trimmed_data directory) then the for loop will be run for the leaf, root and stem data.
    then 
        for R1 in ${RAWFQ}/*R1.fq.gz #For each of read 1 in RAWFQ do the following.
          do
            # Find read 1 file
            R1BASE=$(basename ${R1}) #Extract the basename of the R1 file pathway and hold it in the R1BASE variable.  
            echo -e "Found ${R1BASE}" #Print to stdout that the R1BASE file was found.
            
            # Define read 2 file
            R2=${R1%_R1.fq.gz}_R2.fq.gz #Define read 2 as the R1 variable by chopping off _R1.fq.gz using % from the end of the file name and replace with _R2.fq.gz
            R2BASE=$(basename ${R2}) #Extract the basename of the R2 file pathway and hold it in the R2BASE variable.
            echo -e "Found ${R2BASE}\n" #Print to stdout that the R2BASE file was found.
            
            # Run trimmomatic
            echo -e "Running trimmomatic for ${R1BASE} and ${R2BASE}\n" #Print to stdout that TRIMMOMATIC is being run for the paired end data files.
            trimmomatic PE \
              -threads ${THREADS} \
              -summary ${TRIMMED}/log/${R1BASE%_R1.fq.gz}_trimmomatic.log \
              ${R1} ${R2} \
              ${TRIMMED}/trimmed_data/${R1BASE%_R1.fq.gz}_trimmomatic_trim_R1.fq.gz \
              ${TRIMMED}/discarded/${R1BASE%_R1.fq.gz}_1U.fq.gz \
              ${TRIMMED}/trimmed_data/${R1BASE%_R1.fq.gz}_trimmomatic_trim_R2.fq.gz \
              ${TRIMMED}/discarded/${R1BASE%_R1.fq.gz}_2U.fq.gz \
              ILLUMINACLIP:${CONDA_PREFIX}/share/trimmomatic-0.39-1/adapters/TruSeq3-PE.fa:2:30:10:3:true \
              SLIDINGWINDOW:4:20 \
              LEADING:5 \
              TRAILING:5 \
              MINLEN:30 #Trimmomatic will be run for the paired end data files, output to trimmed_data, failed reads to discarded and a log file will be made.
        done #This will run three time first for the leaf data then for the stem and root data.
  else
   echo -e "Fastp output in ${TRIMMED}/trimmed_data detected from a previous run and will not be rerun.\n" #If the trimmed files already exist then say so and do not run fastp.
fi 

##----------------------------------##
          ##Running FastQC##
##----------------------------------##

#Check for existing FastQC files and only run FastQC if required.

TRIMFQC=$(find "${TRIMMED}/fastqc" -name *zip | head -n1) #Define a variable called TRIMFQC and assign it a command. Find will see if any zip files exist in the ${TRIMMED}/fastqc directory.
if [ -z "${TRIMFQC}" ] #If the TRIMFQC variable is empty (no files were in the ${TRIMMED}/fastqc directory) then FastQC will be run.
  then 
    echo -e "Running FastQC" #Print to Stdout that FastQC is running.
    fastqc -t ${THREADS} -o ${TRIMMED}/fastqc ${TRIMMED}/trimmed_data/*gz #Run FastQC for all .gz files in the ${TRIMMED}/fastqc directory.
  else
    echo -e "FastQC output in ${TRIMMED}/fastqc detected from a previous run and will not be rerun.\n" #If the files already exist then say so and do not run FastQC.
fi

echo -e "Done"
