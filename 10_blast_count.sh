#!/bin/bash

# Define the key directories and files
PROJROOT=~/KS_transcriptome
BLASTOUT=${PROJROOT}/04_assembly_assessment/blast/blastout
FILEOUT1=${PROJROOT}/04_assembly_assessment/blast/Species_count
FILEOUT2=${PROJROOT}/04_assembly_assessment/blast/Protein_count
FASTAINPUT=${PROJROOT}/03_assembly/All_Fasta_Files

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

if [[ -d ${BLASTOUT} ]]; then
	echo -e "Found ${BLASTOUT}\n"
     else
	echo -e "${BLASTOUT} not found. \nExiting with error code 2" >&2
	exit 2
fi

if [[ -d ${FILEOUT1} ]]; then
	echo -e "Found ${FILEOUT1}\n"
     else
	echo -e "${FILEOUT1} not found. \nCreating directory ${FILEOUT1}\n"
	mkdir -p ${FILEOUT1}
fi 

if [[ -d ${FILEOUT2} ]]; then
	echo -e "Found ${FILEOUT2}\n"
     else
	echo -e "${FILEOUT2} not found. \nCreating directory ${FILEOUT2}\n"
	mkdir -p ${FILEOUT2}
fi

#Create txt files that count the total number of unique blast hits from each transcript for each species.
for TSV in ${BLASTOUT}/*.tsv ; 
  do 

  BASE=$(basename ${TSV})
  echo -e "Found ${BASE}."

  echo -e "#" > ${FILEOUT1}/${BASE%.tsv}_Species_count.txt
  
  sort ${BLASTOUT}/${BASE} | uniq -c | sort -nr >> ${FILEOUT1}/${BASE%.tsv}_Species_count.txt
  
done 

#Create table that shows total blast hits, unique blast hits per transcript and number of unique proteins.

echo -e "Assembly\tTotal_Blast_Hits\tUnique_Blast_Hits\tUnique_Proteins\tTotal_Transcripts" > ${FILEOUT2}/Protein_count.tsv

#Assembly 1

A1="01 Trinity leaf 30"
TBH=$(wc -l < ${BLASTOUT}/KS_leaf30_trinity_blastx.outfmt6)
UBH=$(wc -l < ${BLASTOUT}/KS_leaf30_trinity_blastx.outfmt6.w_pct_hit_length)
UP=$(wc -l < ${BLASTOUT}/KS_leaf30_trinity_blastx.outfmt6.hist.list)
TT=$(grep -c '^>' ${FASTAINPUT}/KS_leaf30_trinity.fa)

echo -e "${A1}\t${TBH}\t${UBH}\t${UP}\t${TT}" >> ${FILEOUT2}/Protein_count.tsv

#Assembly 2

A1="02 Trinity root 30"
TBH=$(wc -l < ${BLASTOUT}/KS_root30_trinity_blastx.outfmt6)
UBH=$(wc -l < ${BLASTOUT}/KS_root30_trinity_blastx.outfmt6.w_pct_hit_length)
UP=$(wc -l < ${BLASTOUT}/KS_root30_trinity_blastx.outfmt6.hist.list)
TT=$(grep -c '^>' ${FASTAINPUT}/KS_root30_trinity.fa)

echo -e "${A1}\t${TBH}\t${UBH}\t${UP}\t${TT}" >> ${FILEOUT2}/Protein_count.tsv

#Assembly 3

A1="03 Trinity stem 30"
TBH=$(wc -l < ${BLASTOUT}/KS_stem30_trinity_blastx.outfmt6)
UBH=$(wc -l < ${BLASTOUT}/KS_stem30_trinity_blastx.outfmt6.w_pct_hit_length)
UP=$(wc -l < ${BLASTOUT}/KS_stem30_trinity_blastx.outfmt6.hist.list)
TT=$(grep -c '^>' ${FASTAINPUT}/KS_stem30_trinity.fa)

echo -e "${A1}\t${TBH}\t${UBH}\t${UP}\t${TT}" >> ${FILEOUT2}/Protein_count.tsv

#Assembly 4

A1="04 Trinity All Tissues 30"
TBH=$(wc -l < ${BLASTOUT}/KS_LRS30_trinity_blastx.outfmt6)
UBH=$(wc -l < ${BLASTOUT}/KS_LRS30_trinity_blastx.outfmt6.w_pct_hit_length)
UP=$(wc -l < ${BLASTOUT}/KS_LRS30_trinity_blastx.outfmt6.hist.list)
TT=$(grep -c '^>' ${FASTAINPUT}/KS_LRS30_trinity.fa)

echo -e "${A1}\t${TBH}\t${UBH}\t${UP}\t${TT}" >> ${FILEOUT2}/Protein_count.tsv

#Assembly 5

A1="05 Trinity All Tissues 100"
TBH=$(wc -l < ${BLASTOUT}/KS_LRS100_trinity_blastx.outfmt6)
UBH=$(wc -l < ${BLASTOUT}/KS_LRS100_trinity_blastx.outfmt6.w_pct_hit_length)
UP=$(wc -l < ${BLASTOUT}/KS_LRS100_trinity_blastx.outfmt6.hist.list)
TT=$(grep -c '^>' ${FASTAINPUT}/KS_LRS100_trinity.fa)

echo -e "${A1}\t${TBH}\t${UBH}\t${UP}\t${TT}" >> ${FILEOUT2}/Protein_count.tsv

#Assembly 6

A1="06 Trinity STAR 100"
TBH=$(wc -l < ${BLASTOUT}/KS_trinity_star_blastx.outfmt6)
UBH=$(wc -l < ${BLASTOUT}/KS_trinity_star_blastx.outfmt6.w_pct_hit_length)
UP=$(wc -l < ${BLASTOUT}/KS_trinity_star_blastx.outfmt6.hist.list)
TT=$(grep -c '^>' ${FASTAINPUT}/KS_trinity_star.fa)

echo -e "${A1}\t${TBH}\t${UBH}\t${UP}\t${TT}" >> ${FILEOUT2}/Protein_count.tsv

#Assembly 7

A1="07 Trinity HISAT2 100"
TBH=$(wc -l < ${BLASTOUT}/KS_trinity_hisat2_blastx.outfmt6)
UBH=$(wc -l < ${BLASTOUT}/KS_trinity_hisat2_blastx.outfmt6.w_pct_hit_length)
UP=$(wc -l < ${BLASTOUT}/KS_trinity_hisat2_blastx.outfmt6.hist.list)
TT=$(grep -c '^>' ${FASTAINPUT}/KS_trinity_hisat2.fa)

echo -e "${A1}\t${TBH}\t${UBH}\t${UP}\t${TT}" >> ${FILEOUT2}/Protein_count.tsv

#Assembly 8

A1="08 Stringtie STAR 100"
TBH=$(wc -l < ${BLASTOUT}/KS_LRS100_stringtie_star_blastx.outfmt6)
UBH=$(wc -l < ${BLASTOUT}/KS_LRS100_stringtie_star_blastx.outfmt6.w_pct_hit_length)
UP=$(wc -l < ${BLASTOUT}/KS_LRS100_stringtie_star_blastx.outfmt6.hist.list)
TT=$(grep -c '^>' ${FASTAINPUT}/KS_LRS100_stringtie_star.fa)

echo -e "${A1}\t${TBH}\t${UBH}\t${UP}\t${TT}" >> ${FILEOUT2}/Protein_count.tsv

#Assembly 9

A1="09 Stringtie HISAT2 100"
TBH=$(wc -l < ${BLASTOUT}/KS_LRS100_stringtie_hisat_blastx.outfmt6)
UBH=$(wc -l < ${BLASTOUT}/KS_LRS100_stringtie_hisat_blastx.outfmt6.w_pct_hit_length)
UP=$(wc -l < ${BLASTOUT}/KS_LRS100_stringtie_hisat_blastx.outfmt6.hist.list)
TT=$(grep -c '^>' ${FASTAINPUT}/KS_LRS100_stringtie_hisat.fa)

echo -e "${A1}\t${TBH}\t${UBH}\t${UP}\t${TT}" >> ${FILEOUT2}/Protein_count.tsv