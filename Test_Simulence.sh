#!/bin/bash

#Set variable names
species=$1
dir_name="${species}_test"
my_wd=$(pwd)

if [ "$species" = "human" ] || [ "$species" = "mouse" ] || [ "$species" = "rat" ]; then 
	: # Do nothing
else
	echo "
	Command line argument 1 should be either 'human', 'mouse', or 'rat'. 
	You entered: $1
	"
	exit 1
fi

#Create directories
mkdir -p ${species}
mkdir -p ${species}/starDB
mkdir -p ${species}/star_results

#Run script to get spike sequences. 10 represents 10%, the % of sequences to use as spike-ins 
python bin/FASTA_Randomly-Select-Sequences.py \
	All-Sequences/${species}_tRNAs-and-ncRNAs_relative_cdhit.fa \
	10 \
	${species}/${species}_ncRNAs 

#Concatenate the output of this script: 
cat ${species}/${species}_ncRNAs_* > ${species}/All-${species}-ncRNAs_spiked.fa 

#Simulate reads from the new spiked sequence dataset: 
python bin/Simulence.py \
	FASTQ \
	${species}/All-${species}-ncRNAs_spiked.fa \
	35 \
	100 \
	${species}/All-${species}-ncRNAs_spiked_Simulated-Reads.fq & 

#Build database: 
STAR \
	--runThreadN 5 \
	--runMode genomeGenerate \
	--genomeDir ${my_wd}/${species}/starDB/ \
	--genomeFastaFiles ${my_wd}/${species}/${species}_ncRNAs_not-spike.fa 

wait

#Run STAR: 
STAR \
	--genomeDir ${my_wd}/${species}/starDB/ \
	--runThreadN 10 \
	--readFilesIn ${my_wd}/${species}/All-${species}-ncRNAs_spiked_Simulated-Reads.fq \
	--outFileNamePrefix ${species}/star_results/ 

#Collapse aligned reads and remove multi-mappers: 

bin/./SAMcollapse.sh ${species}/star_results/Aligned.out.sam 20 

mv Collapsed.sam ${species}/

#Run SAM-to-Stats.py: 
python bin/SAM-to-Stats.py ${species}/Collapsed.sam ${species}/All-${species}-ncRNAs_spiked_Simulated-Reads.fq.log ${species}/All-${species}-ncRNAs_spiked_Simulated-Reads_Stats.txt

# Remove crud
rm Log.out
