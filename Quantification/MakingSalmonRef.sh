#!/bin/sh

#SBATCH --partition=peixoto

#SBATCH --job-name=SalmonRef.sh

#SBATCH --error=SalmonRef.err

#SBATCH --output=SalmonRef.out

#SBATCH --cpus-per-task=24

#SBATCH --ntasks=1


# Go to the location you will be creating your index at, for me it was
#### /data/peixoto/WTShank3_SleepDevelopment/SalmonQuant3/Gencode_vM32/SalmonIndexFiles




wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/GRCm39.primary_assembly.genome.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.transcripts.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz



###############################################
########## Files here are for bulk ############
# to index the reference features
###############################################


# 1. Salmon indexing requires the names of the genome targets, which is extractable by using the grep command:

grep "^>" <(gunzip -c GRCm39.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys.txt

sed -i.bak -e 's/>//g' decoys.txt

# 2. Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference

cat gencode.vM32.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz > gentrome_transcripts_mouse.fa.gz

######### At the end here, you should get a total of 6 files.. 1. decoys_mouse.txt 2. gencode.vM32.annotation.gtf.gz 3.gentrome_transcripts_mouse.fa.gz 4.decoys_mouse.txt.bak 5.gencode.vM32.transcripts.fa.gz 6.GRCm39.primary_assembly.genome.fa.gz

