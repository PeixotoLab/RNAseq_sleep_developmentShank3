#!/bin/sh

#SBATCH --partition=peixoto
#SBATCH --job-name=SalmonIndex.sh
#SBATCH --error=SalmonIndexerr
#SBATCH --output=SalmonIndex.out
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1

module load salmon

salmon index -t /data/peixoto/WTShank3_SleepDevelopment/SalmonQuant2/SalmonIndexFiles/gentrome_transcripts_mouse.fa.gz \
             -d /data/peixoto/WTShank3_SleepDevelopment/SalmonQuant2/SalmonIndexFiles/decoys.txt \
             -i /data/peixoto/WTShank3_SleepDevelopment/SalmonQuant2/SalmonIndexFiles/gencode.vM32-salmon-index-withdecoys \
             --gencode -p 4 -k 31 
             
             