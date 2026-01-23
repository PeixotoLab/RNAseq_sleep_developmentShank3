#!/bin/bash
#SBATCH --partition=free
#SBATCH --job-name=build_salmon_index
#SBATCH --error=build_index_salmon_v1.err
#SBATCH --output=build_index_salmon_v1.out
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
module load salmon

###############################################
########## Files here are for bulk ############
# to build salmon index
# Use
# #$ -l mem_free=10G,h_vmem=10G
# Below is for mouse
###############################################
cd /data/peixoto/development_WT_S3_bulk/salmon_index_files
# create salmon index with decoys (this process takes ~1 hour for mouse)
salmon index -t gentrome_transcripts_mouse.fa.gz \
             -d decoys_mouse.txt \
             -i gencode.vM37-salmon-index-v1.0.0-mouse-withdecoys \
             --gencode --threads 4 -k 31

