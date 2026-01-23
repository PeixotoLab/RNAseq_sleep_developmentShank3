#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=2G,h_vmem=2G

###############################################
########## Files here are for bulk ############
# to index the reference features
###############################################

# 1. Salmon indexing requires the names of the genome targets, which is extractable by using the grep command:
grep "^>" <(gunzip -c GRCm39.primary_assembly.genome.fa.gz) | cut -d " " -f 1 > decoys_mouse.txt
sed -i.bak -e 's/>//g' decoys_mouse.txt

# 2. Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
cat gencode.vM37.transcripts.fa.gz GRCm39.primary_assembly.genome.fa.gz > gentrome_transcripts_mouse.fa.gz
