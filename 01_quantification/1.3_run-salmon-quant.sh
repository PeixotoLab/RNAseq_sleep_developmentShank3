#!/bin/bash
#SBATCH --partition=free
#SBATCH --job-name=run-salmon-quant
#SBATCH --error=run-salmon-quant.err
#SBATCH --output=run-salmon-quant.out
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1

module load salmon

a=/data/peixoto/development_WT_S3_bulk

for file in `ls $a/P24_P30_P90_RNAseq_bulk_FASTQ/*.fastq.gz | sed 's/[12].fastq.gz//' | sort -u`;
do 
    fastq_1=`echo $file | sed "s/$/1.fastq.gz/"`
    fastq_2=`echo $file | sed "s/$/2.fastq.gz/"`
    output_file=`basename $file | sed "s/_R//"`
    salmon quant -l A \
        --index $a/salmon_index_files/gencode.vM37-salmon-index-v1.0.0-mouse-withdecoys \
        -1 $fastq_1 \
        -2 $fastq_2 \
        --threads 6 \
        --output $a/01_quantification/bulk/salmon_quants/${output_file}_quant \
        --numBootstraps 30
done


