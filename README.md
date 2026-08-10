# RNAseq_sleep_developmentShank3
This repository contains code for the analysis of RNA-seq and EEG spectral data in male wild-type (WT) and Shank3 mutant (S3) mice at three ages: postnatal day 24 (P24), P30 and adult. The aim of this study is to assess the molecular and physiological response to sleep deprivation and recovery sleep in WT versus S3 mice across development and into adulthood. 
RNA-seq data is for 3 hours of sleep deprivation (SD) in P24/P30 mice, 5 hours of SD in adult mice and 5 hours of SD followed by 2 hours of recovery sleep (SD+RS) in adult mice. Each time point was compared to allowed to sleep controls which were left undisturbed in their home cage for the duration of the experiment (HC3, HC5, HC7). 
WT data was obtained from GSE211301 (Muheim et al., 2023, SD and HC, P24, P30) and GSE113754 (Ingiosi et al., 2019, SD and HC, adult). All S3 data and WT adult SD+RS and HC data are available at GSE306621. 
EEG data is for a baseline recording day and a 3 hour SD + recovery day in WT and S3, P24, P30 and is available at https://sleepdata.org/datasets/medina-2022 (Medina et al., 2022). 

## Authors:
Elliot Wald (elliot.wald@wsu.edu)

Elizabeth Medina (325medinaelizabeth@gmail.com)

Catilin Ottaway (caitlin.ottaway@wsu.edu)

Alex Popescu (alex.popescu@yale.edu)

Katie Ford (kaitlyn.ford@wsu.edu)

Christine Muheim (christine.muheim@wsu.edu)

Stephanie Hicks (shicks19@jhu.edu)

Michael Rempe (michael.rempe@wsu.edu)

## Data
### RNA-seq
All sequencing data have been deposited in the Gene Expression Omnibus Database (GEO) under accession numbers GSE30661 (S3, SD, SD+RS, HC, all ages; WT, SD+RS and HC, adult), GSE211301 (WT, SD and HC, P24, P30) and GSE113754 (WT, SD and HC, adult). 

We already had gzipped FASTQ files for S3 samples and adult WT samples. WT P24/P30 FASTQ files were downloaded from SRA using the SRA toolkit described here: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

First, we identified which samples needed to be downloaded following `01_quantification/download-geo-data.R`

Next, we downloaded the appropriate samples using prefetch and fasterq-dump:
 
```
d=/Users/elliotwald/Desktop
 
prefetch --option-file $d/SRR_files_WTP24_WTP30_060225.txt
```

``` 
a=/Users/elliotwald/Desktop
 
for file in $(<$a/SRR_files_WTP24_WTP30_060225.txt)
do
    echo "$file.sra"
    fasterq-dump -O $a/WT_P24_P30_FASTQ/ \
        -f --threads 6 \
        $a/WT_P24_P30/sra/$file.sra
done
```
 
Once files were downloaded, they were gzipped for use in the Salmon quantification pipeline. 

### EEG 
Raw EEG data is available at https://sleepdata.org/datasets/medina-2022 (Medina et al., 2022) and can be downloaded using NSRR gem as described here: https://github.com/nsrr/nsrr-gem/blob/master/README.md#prerequisites 

## RNA-seq Analysis
### Transcript Quantification
Raw sequencing reads were quantified using Salmon v1.10.0. The most recent Salmon release (July 2026) is v2.4.1. Salmon 2.0 was re-written in Rust as described here: https://github.com/COMBINE-lab/salmon.
The last C++ release (v1.12.0) is described here: https://salmon.readthedocs.io/en/latest/index.html. 

The `01_quantification/1.1_create-decoys-salmon.sh` and `01_quantification/1.2_build-index-salmon.sh` scripts contain code for building a decoy-aware index using GENCODE release M37 and genome assembly version GRCm39.
Files for the most recent GENCODE release can be found at https://www.gencodegenes.org/mouse/. 

The `01_quantification/1.3_run-salmon-quant.sh` script contains code for running Salmon’s mapping-based mode to quantify FASTQ files.
We specified `-l A` to allow Salmon to automatically infer the library type and `–numBootstraps 30` to generate bootstrap samples for Fishpond Swish. 

The package tximeta (v1.24.0) from R/Bioconductor (Love et al., 2020) was used to import quantified outputs from Salmon into RStudio following `01_quantification/1.4_run-tximeta-for-fishpond.R`. 
Included in this script is code to make a linkedTxome which saves metadata information about the transcriptome used for quantification (decoy-aware index file, FASTA file of transcript sequences and GTF annotation file).
The resulting JSON file is `01_quantification/gencode.vM37-salmon-index-v1.0.0-mouse-withdecoys.json`. Two SummarizedExperiment (SE) objects were created – one for transcript-level quantifications and one for gene-level quantifications.
Count information was extracted from the objects and can be found in `01_quantification` as `se_HCSD_WT_Shank3_SleepDevelopment_salmon_063025.txt.gz` and `gse_HCSD_WT_Shank3_SleepDevelopment_salmon_063025.txt` for transcript-level and gene-level respectively.
These files are also available through GSE30661.

### Normalization and Differential Expression Analysis
The `02_differential_expression_analysis/2.2_differential_gene_expression.R` script contains code for RUV normalization (Risso et al., 2014) and Fishpond Swish differential expression analysis (Zhu et al., 2019).
The quant files from Salmon were imported as a SE object and summarized to the gene level, keeping inferential replicates for the Swish method. Genes with less than 10 reads across 5 samples were filtered out of the dataset to avoid zero inflation.
The SE object was then split into juvenile samples (P24/P30) and adult samples, as statistical analysis was performed for each group independently. Filtered (but not normalized) gene lists were extracted from each group to serve as gene background
from the mouse cortex during functional enrichment analysis as described in the following section. These lists can be found in `03_functional_enrichment_analysis` as `adultBackground.txt` and `juvBackground.txt`.
Upper-quartile normalization was completed using `betweenLaneNormalization` from the EDASeq package v2.40.0 to account for differences in library size, then followed by RUVseq normalization using `RUVs` from the RUVSeq
package v1.40.0 to account for biological noise. For RUVs, we used technical replicates to estimate unwanted factors. Parameter k=15 was selected for both adult and juvenile datasets based on positive/negative control recovery in adult WT samples,
sample separation across PC1-3 on principal component analysis plots and differentially expressed gene (DEG) counts. The positive and negative control genes used with adult WT samples were obtained from Popescu & Ottaway et al., 2025 and
mapped to Ensembl release 114 using `02_differential_expression_analysis/2.1_control_gene_annotation.R`.  

Differential gene expression was estimated using the nonparametric expression analysis method from Fishpond Swish (fishpond package v. 2.12.0) with an FDR < 0.05. Comparisons were as follows for each genotype:
HCP24 vs. SDP24, HCP30 vs. SDP30, HCP90 vs. SDP90 and HCP90 vs. SD+RSP90 (a total of 8 comparisons). Lists of all expressed genes and total DEGs were extracted and saved for downstream applications. 

For each comparison, differential gene expression was visualized by plotting the log2 fold change of all expressed genes in WT versus S3 mice and color-coding based on FDR cutoff and group specificity.
Code for making these plots can be found in `02_differential_expression_analysis/2.3_logLog_plot_visualization.R`. The input for this script is the lists of all expressed genes which are available 
at `02_differential_expression_analysis/allExpressed_100125_K15.xlsx`. 

To further prioritize genes that are unique to each group, lists of total DEGs were intersected and visualized using the UpSetR v1.14.0 and ComplexUpset v1.3.3 packages from R/Bioconductor.
Code for making UpSet plots can be found in `02_differential_expression_analysis/2.4_UpSet_plot_visualization.R`. The input for this script is the lists of total DEGs which are available at
`02_differential_expression_analysis/total_DEG_lists_100125_K15.xlsx`. UpSet plot analysis allowed us to define all possible intersections between timepoints and genotypes for upregulated and downregulated DEGs. 
Intersections representing unique genes (genes only associated with one genotype and timepoint) were selected for further analysis. 

Lists of unique DEGs as defined by UpSet plot intersections can be found at `02_differential_expression_analysis/upsetPlot_geneLists_orderedIntersections_0723.xlsx`. These lists were 
assembled as supplements to the paper using `02_differential_expression_analysis/2.5_assemble_geneLists.R` which takes `total_DEG_lists_100125_K15.xlsx` and `upsetPlot_geneLists_orderedIntersections_0723.xlsx` as input.

### Functional Enrichment Analysis
Functional annotation of the unique gene lists was performed using the Database for Annotation, Visualization, and Integrated Discovery v2021 (DAVID) with DAVID Knowledgebase v2025_1 (Sherman et al., 2022). 
The UniProt terms biological process (BP) and molecular function (MF) along with the Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways were selected as the annotation categories. For each functional term, enrichment was defined relative to all expressed genes in the frontal cortex (21,852 after filtering)
using EASE score <0.05. All enriched terms were clustered based on overlap of genes using final group memberships of three and a similarity threshold of 0.20 (at least 20% overlap in genes). DAVID output was saved as `03_quantification/functionalAnnotation_P24-30-90_intersections.xlsx` 
and serves as input for both `03_functional_enrichment_analysis/3.1_functional_enrichment_plot_visualization.R` and `03_functional_enrichment_analysis/3.2_annotate_DAVID_output.R`. 

`03_functional_enrichment_analysis/3.1_functional_enrichment_plot_visualization` contains code for visualizing both clustered and unclustered enriched terms using ggplot2 (v. 3.5.2).

`03_functional_enrichment_analysis/3.2_annotate_DAVID_output.R` contains code for annotating the output from DAVID with gene names, identifying hub genes, and finding recovered positive controls in adult WT samples.
Hub genes were defined by identifying genes that appeared in all or almost all of the individual terms for a functional annotation cluster.

The annotated DAVID output is available at `03_functional_enrichment_analysis/annotatedDAVID_output_uniqueIntersections_0201.xlsx` and serves as input for the `03_functional_enrichment_analysis/3.3_pathway_comparison_across_conditions.R` script. 
This script contains code for visualizing terms that are functionally enriched in at least three experimental groups, an analysis that allowed for the identification of functions that were common across ages and genotypes which may mediate sleep homeostasis deficits in the S3 animals.

## Spectral Analysis
`04_EEG_analysis` contains code for time in state analysis of P24/P30 EEG which is modified from https://github.com/PeixotoLab/MeCP2Z (Al Maghribi et al., 2025). This code was initially written to use sex (male or female) and genotype (WT or S3) as the input variables. In `InputParams_p24_p30_for_spectra.m`, the sex variable was changed to be `S.Sexes    	   = {'P24','P30'}`. For the remaining scripts, when sex is listed as the name of the variable, the actual variable being used is age as defined in  `InputParams_p24_p30_for_spectra.m`. 

An additional note is that this code was written to run all spectral analyses at once rather than specific tests. It takes .mat files as input and produces figures, statistics and tables for time in state, spectral power, bout analysis, etc. Only the time in state analysis was reported in this paper. For more information see https://github.com/PeixotoLab/MeCP2Z and the corresponding paper. 

The requirements for running this code are Matlab, the control system toolbox, the statistics/machine learning toolbox and the bioinformatics toolbox. The starting point for each recording is a .mat file created with `ConvertSleepSignCSVtoMatlabMAT.m`. 
