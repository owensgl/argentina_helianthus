# argentina_helianthus
Scripts used for "Gene flow in Argentinian sunflowers as revealed by genotyping by sequencing data"

## Aligning and calling SNPs
* GBS_fastq_Demultiplexer_v9_2Enzyme2barcode.pl: Demultiplexing the raw fastq files into individual samples and trimming adapter sequence.
* align_process.bash: Initial steps to trim raw data, align using nextgenmap, and format BAM files.
* make_gcvf_gatk.bash: Realigning around indels and creating gvcf files.
* Genotype_gvcf.sh: Genotyping gvcf files to create vcf.
* vcf2vertical_dep_GATK36_full.pl: Converts a vcf file to a flat tab separated genotype file.
## Analyses
* hiest.R: HIest analysis and plotting.
* PCA.R: PCA analysis with Argentina samples.
* PCA_with_other_species.R: PCA analysis with Argentina and USA samples.
* plot_ngsadmix.R: Plotting NGSadmix results.
* multiplot.R: Utility for combining plots.
* Ana_sample_info_prelim.txt: Sample information.
* PCA_ngsadmix_results.tab: Combined results table. 
