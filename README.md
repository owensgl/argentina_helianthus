# argentina_helianthus
Scripts used for Mondon et al., 2017

## Aligning and calling SNPs
* align_process.bash: Initial steps to trim raw data, align using nextgenmap, and format BAM files.
* make_gcvf_gatk.bash: Realigning around indels and creating gvcf files.
* Genotype_gvcf.sh: Genotyping gvcf files to create vcf.
* vcf2vertical_dep_GATK36_full.pl: Converts a vcf file to a flat tab separated genotype file.
## Analyses
* hiest.R: HIest analysis and plotting.
* SNPtable2argentina_diverged_loci.pl: Locates fixed differences between domestic and wild H. annuus and counts the percentage of domestic alleles in each sample.
* domestication_alleles.R: Plotting of domestic allele counts.
* PCA.R: PCA analysis with Argentina samples.
* PCA_with_other_species.R: PCA analysis with Argentina and USA samples.
* plot_ngsadmix.R: Plotting NGSadmix results.
* multiplot.R: Utility for combining plots.
* Ana_sample_info_prelim.txt: Sample information.
* PCA_ngsadmix_results.tab: Combined results table. 
