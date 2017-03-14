library(ggplot2)
library(SNPRelate)
library("bigmemory")
library(dplyr)
library("mapdata")
library("maps")
library("maptools")
library("rworldmap")
library(ggmap)
library(RColorBrewer)

filename <- c("/home/owens/working/Ana/Ana.GATK.2016.plusUSA.80.5.ntab")
gdsname <- c("/home/owens/working/Ana/Ana.GATK.2016.plusUSA.80.5.gds")
#Load genotype matrix
geno_matrix <- read.big.matrix(filename, sep="\t", header = TRUE)
geno_matrix  <- geno_matrix[,-c(1:2)]


genotype_id_df <- read.table(pipe(paste("cut -f1,2", filename,sep=" ")), header = TRUE)
genotype_id_df$id <- paste(genotype_id_df$CHROM,genotype_id_df$POS,sep="_")
genotype_id_df$CHROM <- gsub('Ha', '', genotype_id_df$CHROM) 
genotype_id_df$CHROM <- gsub('_73Ns', '', genotype_id_df$CHROM) 
genotype_id_df$CHROM <- as.numeric(genotype_id_df$CHROM)
sample_id <- colnames(geno_matrix)

snpgdsCreateGeno(gdsname, genmat = geno_matrix, 
                 sample.id = sample_id, snpfirstdim = TRUE, 
                 snp.id = genotype_id_df$id, snp.chromosome = genotype_id_df$CHROM, 
                 snp.position = genotype_id_df$POS)
rm(list=c("geno_matrix"))


##################
#Load gds file
###################

genofile_all <- snpgdsOpen(gdsname, allow.duplicate = TRUE)
snpset_pruned <- snpgdsLDpruning(genofile_all)
snpset.id <- unlist(snpset_pruned)



#############
#Run PCA for Niveus
###############
labels <- read.table("Ana_sample_info_prelim.txt",header=T)
labels$Name <- gsub(".bam", "", labels$Name)
colnames(labels) <- c("sample.id","location","type","plate","category")
labels %>% filter(location != "NA") %>% filter(category != "NA") %>%
  filter(category != "wild_ann", category != "dom") %>%
  filter(category == "argentina" | category == "niv") %>%
  select(sample.id) %>%unlist(.) %>% as.vector(.) -> pca_samples 

pca <- snpgdsPCA(genofile_all, num.thread = 2, eigen.cnt = 16, sample.id = pca_samples, missing.rate = 0.10, maf = 0.05)

# variance proportion (%)
niv.pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


# make a data.frame
tab.niv <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)

tab.niv <- merge(tab.niv, labels)


#Plot PCA
##########

ggplot(data=tab.niv,aes(EV1,EV2)) + geom_point(aes(color=category)) + ylab("Principal component 2") + xlab("Principal component 1") +
  theme_classic() + theme(legend.position="top",
                          legend.text.align = 0,
                          legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")) +
  scale_color_brewer(name="Type",palette = "Dark2",
                labels=c(expression(paste("Argentinian ",italic("Helianthus"))),
                         expression(paste(italic("H. niveus"))))) +
  guides(col = guide_legend(ncol=1)) +
  theme(panel.border = element_rect(fill = NA, colour = "grey50"))-> niv.plot

#############
#Run PCA for praecox
###############
labels <- read.table("Ana_sample_info_prelim.txt",header=T)
labels$Name <- gsub(".bam", "", labels$Name)
colnames(labels) <- c("sample.id","location","type","plate","category")
labels %>% filter(location != "NA") %>% filter(category != "NA") %>%
  filter(category != "wild_ann", category != "dom") %>%
  filter(category == "argentina" | category == "pra") %>%
  select(sample.id) %>%unlist(.) %>% as.vector(.) -> pca_samples 

pca <- snpgdsPCA(genofile_all, num.thread = 2, eigen.cnt = 16, sample.id = pca_samples, missing.rate = 0.10, maf = 0.05)

# variance proportion (%)
pra.pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


# make a data.frame
tab.pra <- data.frame(sample.id = pca$sample.id,
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      EV3 = pca$eigenvect[,3],
                      EV4 = pca$eigenvect[,4],
                      stringsAsFactors = FALSE)

tab.pra <- merge(tab.pra, labels)


#Plot PCA
##########

ggplot(data=tab.pra,aes(EV1,EV2)) + geom_point(aes(color=category)) + ylab("Principal component 2") + xlab("Principal component 1") +
  theme_classic() + theme(legend.position="top",
                          legend.text.align = 0,
                          legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")) +
  scale_color_brewer(name="Type",palette = "Dark2",
                     labels=c(expression(paste("Argentinian ",italic("Helianthus"))),
                              expression(paste(italic("H. praecox"))))) +
  guides(col = guide_legend(ncol=1)) +
  theme(panel.border = element_rect(fill = NA, colour = "grey50")) -> pra.plot

#############
#Run PCA for debilis
###############
labels <- read.table("Ana_sample_info_prelim.txt",header=T)
labels$Name <- gsub(".bam", "", labels$Name)
colnames(labels) <- c("sample.id","location","type","plate","category")
labels %>% filter(location != "NA") %>% filter(category != "NA") %>%
  filter(category != "wild_ann", category != "dom") %>%
  filter(category == "argentina" | category == "deb") %>%
  select(sample.id) %>%unlist(.) %>% as.vector(.) -> pca_samples 

pca <- snpgdsPCA(genofile_all, num.thread = 2, eigen.cnt = 16, sample.id = pca_samples, missing.rate = 0.10, maf = 0.05)

# variance proportion (%)
deb.pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


# make a data.frame
tab.deb <- data.frame(sample.id = pca$sample.id,
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      EV3 = pca$eigenvect[,3],
                      EV4 = pca$eigenvect[,4],
                      stringsAsFactors = FALSE)

tab.deb <- merge(tab.deb, labels)


#Plot PCA
##########

ggplot(data=tab.deb,aes(EV1,EV2)) + geom_point(aes(color=category)) + ylab("Principal component 2") + xlab("Principal component 1") +
  theme_classic() + theme(legend.position="top",
                          legend.text.align = 0,
                          legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")) +
  scale_color_brewer(name="Type",palette = "Dark2",
                     labels=c(expression(paste("Argentinian ",italic("Helianthus"))),
                              expression(paste(italic("H. debilis"))))) +
  guides(col = guide_legend(ncol=1)) +
  theme(panel.border = element_rect(fill = NA, colour = "grey50")) -> deb.plot

#############
#Run PCA for petiolaris
###############
labels <- read.table("Ana_sample_info_prelim.txt",header=T)
labels$Name <- gsub(".bam", "", labels$Name)
colnames(labels) <- c("sample.id","location","type","plate","category")
labels %>% filter(location != "NA") %>% filter(category != "NA") %>%
  filter(category != "wild_ann", category != "dom") %>%
  filter(category == "argentina" | category == "petpet" | category == "petfal") %>%
  select(sample.id) %>%unlist(.) %>% as.vector(.) -> pca_samples 

pca <- snpgdsPCA(genofile_all, num.thread = 2, eigen.cnt = 16, sample.id = pca_samples, missing.rate = 0.10, maf = 0.05)

# variance proportion (%)
pet.pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


# make a data.frame
tab.pet <- data.frame(sample.id = pca$sample.id,
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      EV3 = pca$eigenvect[,3],
                      EV4 = pca$eigenvect[,4],
                      stringsAsFactors = FALSE)

tab.pet <- merge(tab.pet, labels)


#Plot PCA
##########

ggplot(data=tab.pet,aes(EV1,EV2)) + geom_point(aes(color=category)) + ylab("Principal component 2") + xlab("Principal component 1") +
  theme_classic() + theme(legend.position="top",
                          legend.text.align = 0,
                          legend.background = element_rect(fill="gray90", size=.5, linetype="dotted")) +
  scale_color_brewer(name="Type",palette = "Dark2",
                     labels=c(expression(paste("Argentinian ",italic("Helianthus"))),
                              expression(paste(italic("H. petiolaris fallax"))),
                              expression(paste(italic("H. petiolaris petiolaris"))))) +
  guides(col = guide_legend(ncol=1)) +
  theme(panel.border = element_rect(fill = NA, colour = "grey50")) -> pet.plot

####
#Plot all four PCA together
pdf("Ana.multispecies.PCA.pdf")
multiplot(pet.plot,deb.plot,pra.plot,niv.plot,cols=2)
dev.off()
