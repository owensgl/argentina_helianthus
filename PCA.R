library(ggplot2)
library(SNPRelate)
library("bigmemory")
library(dplyr)
library("mapdata")
library("maps")
library("maptools")
library("rworldmap")
library(ggmap)
#Load genotype matrix
geno_matrix <- read.big.matrix("/media/owens/Childs1/Ana/Ana.GATK.2016.80.5.ntab", sep="\t", header = TRUE)
#geno_matrix <- read.big.matrix("/media/owens/Childs1/Ana/Ana.GATK.2016.plusEST.depth10.90.ntab", sep="\t", header = TRUE)
geno_matrix  <- geno_matrix[,-c(1:2)]


genotype_id_df <- read.table(pipe("cut -f1,2 /media/owens/Childs1/Ana/Ana.GATK.2016.80.5.ntab"), header = TRUE)
#genotype_id_df <- read.table(pipe("cut -f1,2 /media/owens/Childs1/Ana/Ana.GATK.2016.plusEST.depth10.90.ntab"), header = TRUE)
genotype_id_df$id <- paste(genotype_id_df$CHROM,genotype_id_df$POS,sep="_")
sample_id <- colnames(geno_matrix)

snpgdsCreateGeno("/media/owens/Childs1/Ana/Ana.GATK.2016.80.5.gds", genmat = geno_matrix, 
#snpgdsCreateGeno("/media/owens/Childs1/Ana/Ana.GATK.2016.plusEST.depth10.90.gds", genmat = geno_matrix, 
                sample.id = sample_id, snpfirstdim = TRUE, 
                 snp.id = genotype_id_df$id, snp.chromosome = genotype_id_df$CHROM, 
                 snp.position = genotype_id_df$POS)

rm(list=c("geno_matrix"))



###################
#Load sample info
###################
labels <- read.table("/home/owens/bin/Ana/Ana_sample_info_prelim.txt",header=T)
labels <- labels %>% filter(plate == "1" | plate == "2") %>% filter(Location != "<NA>")
samplelist <- labels$Name
##################
#Load gds file
###################

genofile_all <- snpgdsOpen("/media/owens/Childs1/Ana/Ana.GATK.2016.80.5.gds", allow.duplicate = TRUE)
#genofile_all <- snpgdsOpen("/media/owens/Childs1/Ana/Ana.GATK.2016.plusEST.depth10.90.gds", allow.duplicate = TRUE)

snpset_pruned <- snpgdsLDpruning(genofile_all)
snpset.id <- unlist(snpset_pruned)
#############
#Run PCA
###############

pca_samples <- sample_id[!grepl(bad_samples, sample_id)]

pca <- snpgdsPCA(genofile_all, num.thread = 2, eigen.cnt = 16, snp.id=snpset.id, sample.id = samplelist, missing.rate = 0.10, maf = 0.05)

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Gather sample information

colnames(labels) <- c("sample.id","location","type","plate")
locations <- read.table("Ana_locations.txt",header=T)

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)

tab <- merge(tab, labels)
tab <- merge(tab,locations)

###########
#Plot PCA
##########
pdf("Ana.GATK.2016.80.5.PCA.Nov4.pdf")
ggplot(data=tab,aes(EV1,EV2)) + geom_point(aes(color=location)) + 
  ylab(paste("Principal component 2 (",round(pc.percent[2],2),"% PVE)",sep="")) + 
  xlab(paste("Principal component 1 (",round(pc.percent[1],2),"% PVE)",sep="")) +
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_label(aes(label = "H. annuus", x = -.12, y = -.15), hjust = 0, vjust = 0, color="black", size=3.5,fontface = "italic") +
  geom_label(aes(label = "H. petiolaris", x = .06, y = -.15), hjust = 0, vjust = 0, color="black", size=3.5,fontface = "italic") +
  geom_label(aes(label = "Off-type", x = -.02, y = -.15), hjust = 0, vjust = 0, color="black", size=3.5) +
  ylim(c(-.15,.2))
  
dev.off()

#####
#Plot pure petiolaris on a map. 
pet.pops <- tab %>% filter(EV1 > 0.05) %>% group_by(location) %>% summarise(lat=mean(lat),long=mean(long),EV2=mean(EV2)) %>%
  mutate(type = ifelse(EV2 > 0, 1,2))
map <- get_map(location = c(lon = -64, lat = -36), zoom=7)

largermap <- get_map(location = c(lon = -64, lat = -36), zoom=5)
bigmap <- ggmap(largermap) +geom_rect(aes(xmin=-67,xmax=-61,ymin=-39,ymax=-33),fill=NA,color="black") +
  xlab("Longitude") + ylab("Latitude")
smallmap <- ggmap(map) + geom_point(data=pet.pops, aes(y=lat,x=long, color=as.factor(type)),size=4,alpha=0.8)+ 
  scale_color_manual(values=c("red", "blue"), 
                    name="Population group",
                    breaks=c("1", "2" ),
                    labels=c("Population 1", "Population 2")) +
  xlab("Longitude") + ylab("Latitude")

pdf("Ana.GATK.2016.80.5.PC2onmap.Nov4.pdf",height = 6,width=10)
multiplot(bigmap,smallmap,cols=2)
dev.off()

####Determine pure petiolaris and hybrid groups
#Pure Petiolaris 1
tab %>% filter(EV1 > 0.05, EV2 > 0) %>% select(sample.id)

#Hybrid Petiolaris 1
tab %>% filter(EV1 < 0.05, EV1 > -0.05,EV2 > 0) %>% select(sample.id)

#Pure Petiolaris 2
tab %>% filter(EV1 > 0.05, EV2 < 0) %>% select(sample.id)

#Hybrid Petiolaris 2
tab %>% filter(EV1 < 0.05, EV1 > -0.05,EV2 < 0) %>% select(sample.id) %>% write.table(., file="")

tab$hiest_type <- NULL
tab$pet_clade <- NULL
#label tab file
for (i in 1:nrow(tab)){
  if (tab$EV1[i] > 0.05 & tab$EV2[i] > 0){
    tab$hiest_type[i] <- "P2"
    tab$pet_clade[i] <- "1"
  }else if(tab$EV1[i] < 0.05 & tab$EV1[i] > -0.05 & tab$EV2[i] > 0){
    tab$hiest_type[i] <- "H"
    tab$pet_clade[i] <- "1"
  }else if(tab$EV1[i] > 0.05 & tab$EV2[i] < 0){
    tab$hiest_type[i] <- "P2"
    tab$pet_clade[i] <- "2"
  }else if (tab$EV1[i] < 0.05 & tab$EV1[i] > -0.05 & tab$EV2[i] < 0){
    tab$hiest_type[i] <- "H"
    tab$pet_clade[i] <- "2"
  }else if (tab$EV1[i] < -0.05){
    tab$hiest_type[i] <- "P1"
    tab$pet_clade[i] <- "0"
  }
}

write.table(tab, file="/media/owens/Childs1/Ana/Ana.PCAdata.full.txt",quote=F, row.names=F)
tab %>% filter(hiest_type == "P1" | pet_clade == 1) %>% select(sample.id, hiest_type) %>% 
  write.table(., file="/media/owens/Childs1/Ana/Ana.hiestlist.pet1.txt",quote=F,sep="\t",row.names=F, col.names=F)

tab %>% filter(hiest_type == "P1" | pet_clade == 2) %>% select(sample.id, hiest_type) %>% 
  write.table(., file="/media/owens/Childs1/Ana/Ana.hiestlist.pet2.txt",quote=F,sep="\t",row.names=F, col.names=F)
