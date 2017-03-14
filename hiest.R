library(HIest)
library(dplyr)

G <-read.table('/home/owens/working/Ana/Ana.GATK.2016.depth10.fixeddif.hiest.G', sep = '\t',header=T, row.names=1)
P <-read.table('/home/owens/working/Ana/Ana.GATK.2016.depth10.fixeddif.hiest.P', sep = '\t', header=T)

est <- HIest(G,P,method="SANN",type="allele.count",iterations=1000,surf=F,startgrid=99)

est$sample <- rownames(G)
labels <- read.table("Ana_sample_info_prelim.txt",header=T)
labels$sample <- gsub(".bam", "", labels$Name)
est.labels <- merge(est, labels)
brewercolor <- brewer.pal(n = 1, name = "Set1")

pdf("Ana.GATK.2016.fixeddif.hiest.Jul13.pdf")
est.labels %>% filter(Location != "NA") %>%
ggplot(., aes(x=S, y=H)) +geom_segment(aes(x = 0, y = 0, xend = .5, yend = 1)) +
  geom_segment(aes(x=0,y=0,xend=1,yend=0)) + ggtitle("Argentinian sunflowers") +
  geom_segment(aes(x = 1, y = 0, xend = .5, yend = 1)) +geom_point(size=5,alpha=0.5) +
  theme_bw() + geom_point(aes(x=0.54, y=0.64),color="#E41A1C",size=4) +
  geom_label(aes(x = 0.1, y = -.03, label= "H. petiolaris", fontface="italic")) +
  geom_label(aes(x = 0.9, y = -.03, label= "H. annuus", fontface="italic")) +
  geom_label(aes(x = 0.9, y = 0.64, label= "F1"),color="#E41A1C") +
  geom_segment(aes(x = 0.8, y = 0.64, xend = 0.6, yend = 0.64),color="#E41A1C",
               arrow = arrow(length = unit(0.03, "npc")))

dev.off()

expression(paste("H. petiolaris"))