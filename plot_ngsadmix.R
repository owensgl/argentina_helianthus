library(ggplot2)
library(dplyr)
library(reshape)
source("multiplot.R")
data2 <- read.table("/home/owens/working/Ana/Ana.ngsadmix.2.qopt",colClasses="numeric")
colnames(data2) <- c("G1","G2")
labels <- read.table("/home/owens/bin/Ana/Ana_sample_info_prelim.txt",header=T)
labels <- labels %>% filter(plate == "1" | plate == "2")

data2 <- cbind(data2, labels)

data2 <- melt(data2, id.var=c("Name","Location","Type","plate"))
data2$value <- as.numeric(data2$value)

data2$Name <- factor(data2$Name, 
                         levels=data2[order(data2$Type),]$Name)


levels(data2$Location)
levels(data2$Location) <- c("BAR", "CAT","CHU","CZ","HIL","SAL","SAN","UNI","Deb","DIA","Dom","Neg","Niv","NivCan","NivTep","Pet","PetFal","PetPet","Pra","RCU","USAwild","WIN")

pdf("Ana_NGSadmix_k2.pdf",width=10,height=6)
data2 %>%
  filter(Location != "NA") %>%
  group_by(Name) %>%
  mutate(major_k_qval = max(value)) %>%
  mutate(G1_qval = value[variable == "G1"]) %>%
  mutate(major_k = variable[value == max(value)]) %>% 
  ungroup %>%
  arrange(major_k, desc(G1_qval)) %>%
  ungroup %>%
  mutate(Name = factor(Name, levels = as.character(Name))) %>%
  ggplot(aes(x = Name, y = value, fill = factor(variable)))+
  geom_bar(stat = "identity", width = 1.1)+
  facet_grid(~Location, scales = "free", space = "free",switch="x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        panel.margin = unit(0.1, "lines"), 
        strip.background = element_blank(),
        strip.switch.pad.grid = unit(0.0, "lines"))+
  theme(plot.margin = unit(c(0,0.5,0.5,0.0), "cm")) +
  ylab("q-value")+
  scale_fill_manual(values=c("#123eabff","#ff0000ff"),
                    name="Groups",
                    labels=c("Annuus (1)", "Petiolaris (2)")) +
  theme(strip.text.x = element_text(angle=90,size=6)) -> k2
dev.off()

#####plot with k=3
data <- read.table("/home/owens/working/Ana/Ana.ngsadmix.3.qopt",colClasses="numeric")
colnames(data) <- c("G1","G2","G3")
labels <- read.table("Ana_sample_info_prelim.txt",header=T)
labels <- labels %>% filter(plate == "1" | plate == "2")

data <- cbind(data, labels)

data <- melt(data, id.var=c("Name","Location","Type","plate"))
data$value <- as.numeric(data$value)

data$Name <- factor(data$Name, 
                    levels=data[order(data$value),]$Name)

levels(data$Location)
levels(data$Location) <- c("BAR", "CAT","CHU","CZ","HIL","SAL","SAN","UNI","Deb","DIA","Dom","Neg","Niv","NivCan","NivTep","Pet","PetFal","PetPet","Pra","RCU","USAwild","WIN")


pdf("Ana_NGSadmix_k3.pdf",width=10,height=6)
data %>%
  filter(Location != "NA") %>%
  group_by(Name) %>%
  mutate(major_k_qval = max(value)) %>%
  mutate(G1_qval = value[variable == "G1"]) %>%
  mutate(G2_qval = value[variable == "G2"]) %>%
  mutate(major_k = variable[value == max(value)]) %>% 
  ungroup %>%
  arrange(desc(G1_qval),desc(G2_qval)) %>%
  ungroup %>%
  mutate(Name = factor(Name, levels = as.character(Name))) %>%
  ggplot(aes(x = Name, y = value, fill = factor(variable)))+
  geom_bar(stat = "identity", width = 1.1)+
  facet_grid(~Location, scales = "free", space = "free", switch = "x") +
  theme_classic()+
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        panel.margin = unit(0.1, "lines"), 
        strip.background = element_blank(),
        strip.switch.pad.grid = unit(0.0, "lines"))+
  theme(plot.margin = unit(c(0,0.5,0.5,0.0), "cm")) +
  ylab("q-value")+
  scale_fill_manual(values=c("#123eabff","#ff0000ff","dark green"),
                      name="Groups",
                    labels=c("Annuus (1)", "Petiolaris (2)", "Petiolaris (3)")) +
  theme(strip.text.x = element_text(angle=90,size=6)) -> k3

dev.off()



pdf("Ana_NGSadmix.pdf",width=10,height=6)
multiplot(k2,k3, cols=1)
dev.off()

