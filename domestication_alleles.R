########
#Add in domestication allele data
dom <- read.delim("/home/owens/working/Ana/Ana.GATK.2016.plusEST.depth10.domesticationfixed.countdiverged.txt",header=T)
colnames(dom) <- c("Name","P1","P2","total_diagnostic","percent_dom")

dom <- merge(labels, dom)
dom <- dom %>% filter(plate == "1" | plate == "2")
colnames(dom)[1] <- c("Name")
dom$P1 <- NULL
dom$P2 <- NULL
dom$total_diagnostic <- NULL

data.dom <- read.table("/home/owens/working/Ana/Ana.ngsadmix.2.qopt",colClasses="numeric")
colnames(data.dom) <- c("G1","G2")
data.dom <- cbind(data.dom, dom)

data.dom <- melt(data.dom, id.var=c("Name","Location","Type","plate","percent_dom"))
data.dom$value <- as.numeric(data.dom$value)

data.dom$Name <- factor(data.dom$Name, 
                        levels=data.dom[order(data.dom$Type),]$Name)


pdf("Ana_NGSadmix_k2_plusdom.pdf",width=10,height=6)
data %>%
  group_by(Name) %>%
  mutate(major_k_qval = max(value)) %>%
  mutate(major_k = variable[value == max(value)]) %>% 
  ungroup %>%
  arrange(major_k, desc(major_k_qval)) %>%
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
  theme(plot.margin = unit(c(0,0.5,0.0,0.25), "cm")) +
  ylab("q-value")+
  scale_fill_manual(values=c("#123eabff","#ff0000ff"),
                    name="K values",
                    labels=c("Annuus (1)", "Petiolaris (2)")) +
  geom_point(aes(x=Name,y=percent_dom))

dev.off()




###########Plot domestic alleles by type

dom <- read.delim("/home/owens/working/Ana/Ana.GATK.2016.plusEST.depth10.domesticationfixed.countdiverged.txt",header=T)
colnames(dom) <- c("Name","P1","P2","total_diagnostic","percent_dom")

dom <- merge(labels, dom)
dom <- dom %>% filter(plate == "1" | plate == "2")
colnames(dom)[1] <- c("Name")
dom$P1 <- NULL
dom$P2 <- NULL
dom$total_diagnostic <- NULL

data.dom <- read.table("/home/owens/working/Ana/Ana.ngsadmix.2.qopt",colClasses="numeric")
colnames(data.dom) <- c("G1","G2")
data.dom <- cbind(data.dom, dom)


data.dom$genetictype <- NULL
for (i in 1:nrow(data.dom)){
  if (data.dom$G1[i] > 0.75){
    data.dom$genetictype[i] <- "A"
  }else if (data.dom$G1[i] < 0.25){
    data.dom$genetictype[i] <- "P"
  }else{
    data.dom$genetictype[i] <- "H"
  }
}

levels(data.dom$Location)
levels(data.dom$Location) <- c("BAR", "CAT","CHU","CZ","HIL","SAL","SAN","UNI","Deb","DIA","Dom","Neg","Niv","NivCan","NivTep","Pet","PetFal","PetPet","Pra","RCU","USAwild","WIN")

pdf("Ana_domestic_alleles.pdf")
data.dom %>% filter(Location != "NA") %>%
  ggplot(.,aes(x=Location,y=percent_dom,color=genetictype)) + geom_boxplot(fill="grey") + theme_classic() +
  ylab("Proportion of domestic alleles") +
  theme(legend.text.align = 0) +
  scale_color_manual(values=c("#123eabff","purple", "#ff0000ff"), 
                     name="Genetic group",
                     breaks=c("A", "H", "P"),
                     labels=c(expression(paste(italic("H. Annuus"))),
                              "Hybrid",
                              expression(paste(italic("H. Petiolaris")))))
dev.off()


data.dom %>% filter(Location == "CZ") %>%group_by(genetictype)%>% 
  summarise(mean(percent_dom))



######
dom %>% filter(Location != "NA") -> dom.used
range(dom.used$total_diagnostic)
mean(dom.used$total_diagnostic)
