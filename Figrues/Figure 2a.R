library(tidyverse)
evec=read.table("ChickenGTEx_wgs_rnaseq_total.filtered.evec")
evec=evec[,1:3]
names(evec)=c("BioSample","PC1","PC2")
rnaseq_samples=read.table("rnaseq_breeds.txt", sep="\t", header=T)
rnaseq_evec=subset(evec, BioSample %in% as.vector(rnaseq_samples$BioSample))
rnaseq_samples %>% mutate(breed=case_when(Breed=="Red jungle fowl" ~ "Red jungle fowl", 
                                           Breed=="Leghorn" ~ "Leghorn", 
                                           Breed=="Rhode Island Red" ~ "Rhode Island Red", 
                                           Breed=="COBB"~"COBB",
                                           Breed=="Tibetan"~"Tibetan",
                                           TRUE   ~ "Other")) -> rnaseq_samples
rnaseq_df=merge(rnaseq_evec, rnaseq_samples, by = "BioSample")
rnaseq_df$breed=factor(rnaseq_df$breed, levels=c("Red jungle fowl","Leghorn", "Rhode Island Red", "COBB", "Tibetan", "Other"))


ggplot(my_df %>% filter(a != "11")) +
  geom_point(aes(x = b, y = c), size = 10, color = "azure3") +
  geom_point(data = my_df %>% filter(a == "11"), aes(x = b, y = c), size = 10, color = "blue")+
  scale_size_manual(values =c(1, 5))+
  theme(legend.position = "none")


p=ggplot(rnaseq_df %>% filter(breed=="Other"))+
geom_point(aes(x=PC1,y=PC2), color = "#BDBDBD", shape=3, size=2)+
geom_point(data=rnaseq_df %>% filter(breed != "Other"),mapping=aes(x=PC1,y=PC2, color = breed, size=breed, shape=breed))+
xlab("PC1 (16.5%)")+ylab("PC2 (6.9%)")+
scale_color_manual(values=c("#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#E6AB02"))+
scale_size_manual(values=c(2,2,2,2,2))+
scale_shape_manual(values=c(19,19,19,19,19))+
theme_linedraw(base_size = 30, base_line_size = 1.5)+
theme(axis.text.x = element_text(color="black", size=25),
        axis.text.y = element_text(color="black", size=25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
pdf("rnaseq_samples_pca.pdf", width=10,height=6)
print(p)
dev.off()

wgs_samples=read.table("wgs_breeds.txt", sep="\t", header=T)
wgs_evec=subset(evec, BioSample %in% as.vector(wgs_samples$BioSampl))
wgs_samples %>% mutate(breed=case_when(Breed=="Red Jungle Fowl" ~ "Red jungle fowl", 
                                           Breed=="Leghorn" ~ "Leghorn", 
                                           Breed=="Rhode Island Red" ~ "Rhode Island Red", 
                                           Breed=="COBB"~"COBB",
                                           Breed=="Tibetan"~"Tibetan",
                                           TRUE   ~ "Other")) -> wgs_samples
wgs_df=merge(wgs_evec, wgs_samples, by = "BioSample")
wgs_df$breed=factor(wgs_df$breed, levels=c("Red jungle fowl","Leghorn", "Rhode Island Red", "COBB", "Tibetan", "Other"))
p=ggplot(wgs_df %>% filter(breed=="Other"))+
geom_point(aes(x=PC1,y=PC2), color = "#BDBDBD", shape=3, size=2)+
geom_point(data=wgs_df %>% filter(breed != "Other"),mapping=aes(x=PC1,y=PC2, color = breed, size=breed, shape=breed))+
xlab("PC1 (16.5%)")+ylab("PC2 (6.9%)")+
scale_color_manual(values=c("#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#E6AB02"))+
scale_size_manual(values=c(2,2,2,2,2))+
scale_shape_manual(values=c(19,19,19,19,19))+
theme_linedraw(base_size = 30, base_line_size = 1.5)+
theme(axis.text.x = element_text(color="black", size=25),
        axis.text.y = element_text(color="black", size=25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf("wgs_samples_pca.pdf", width=10,height=6)
print(p)
dev.off()