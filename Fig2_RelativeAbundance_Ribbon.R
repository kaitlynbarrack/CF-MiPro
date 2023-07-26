---
title: "RelAbun_ribbon"
author: Kaitlyn Barrack
output: html_document
date: '2023-07-26'
---

```{r setup, include=FALSE}

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq", force=TRUE)



library("phyloseq")
library("vegan")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("gplots")
library("ggpubr")
library("dplyr")

```

```{r}

load("CFMiPro_phyloseq")
load("CFMiPro_metadata")

```

```{r}
#convert everything to % relative abundance
CFMiPro_physeq = transform_sample_counts(CFMiPro_phyloseq, function(x) (x / sum(x))*100 )

```

```{r}
#calculate relative abundance of each Phylum 
phyla_counts_tab <- otu_table(tax_glom(CFMiPro_physeq, taxrank="Phylum"))

# making a vector of Phylum names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(CFMiPro_physeq, taxrank="Phylum"))[,2])
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

dim(phyla_counts_tab)

#make a different file for manipulating
major_taxa_for_plot <- as.data.frame(phyla_counts_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Sample"=row.names(CFMiPro_metadata),
                                  "Genotype"=CFMiPro_metadata$Genotype,
                         "Day_of_passage"=CFMiPro_metadata$Day_of_passage,
                                  "Patient" =CFMiPro_metadata$Patient,
                                  "Medium"=CFMiPro_metadata$Medium,
                                  "Replicate"=CFMiPro_metadata$Replicate,
                                  stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)
  
#########################################

#now for Family
family_counts_tab <- otu_table(tax_glom(CFMiPro_physeq, taxrank= "Family"))

fam_tax_vec <- as.vector(tax_table(tax_glom(CFMiPro_physeq, taxrank="Family"))[,5])
rownames(family_counts_tab) <- as.vector(fam_tax_vec)

head(family_counts_tab)

major_taxa_for_plot_fam <- as.data.frame(family_counts_tab)

major_taxa_for_plot_fam$MTaxa <- row.names(major_taxa_for_plot_fam)
major_taxa_for_plot_fam.g <- gather(major_taxa_for_plot_fam, Sample, Proportion, -MTaxa)

head(major_taxa_for_plot_fam.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Genotype"=CFMiPro_metadata$Genotype,
                                  "Sample"=row.names(CFMiPro_metadata), 
                         "Day_of_passage"=CFMiPro_metadata$Day_of_passage, 
                                  "Medium"=CFMiPro_metadata$Medium, 
                                  "Age" =CFMiPro_metadata$Age, 
                                  "Patient" =CFMiPro_metadata$Patient,
                                  "Replicate" =CFMiPro_metadata$Replicate,
                                  stringsAsFactors=F)

major_taxa_for_plot_fam.g2 <- merge(major_taxa_for_plot_fam.g, sample_info_for_merge)
dim(major_taxa_for_plot_fam.g2)

#And genus...

genus_counts_tab <- otu_table(tax_glom(CFMiPro_physeq, taxrank= "Genus"))

genus_tax_vec <- as.vector(tax_table(tax_glom(CFMiPro_physeq, taxrank="Genus"))[,6])
rownames(genus_counts_tab) <- as.vector(genus_tax_vec)

head(genus_counts_tab)

major_taxa_for_plot_genus <- as.data.frame(genus_counts_tab)

major_taxa_for_plot_genus$MTaxa <- row.names(major_taxa_for_plot_genus)
major_taxa_for_plot_genus.g <- gather(major_taxa_for_plot_genus, Sample, Proportion, -MTaxa)

head(major_taxa_for_plot_genus.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Genotype"=CFMiPro_metadata$Genotype,
                                  "Sample"=row.names(CFMiPro_metadata), 
                         "Day_of_passage"=CFMiPro_metadata$Day_of_passage, 
                                  "Medium"=CFMiPro_metadata$Medium, 
                                  "Age" =CFMiPro_metadata$Age, 
                                  "Patient" =CFMiPro_metadata$Patient,
                                  "Replicate" =CFMiPro_metadata$Replicate,
                                  stringsAsFactors=F)


major_taxa_for_plot_genus.g2 <- merge(major_taxa_for_plot_genus.g, sample_info_for_merge)
dim(major_taxa_for_plot_genus.g2)



```

```{r}
#Phylum level line graph 

#subset on phyla of interest
new_gplot <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa =="Proteobacteria" |
                                      major_taxa_for_plot.g2$MTaxa =="Firmicutes" |
                                      major_taxa_for_plot.g2$MTaxa =="Bacteroidota" |
                                      major_taxa_for_plot.g2$MTaxa =="Actinobacteriota" |
                                      major_taxa_for_plot.g2$MTaxa =="Verrucomicrobiota",]

library(tibble)
install.packages("readxl")
install.packages("writexl")
library(readxl)
library(writexl)
install.packages("openxlsx")
library(openxlsx)
install.packages("viridis")
library(viridis)

write.xlsx(new_gplot, 'CFMiPro_16s_relativeabundance_Phylum.xlsx')

#Day 0 CF - Table S3
CF_phylum <- subset(new_gplot, Genotype == 'CF')
CF0_phylum <- subset(CF_phylum, Day_of_passage == '0')
CF0_phylum_agr <- aggregate(Proportion ~ MTaxa+Medium, data=CF0_phylum, FUN=mean)

#Day 1-5 CF - Table S3
CF15_phylum <- subset(CF_phylum, Day_of_passage %in% c('1','2','3','4','5'))
CF15_phylum_agr <- aggregate(Proportion ~ MTaxa+Medium, data=CF15_phylum, FUN=mean)

#Day 0 nonCF - Table S3
nonCF_phylum <- subset(new_gplot, Genotype == 'nonCF')
nonCF0_phylum <- subset(nonCF_phylum, Day_of_passage == '0')
nonCF0_phylum_agr <- aggregate(Proportion ~ MTaxa+Medium, data=nonCF0_phylum, FUN=mean)

#Day 1-5 nonCF - Table S3
nonCF15_phylum <- subset(nonCF_phylum, Day_of_passage %in% c('1','2','3','4','5'))
nonCF15_phylum_agr <- aggregate(Proportion ~ MTaxa+Medium, data=nonCF15_phylum, FUN=mean)

#ribbon plot

Phylum_ribbon <- ggplot(new_gplot, aes(x=Day_of_passage, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), linewidth =2) + theme_bw() +
  labs(x = "Day of passage", y = "Relative Abundance (%)") +
  xlim(0,5) +
  scale_color_viridis(discrete=TRUE, name = "Phylum", labels = c("Actinobacteria",
                                                                 "Bacteroidetes",
                                                                 "Firmicutes",
                                                                 "Proteobacteria",
                                                                 "Verrucomicrobiota")) +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=10,face="bold"), 
        legend.text=element_text(size=10), 
        legend.title=element_text(size=10, face="bold"),
        legend.position = ,
        plot.margin=unit(c(0.4,1,0.4,0.4),"cm"),
        legend.key.size = unit(0.7, 'cm'))+
  facet_wrap(~Genotype+Medium)

Phylum_ribbon #Figure 2A


#Subset to see just MiPro among all patients - Figure S4
new_gplot_MiPro <- subset(new_gplot, Medium == 'MiPro')

Phylum_MiPro_ribbon <- ggplot(new_gplot_MiPro, aes(x=Day_of_passage, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), linewidth =2) + theme_bw() 
  labs(x = "Day of passage", y = "Relative Abundance (%)") +
  xlim(0,5) +
  scale_color_viridis(discrete=TRUE, name = "Phylum", labels = c("Actinobacteria",
                                                                 "Bacteroidetes",
                                                                 "Firmicutes",
                                                                 "Proteobacteria",
                                                                 "Verrucomicrobiota")) +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=10,face="bold"), 
        legend.text=element_text(size=10), 
        legend.title=element_text(size=10, face="bold"),
        legend.position = ,
        plot.margin=unit(c(0.4,1,0.4,0.4),"cm"),
        legend.key.size = unit(0.7, 'cm'))+
  facet_wrap(~Genotype+Patient)

Phylum_MiPro_ribbon



#Family level line graph

new_gplot_fam <- major_taxa_for_plot_fam.g2[major_taxa_for_plot_fam.g2$MTaxa =="Enterobacteriaceae" |
                                      major_taxa_for_plot_fam.g2$MTaxa =="Clostridiaceae" |
                                      major_taxa_for_plot_fam.g2$MTaxa =="Lactobacillaceae" |
                                      major_taxa_for_plot_fam.g2$MTaxa =="Akkermansiaceae" |
                                      major_taxa_for_plot_fam.g2$MTaxa =="Pseudomonadaceae" |
                                      major_taxa_for_plot_fam.g2$MTaxa =="Lachnospiraceae" |
                                      major_taxa_for_plot_fam.g2$MTaxa =="Bifidobacteriaceae" |
                                      major_taxa_for_plot_fam.g2$MTaxa =="Ruminococcaceae" |
                                      major_taxa_for_plot_fam.g2$MTaxa =="Bacteroidaceae" |
                                      major_taxa_for_plot_fam.g2$MTaxa =="Veillonellaceae",]


write.xlsx(new_gplot_fam, 'CFMiPro_16s_relativeabundance_top10Family.xlsx')

#Day 0 CF - Table S8
CF_fam <- subset(new_gplot_fam, Genotype == 'CF')
CF0_fam <- subset(CF_fam, Day_of_passage == '0')
CF0_fam_agr <- aggregate(Proportion ~ MTaxa+Medium, data=CF0_fam, FUN=mean)

#Day 1-5 CF - Table S8
CF15_fam <- subset(CF_fam, Day_of_passage %in% c('1','2','3','4','5'))
CF15_fam_agr <- aggregate(Proportion ~ MTaxa+Medium, data=CF15_fam, FUN=mean)

#Day 0 nonCF - Table S8
nonCF_fam <- subset(new_gplot_fam, Genotype == 'nonCF')
nonCF0_fam <- subset(nonCF_fam, Day_of_passage == '0')
nonCF0_fam_agr <- aggregate(Proportion ~ MTaxa+Medium, data=nonCF0_fam, FUN=mean)

#Day 1-5 nonCF - Table S8
nonCF15_fam <- subset(nonCF_fam, Day_of_passage %in% c('1','2','3','4','5'))
nonCF15_fam_agr <- aggregate(Proportion ~ MTaxa+Medium, data=nonCF15_fam, FUN=mean)


Family_ribbon <- ggplot(new_gplot_fam, aes(x=Day_of_passage, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), linewidth =2) + theme_bw() +
  labs(x = "Day of passage", y = "Relative Abundance (%)") +
  xlim(0,5) +
  ylim(-25,100)+
  scale_color_viridis(discrete=TRUE, name = "Family") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=10,face="bold"), 
        legend.text=element_text(size=10), 
        legend.title=element_text(size=10, face="bold"),
        legend.position = ,
        plot.margin=unit(c(0.4,1,0.4,0.4),"cm"),
        legend.key.size = unit(0.7, 'cm')) +
  facet_wrap(~Genotype+Medium)

Family_ribbon #Figure 2B

#Genus level line graph

new_gplot_CFgenus <- major_taxa_for_plot_genus.g2[major_taxa_for_plot_genus.g2$MTaxa =="Escherichia.Shigella" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Clostridium.sensu.stricto.1" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Lactobacillus" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Klebsiella" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Akkermansia" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Citrobacter" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Veillonella" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Bacteroides" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Bifidobacterium" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Streptococcus" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Enterococcus" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Blautia" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Pseudomonas" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Roseburia" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="X.Ruminococcus..gnavus.group" |
                                      major_taxa_for_plot_genus.g2$MTaxa =="Faecalibacterium",]


write.xlsx(new_gplot_CFgenus, 'CFMiPro_16s_relativeabundance_16CFGenera.xlsx')

#Day 0 CF - Table S9
Day0_genus <- subset(new_gplot_CFgenus, Day_of_passage == '0')
CF_Day0_genus <- subset(Day0_genus, Genotype == 'CF')
CF_Day0_genus_agr <- aggregate(Proportion ~ MTaxa+Medium, data=CF_Day0_genus, FUN=mean)

#Day 1-5 CF - Table S9
CF_genus <- subset(new_gplot_CFgenus, Genotype == 'CF')
CF15_genus <- subset(CF_genus, Day_of_passage %in% c('1','2','3','4','5'))
CF15_genus_agr <- aggregate(Proportion ~ MTaxa+Medium, data=CF15_genus, FUN=mean)

#Day 0 CF - Table S9
nonCF_genus <- subset(new_gplot_CFgenus, Genotype == 'nonCF')
nonCF_Day0_genus <- subset(nonCF_genus, Day_of_passage == '0')
nonCF_Day0_genus_agr <- aggregate(Proportion ~ MTaxa+Medium, data=nonCF_Day0_genus, FUN=mean)

#Day 1-5 CF - Table S9
nonCF15_genus <- subset(nonCF_genus, Day_of_passage %in% c('1','2','3','4','5'))
nonCF15_genus_agr <- aggregate(Proportion ~ MTaxa+Medium, data=nonCF15_genus, FUN=mean)


#Genus ribbon - Figure S5

CFGenus_ribbon <- ggplot(new_gplot_CFgenus, aes(x=Day_of_passage, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), linewidth =2) + theme_bw() +
  labs(x = "Day of passage", y = "Relative Abundance (%)") +
  xlim(0,5) +
  ylim(-25,100) +
  scale_color_viridis(discrete=TRUE, name = "Genus") +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=10,face="bold"), 
        legend.text=element_text(size=10), 
        legend.title=element_text(size=10, face="bold"),
        legend.position = ,
        plot.margin=unit(c(0.4,1,0.4,0.4),"cm"),
        legend.key.size = unit(0.7, 'cm')) +
  facet_wrap(~Genotype+Medium)

CFGenus_ribbon

######

