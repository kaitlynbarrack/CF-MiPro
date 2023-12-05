---
title: "RelAbun_ribbon"
author: Kaitlyn Barrack, adapted from Courtney Price (DOI: https://doi.org/10.1128/msphere.00046-23)
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


#Subset to see just MiPro among all patients - Figure S6
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

#Day 0 CF - Table S11
Day0_genus <- subset(new_gplot_CFgenus, Day_of_passage == '0')
CF_Day0_genus <- subset(Day0_genus, Genotype == 'CF')
CF_Day0_genus_agr <- aggregate(Proportion ~ MTaxa+Medium, data=CF_Day0_genus, FUN=mean)

#Day 1-5 CF - Table S11
CF_genus <- subset(new_gplot_CFgenus, Genotype == 'CF')
CF15_genus <- subset(CF_genus, Day_of_passage %in% c('1','2','3','4','5'))
CF15_genus_agr <- aggregate(Proportion ~ MTaxa+Medium, data=CF15_genus, FUN=mean)

#Day 0 CF - Table S11
nonCF_genus <- subset(new_gplot_CFgenus, Genotype == 'nonCF')
nonCF_Day0_genus <- subset(nonCF_genus, Day_of_passage == '0')
nonCF_Day0_genus_agr <- aggregate(Proportion ~ MTaxa+Medium, data=nonCF_Day0_genus, FUN=mean)

#Day 1-5 CF - Table S11
nonCF15_genus <- subset(nonCF_genus, Day_of_passage %in% c('1','2','3','4','5'))
nonCF15_genus_agr <- aggregate(Proportion ~ MTaxa+Medium, data=nonCF15_genus, FUN=mean)


#Genus ribbon - Figure S7

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


#Linear model for each phylum:

install.packages("readr")
library(readr)

#Actinobacteria

Actino_phylum <- read_csv("7-20-23_CFMiPro_16S_PhylumRA_Actinobacteriota.csv") 
Actino_phylum$Medium <- factor(Actino_phylum$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Actino_phylum_CF <- subset(Actino_phylum, Genotype == 'CF')
Actino_phylum_nonCF <- subset(Actino_phylum, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Actino_phylum_CF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)       7.1331     1.2367  37.9799   5.768 1.19e-06 ***
#Mediumlow_CF      0.3170     1.1133 150.0000   0.285    0.776    
#Mediummedian_CF   0.2192     1.1133 150.0000   0.197    0.844    
#Day_of_passage   -1.5444     0.2661 150.0000  -5.803 3.74e-08 ***

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Actino_phylum_nonCF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)       6.9319     1.3511  22.1478   5.131 3.77e-05 ***
#Mediumlow_CF      0.3147     0.8395 218.0000   0.375    0.708    
#Mediummedian_CF   0.7484     0.8395 218.0000   0.892    0.374    
#Day_of_passage   -1.3932     0.2007 218.0000  -6.942 4.36e-11 ***

#Bacteroidota

Bact_phylum <- read_csv("7-20-23_CFMiPro_16S_PhylumRA_Bacteroidota.csv") 
Bact_phylum$Medium <- factor(Bact_phylum$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Bact_phylum_CF <- subset(Bact_phylum, Genotype == 'CF')
Bact_phylum_nonCF <- subset(Bact_phylum, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Bact_phylum_CF))

#Fixed effects:
                # Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)       5.11714    2.27821   9.58604   2.246 0.049594 *  
#Mediumlow_CF     -1.89005    0.80721 150.00000  -2.341 0.020525 *  
#Mediummedian_CF  -2.78199    0.80721 150.00000  -3.446 0.000737 ***
#Day_of_passage    0.01092    0.19296 150.00000   0.057 0.954943 

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Bact_phylum_nonCF))

#Fixed effects:
                 #Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)      14.09497    1.66672  22.83122   8.457 1.73e-08 ***
#Mediumlow_CF     -0.02252    1.05749 218.00000  -0.021    0.983    
#Mediummedian_CF  -4.33624    1.05749 218.00000  -4.100 5.82e-05 ***
#Day_of_passage   -0.28312    0.25279 218.00000  -1.120    0.264

#Firmicutes

Firm_phylum <- read_csv("7-20-23_CFMiPro_16S_PhylumRA_Firmicutes.csv") 
Firm_phylum$Medium <- factor(Firm_phylum$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Firm_phylum_CF <- subset(Firm_phylum, Genotype == 'CF')
Firm_phylum_nonCF <- subset(Firm_phylum, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Firm_phylum_CF))

#Fixed effects:
               # Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)      72.7637     9.0237  10.7431   8.064 7.03e-06 ***
#Mediumlow_CF    -11.9189     4.0289 150.0000  -2.958   0.0036 ** 
#Mediummedian_CF  -1.8411     4.0289 150.0000  -0.457   0.6484    
#Day_of_passage   -0.1176     0.9631 150.0000  -0.122   0.9030   

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Firm_phylum_nonCF))

#Fixed effects:
               # Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)      51.3741     4.9277  26.6951  10.426 6.56e-11 ***
#Mediumlow_CF    -11.5477     3.4279 218.0000  -3.369 0.000893 ***
#Mediummedian_CF -10.3009     3.4279 218.0000  -3.005 0.002967 ** 
#Day_of_passage   -4.3670     0.8194 218.0000  -5.329 2.45e-07 ***

#Proteobacteria

Prot_phylum <- read_csv("7-20-23_CFMiPro_16S_PhylumRA_Proteobacteria.csv") 
Prot_phylum$Medium <- factor(Prot_phylum$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Prot_phylum_CF <- subset(Prot_phylum, Genotype == 'CF')
Prot_phylum_nonCF <- subset(Prot_phylum, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Prot_phylum_CF))

#Fixed effects:
               # Estimate Std. Error      df t value Pr(>|t|)   
#(Intercept)       13.731      8.268  11.837   1.661   0.1230   
#Mediumlow_CF      14.007      4.209 150.000   3.328   0.0011 **
#Mediummedian_CF    4.936      4.209 150.000   1.173   0.2428   
#Day_of_passage     1.833      1.006 150.000   1.822   0.0705 .  

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Prot_phylum_nonCF))

#Fixed effects:
               # Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)      19.6395     4.5022  41.3792   4.362 8.37e-05 ***
#Mediumlow_CF     19.9110     3.7265 218.0000   5.343 2.30e-07 ***
#Mediummedian_CF  23.1570     3.7265 218.0000   6.214 2.59e-09 ***
#Day_of_passage    5.2769     0.8908 218.0000   5.924 1.21e-08 ***

#Verrucomicrobiota

Verr_phylum <- read_csv("7-20-23_CFMiPro_16S_PhylumRA_Verrucomicrobiota.csv") 
Verr_phylum$Medium <- factor(Verr_phylum$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Verr_phylum_CF <- subset(Verr_phylum, Genotype == 'CF')
Verr_phylum_nonCF <- subset(Verr_phylum, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Verr_phylum_CF))

#Fixed effects:
                 #Estimate Std. Error        df t value Pr(>|t|)  
#(Intercept)       1.07852    0.42549  15.84952   2.535   0.0222 *
#Mediumlow_CF     -0.44545    0.27706 150.00000  -1.608   0.1100  
#Mediummedian_CF  -0.44337    0.27706 150.00000  -1.600   0.1116  
#Day_of_passage   -0.16572    0.06623 150.00000  -2.502   0.0134 * 

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Verr_phylum_nonCF))

#Fixed effects:
               # Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)       6.7197     1.6280  37.0227   4.128  0.00020 ***
#Mediumlow_CF     -8.3626     1.2998 218.0000  -6.434 7.78e-10 ***
#Mediummedian_CF  -8.8991     1.2998 218.0000  -6.847 7.58e-11 ***
#Day_of_passage    0.9288     0.3107 218.0000   2.989  0.00312 ** 


#Linear model for each family:

install.packages("readr")
library(readr)

#Akkermansiaceae

Akk_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Akkermansia.csv") 
Akk_family$Medium <- factor(Akk_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Akk_family_CF <- subset(Akk_family, Genotype == 'CF')
Akk_family_nonCF <- subset(Akk_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Akk_family_CF))

#Fixed effects:
                # Estimate Std. Error        df t value Pr(>|t|)  
#(Intercept)       1.07850    0.42549  15.84921   2.535   0.0222 *
#Mediumlow_CF     -0.44560    0.27706 150.00000  -1.608   0.1099  
#Mediummedian_CF  -0.44337    0.27706 150.00000  -1.600   0.1116  
#Day_of_passage   -0.16571    0.06623 150.00000  -2.502   0.0134 *

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Akk_family_nonCF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)       6.6816     1.6293  36.9861   4.101 0.000217 ***
#Mediumlow_CF     -8.3471     1.3004 218.0000  -6.419 8.45e-10 ***
#Mediummedian_CF  -8.8876     1.3004 218.0000  -6.835 8.12e-11 ***
#Day_of_passage    0.9314     0.3108 218.0000   2.996 0.003049 **

#Bacteroidaceae

Bacter_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Bacteroidaceae.csv") 
Bacter_family$Medium <- factor(Bacter_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Bacter_family_CF <- subset(Bacter_family, Genotype == 'CF')
Bacter_family_nonCF <- subset(Bacter_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Bacter_family_CF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)       4.2442     2.1443   9.3656   1.979   0.0779 .  
#Mediumlow_CF     -1.7333     0.7112 150.0000  -2.437   0.0160 *  
#Mediummedian_CF  -2.9293     0.7112 150.0000  -4.119 6.27e-05 ***
#Day_of_passage    0.1759     0.1700 150.0000   1.035   0.3026  

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Bacter_family_nonCF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)      11.6226     1.4119  26.0547   8.232 1.02e-08 ***
#Mediumlow_CF      0.5203     0.9697 218.0000   0.537 0.592094    
#Mediummedian_CF  -3.7881     0.9697 218.0000  -3.907 0.000125 ***
#Day_of_passage   -0.1222     0.2318 218.0000  -0.527 0.598621    

#Bifidobacteriaceae

Bif_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Bifido.csv") 
Bif_family$Medium <- factor(Bif_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Bif_family_CF <- subset(Bif_family, Genotype == 'CF')
Bif_family_nonCF <- subset(Bif_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Bif_family_CF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)       3.6244     1.0295  42.2094   3.520 0.001048 ** 
#Mediumlow_CF      0.9751     0.9486 150.0000   1.028 0.305673    
#Mediummedian_CF   0.9217     0.9486 150.0000   0.972 0.332834    
#Day_of_passage   -0.9026     0.2268 150.0000  -3.980 0.000107 ***

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Bif_family_nonCF))

#Fixed effects:
               # Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)       4.7780     1.3219  21.6400   3.615  0.00157 ** 
#Mediumlow_CF      0.7400     0.8077 218.0000   0.916  0.36057    
#Mediummedian_CF   1.6068     0.8077 218.0000   1.989  0.04791 *  
#Day_of_passage   -1.1434     0.1931 218.0000  -5.922 1.23e-08 ***

#Clostridiaceae

Clost_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Clost.csv") 
Clost_family$Medium <- factor(Clost_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Clost_family_CF <- subset(Clost_family, Genotype == 'CF')
Clost_family_nonCF <- subset(Clost_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Clost_family_CF))

#Fixed effects:
                 #Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)      37.41135    8.82358  17.08520   4.240 0.000546 ***
#Mediumlow_CF     -9.24874    6.00414 150.00000  -1.540 0.125571    
#Mediummedian_CF   0.02504    6.00414 150.00000   0.004 0.996678    
#Day_of_passage    3.53774    1.43526 150.00000   2.465 0.014835 *

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Clost_family_nonCF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)   
#(Intercept)       5.2682     5.2921  14.9844   0.995  0.33531   
#Mediumlow_CF      2.6689     2.0681 218.0000   1.290  0.19825   
#Mediummedian_CF   1.2844     2.0681 218.0000   0.621  0.53521   
#Day_of_passage    1.3487     0.4944 218.0000   2.728  0.00689 **

#Enterobacteriaceae

Ent_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Entero.csv") 
Ent_family$Medium <- factor(Ent_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Ent_family_CF <- subset(Ent_family, Genotype == 'CF')
Ent_family_nonCF <- subset(Ent_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Ent_family_CF))

#Fixed effects:
               # Estimate Std. Error      df t value Pr(>|t|)   
#(Intercept)       11.318      7.901  12.897   1.432  0.17580   
#Mediumlow_CF      14.360      4.398 150.000   3.265  0.00136 **
#Mediummedian_CF    5.455      4.398 150.000   1.240  0.21679   
#Day_of_passage     2.291      1.051 150.000   2.179  0.03088 * 

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Ent_family_nonCF))

#Fixed effects:
               # Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)       14.486      4.316  48.633   3.356  0.00154 ** 
#Mediumlow_CF      18.904      3.740 218.000   5.055 9.12e-07 ***
#Mediummedian_CF   25.228      3.740 218.000   6.746 1.35e-10 ***
#Day_of_passage     5.550      0.894 218.000   6.208 2.67e-09 ***

#Lachnospiraceae

Lach_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Lachno.csv") 
Lach_family$Medium <- factor(Lach_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Lach_family_CF <- subset(Lach_family, Genotype == 'CF')
Lach_family_nonCF <- subset(Lach_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Lach_family_CF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)      26.9144     3.4098  30.6494   7.893 7.11e-09 ***
#Mediumlow_CF     -5.8155     2.9088 150.0000  -1.999   0.0474 *  
#Mediummedian_CF  -6.8195     2.9088 150.0000  -2.344   0.0204 *  
#Day_of_passage   -5.0439     0.6953 150.0000  -7.254 2.01e-11 ***

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Lach_family_nonCF))

#Fixed effects:
               # Estimate Std. Error      df t value Pr(>|t|)    
#(Intercept)       25.415      2.166  96.463  11.736  < 2e-16 ***
#Mediumlow_CF      -5.961      2.159 218.000  -2.761  0.00625 ** 
#Mediummedian_CF   -6.393      2.159 218.000  -2.962  0.00340 ** 
#Day_of_passage    -4.851      0.516 218.000  -9.399  < 2e-16 ***

#Lactobacillaceae

Lact_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Lacto.csv") 
Lact_family$Medium <- factor(Lact_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Lact_family_CF <- subset(Lact_family, Genotype == 'CF')
Lact_family_nonCF <- subset(Lact_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Lact_family_CF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)      -5.5689     5.5118  13.8563  -1.010 0.329645    
#Mediumlow_CF      8.8619     3.2646 150.0000   2.715 0.007415 ** 
#Mediummedian_CF  12.0666     3.2646 150.0000   3.696 0.000306 ***
#Day_of_passage    2.3985     0.7804 150.0000   3.073 0.002514 ** 

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Lact_family_nonCF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)   
#(Intercept)       0.4349     0.7280  67.1735   0.597  0.55223   
#Mediumlow_CF      0.1820     0.6795 218.0000   0.268  0.78905   
#Mediummedian_CF   2.1967     0.6795 218.0000   3.233  0.00142 **
#Day_of_passage   -0.1217     0.1624 218.0000  -0.749  0.45465   

#Pseudomonadaceae

Pseudo_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Pseudo.csv") 
Pseudo_family$Medium <- factor(Pseudo_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Pseudo_family_CF <- subset(Pseudo_family, Genotype == 'CF')
Pseudo_family_nonCF <- subset(Pseudo_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Pseudo_family_CF))

#Fixed effects:
                  #Estimate Std. Error         df t value Pr(>|t|)   
#(Intercept)      1.621e+00  7.434e-01  3.599e+01   2.181  0.03583 * 
#Mediumlow_CF     5.518e-04  6.609e-01  1.500e+02   0.001  0.99933   
#Mediummedian_CF -1.919e-02  6.609e-01  1.500e+02  -0.029  0.97688   
#Day_of_passage  -4.305e-01  1.580e-01  1.500e+02  -2.725  0.00719 **

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Pseudo_family_nonCF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)  
#(Intercept)       0.8243     0.8443  55.0044   0.976   0.3332  
#Mediumlow_CF      1.4048     0.7544 218.0000   1.862   0.0639 .
#Mediummedian_CF  -0.1953     0.7544 218.0000  -0.259   0.7960  
#Day_of_passage   -0.1444     0.1803 218.0000  -0.801   0.4240   

#Ruminococcaceae

Rum_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Rumino.csv") 
Rum_family$Medium <- factor(Rum_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Rum_family_CF <- subset(Rum_family, Genotype == 'CF')
Rum_family_nonCF <- subset(Rum_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Rum_family_CF))

#Fixed effects:
                 #Estimate Std. Error        df t value Pr(>|t|)   
#(Intercept)       0.96002    0.27927  40.56226   3.438  0.00137 **
#Mediumlow_CF     -0.18020    0.25514 150.00000  -0.706  0.48111   
#Mediummedian_CF  -0.19388    0.25514 150.00000  -0.760  0.44851   
#Day_of_passage   -0.20103    0.06099 150.00000  -3.296  0.00122 **

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Rum_family_nonCF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)       4.2613     0.5351  88.1216   7.963 5.46e-12 ***
#Mediumlow_CF     -0.1791     0.5254 218.0000  -0.341    0.734    
#Mediummedian_CF  -0.2134     0.5254 218.0000  -0.406    0.685    
#Day_of_passage   -1.0939     0.1256 218.0000  -8.710 7.67e-16 ***   

#Veillonellaceae

Veil_family <- read_csv("7-20-23_CFMiPro_16S_relativeabundance_Family_Veillonellaceae.csv") 
Veil_family$Medium <- factor(Veil_family$Medium, levels=c("MiPro", "low_CF", "median_CF"))

Veil_family_CF <- subset(Veil_family, Genotype == 'CF')
Veil_family_nonCF <- subset(Veil_family, Genotype == 'nonCF')

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Veil_family_CF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)    
#(Intercept)       1.1091     0.2502  58.9241   4.433 4.11e-05 ***
#Mediumlow_CF     -0.1985     0.2456 150.0000  -0.808 0.420213    
#Mediummedian_CF  -0.2483     0.2456 150.0000  -1.011 0.313563    
#Day_of_passage   -0.2120     0.0587 150.0000  -3.612 0.000413 ***

summary(lmer(Proportion ~ Medium + Day_of_passage + (1 | Patient), data=Veil_family_nonCF))

#Fixed effects:
                #Estimate Std. Error       df t value Pr(>|t|)   
#(Intercept)       0.2188     1.2118  30.2647   0.181  0.85794   
#Mediumlow_CF      0.5450     0.8949 218.0000   0.609  0.54317   
#Mediummedian_CF   2.4719     0.8949 218.0000   2.762  0.00623 **
#Day_of_passage    0.4321     0.2139 218.0000   2.020  0.04463 * 
                                         

