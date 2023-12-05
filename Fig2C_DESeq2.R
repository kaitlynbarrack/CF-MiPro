---
title: "Figure 2C: Differentially-abundant genera between Media"
author: "Kaitlyn Barrack, adapted from Courtney Price (DOI: https://doi.org/10.1128/msphere.00046-23)"
date: '2023-07-26'
output: html_document
---

```{r}
#This script first removes ASVs with <1% prevalence
#Then combines taxa at the genus level 
#Then uses DESeq2 to call differential abundance at the genus level
#modified from preprocess_vst.R

library("phyloseq")
library("ggplot2")
library("plyr")
library("DESeq2")
library("vegan")
library("viridis")

##load in phyloseq data 
load("CFMiPro_phyloseq")
load("CFMiPro_metadata")

#Subset phyloseq objects by genotype:
load("phyloseq_nonCF")
load("phyloseq_CF")

#pull out Days 1-5 in each genotype 
phyloseq_nonCF15 <- subset_samples(phyloseq_nonCF, Day_of_passage %in% c('1','2','3', '4', '5'))

phyloseq_CF15 <- subset_samples(phyloseq_CF, Day_of_passage %in% c('1','2','3', '4', '5'))

CFMiPro_phyloseq_15 <- subset_samples(CFMiPro_phyloseq, Day_of_passage %in% c('1','2','3', '4', '5'))



##############################################################################################################################################
#Run this for nonCF first:

#Preprocessing for DESEq
# Compute prevalence (how many samples each taxa appears in) of each feature, store as data.frame
prevdf = apply(X = otu_table(phyloseq_nonCF15),
               MARGIN = ifelse(taxa_are_rows(phyloseq_nonCF15), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(phyloseq_nonCF15),
                    tax_table(phyloseq_nonCF15))

#a plot of ASV taxa abundances 
prevdf = subset(prevdf, Phylum %in% get_taxa_unique(phyloseq_nonCF15, "Phylum"))
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(phyloseq_nonCF15),color=Phylum)) + 
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + 
  facet_wrap(~Phylum) + theme_bw() + guides(color = FALSE, size = FALSE)

# Define prevalence threshold as 5% of total samples
prevalenceThreshold = 0.05 * nsamples(phyloseq_nonCF15)
prevalenceThreshold

# Remove taxa with <1% prevalence
keepTaxa = rownames(prevdf)[(prevdf$Prevalence >= prevalenceThreshold)]
phyloseq_2 = prune_taxa(keepTaxa, phyloseq_nonCF15)
##########################################################################################################
#significance testing at the genera level 
#merge ASV_physeq_2 at the genus level 
CFMiPro_phyloseq_3 <- tax_glom(phyloseq_2, taxrank="Genus")

#convert to deseq object
#This will model using non-continuous Medium, but will control for patient
#Note to self: this is right. If you don't believe it, check the bioconductor vingette 'analyzing RNA-seq data with DESeq2', 'Multi-factor designs' section
ASV_deseq_genus <- phyloseq_to_deseq2(CFMiPro_phyloseq_3, ~Patient + Medium)

#workaround to deal with 0s
cts <- counts(ASV_deseq_genus)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
ASV_deseq_genus <- estimateSizeFactors(ASV_deseq_genus, geoMeans=geoMeans)

#deseq standard analysis
ASV_deseq_genus <- DESeq(ASV_deseq_genus)

##########################################################################################################

# pulling out our results table, we specify the object, the p-value we are going to use to filter our results 
#and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_Medium_nonCF15 <- results(ASV_deseq_genus_nonCF15, alpha=0.05, contrast=c("Medium","median_CF", "MiPro"))

# we can get a glimpse at what this table currently holds with the summary command
summary(deseq_res_Medium_nonCF15)

# stitch together with the ASV's taxonomic annotations for a quick look at both together
deseq_res_Medium_tax_nonCF15 <- cbind(as(deseq_res_Medium_nonCF15, "data.frame"), as(tax_table(phyloseq_nonCF15_2)[row.names(deseq_res_Medium_nonCF15), ], "matrix"))

# and sort that table by the baseMean column
sum_Medium_nonCF15 <- deseq_res_Medium_tax_nonCF15[order(deseq_res_Medium_tax_nonCF15$baseMean, decreasing=T), ]

write.csv(sum_Medium_nonCF15, "sum_Medium_genus_5%p_nonCF15.csv")

# subset this table to only include these that pass our specified significance level
sum_Medium_sig <- sum_Medium[which(sum_Medium$padj < 0.05), ]
write.csv(sum_Medium_nonCF15_sig, "sum_Medium_5%p_0.05pval_genus_nonCF15.csv")

##########################################################################################################
#I did some manual work here to identify whether significantly altered taxa were relevant to CF
#Under CFChange column
#classified taxa as CF-associated ("CF"), or not CF-associated ("None")

##########################################################################################################
#Making graphs
#DESeq2 visualization methods from phyloseq tutorial
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

#read in manually annotated csv
medCFvMi_plot <- read.csv("sum_Medium_5%p_0.05pval_genus_nonCF15_updated.csv", row.names = 1, header=TRUE)
#Assign column names as the titles in row 1
colnames(medCFvMi_plot) <- medCFvMi_plot[1, ]  
#Drop row 1
medCFvMi_plot <- medCFvMi_plot[-1,]


# Genus order medCFvMi
x = tapply(medCFvMi_plot$log2FoldChange, medCFvMi_plot$Genus, function(x) max(x))
x = sort(x, TRUE)
medCFvMi_plot$Genus = factor(as.character(medCFvMi_plot$Genus), levels=names(x))

gen_medCFvMi <- ggplot(medCFvMi_plot, aes(x=log2FoldChange, y=Genus, color=CFChange)) + 
  geom_point(size=5) + 
  scale_color_viridis(discrete=TRUE) +
  labs(x="Log2 Fold Change", y="Genus",color="CF-Associated Taxa") +
  geom_vline(xintercept = 0) +
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.title=element_text(size=12,face="bold"), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12, face="bold"))+
  geom_label(label="Decreased in CF-MiPro", x=-1.5, y=9, label.size = 0.4, color = "black", size =3.5)+
  geom_label(label="Increased in CF-MiPro", x= 1.1, y=9, label.size = 0.4, color = "black", size =3.5)

#Figure 2C
gen_medCFvMi

#Comparing raw cultures (Day0) to Day 5 of MiPro or Day 5 of medCF
#Separate by MiPro and median-CF-MiPro, all days (use phyloseq object that has Day as characters)
load(CF_phyloseq_chday)

phyloseq_chday_nonCF <- subset_samples(CF_phyloseq_chday, Genotype == 'nonCF')
phyloseq_chday_CF <- subset_samples(CF_phyloseq_chday, Genotype == 'CF')

phyloseq_nonCF_MiPro <- subset_samples(phyloseq_chday_nonCF, Medium == 'MiPro')
phyloseq_nonCF_medCF <- subset_samples(phyloseq_chday_nonCF, Medium == 'median_CF')

phyloseq_CF_MiPro <- subset_samples(phyloseq_chday_CF, Medium == 'MiPro')
phyloseq_CF_medCF <- subset_samples(phyloseq_chday_CF, Medium == 'median_CF')

#Run this for nonCF MiPro first:
#Preprocessing for DESEq
# Compute prevalence (how many samples each taxa appears in) of each feature, store as data.frame
prevdf_nonCFMiPro = apply(X = otu_table(phyloseq_nonCF_MiPro),
               MARGIN = ifelse(taxa_are_rows(phyloseq_nonCF_MiPro), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf_nonCFMiPro = data.frame(Prevalence = prevdf_nonCFMiPro,
                    TotalAbundance = taxa_sums(phyloseq_nonCF_MiPro),
                    tax_table(phyloseq_nonCF_MiPro))

#a plot of ASV taxa abundances 
prevdf_nonCFMiPro = subset(prevdf_nonCFMiPro, Phylum %in% get_taxa_unique(phyloseq_nonCF_MiPro, "Phylum"))

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(phyloseq_nonCF_MiPro),color=Phylum)) + 
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + 
  facet_wrap(~Phylum) + theme_bw() + guides(color = FALSE, size = FALSE)

# Define prevalence threshold as 5% of total samples
prevalenceThreshold_nonCF_MiPro = 0.05 * nsamples(phyloseq_nonCF_MiPro)
prevalenceThreshold_nonCF_MiPro #3.9

# Remove taxa with <1% prevalence
keepTaxa_nonCF_MiPro = rownames(prevdf_nonCFMiPro)[(prevdf_nonCFMiPro$Prevalence >= prevalenceThreshold_nonCF_MiPro)]
phyloseq_nonCF_MiPro_2 = prune_taxa(keepTaxa_nonCF_MiPro, phyloseq_nonCF_MiPro)
##########################################################################################################
#significance testing at the genera level 
#merge ASV_physeq_2 at the genus level 
phyloseq_nonCF_MiPro_3 <- tax_glom(phyloseq_nonCF_MiPro_2, taxrank="Genus")

#convert to deseq object
#This will model using non-continuous Medium, but will control for patient
#Note to self: this is right. If you don't believe it, check the bioconductor vingette 'analyzing RNA-seq data with DESeq2', 'Multi-factor designs' section
ASV_deseq_genus_nonCF_MiPro <- phyloseq_to_deseq2(phyloseq_nonCF_MiPro_3, ~Patient + ch_Day_of_passage)

#workaround to deal with 0s
cts_nonCF_MiPro <- counts(ASV_deseq_genus_nonCF_MiPro)
geoMeans_nonCF_MiPro <- apply(cts_nonCF_MiPro, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
ASV_deseq_genus_nonCF_MiPro <- estimateSizeFactors(ASV_deseq_genus_nonCF_MiPro, geoMeans=geoMeans_nonCF_MiPro)

#deseq standard analysis
ASV_deseq_genus_nonCF_MiPro <- DESeq(ASV_deseq_genus_nonCF_MiPro)

##########################################################################################################

# pulling out our results table, we specify the object, the p-value we are going to use to filter our results 
#and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_Day_nonCF_MiPro <- results(ASV_deseq_genus_nonCF_MiPro, alpha=0.05, contrast=c("ch_Day_of_passage","Day_5", "Day_0"))

# we can get a glimpse at what this table currently holds with the summary command
summary(deseq_res_Day_nonCF_MiPro)

# stitch together with the ASV's taxonomic annotations for a quick look at both together
deseq_res_Day_tax_nonCF_MiPro <- cbind(as(deseq_res_Day_nonCF_MiPro, "data.frame"), as(tax_table(phyloseq_nonCF_MiPro_2)[row.names(deseq_res_Day_nonCF_MiPro), ], "matrix"))

# and sort that table by the baseMean column
sum_Day_nonCF_MiPro <- deseq_res_Day_tax_nonCF_MiPro[order(deseq_res_Day_tax_nonCF_MiPro$baseMean, decreasing=T), ]

write.csv(sum_Day_nonCF_MiPro, "sum_Day_genus_5%p_nonCF_MiPro.csv")

# subset this table to only include these that pass our specified significance level
sum_Day_nonCF_MiPro_sig <- sum_Day_nonCF_MiPro[which(sum_Day_nonCF_MiPro$padj < 0.05), ]
write.csv(sum_Day_nonCF_MiPro_sig, "sum_Day_5%p_0.05pval_nonCF_MiPro_genus_all.csv")

##########################################################################################################
#I did some manual work here to identify whether significantly altered taxa were relevant to CF
#Under CFChange column
#classified taxa as CF-associated ("CF"), or not CF associated ("None")

##########################################################################################################
#Making graphs
#DESeq2 visualization methods from phyloseq tutorial
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

#read in manually annotated csv
nonCF_MiPro_day_plot <- read.csv("sum_Day_5%p_0.05pval_nonCF_MiPro_genus_all.csv", row.names = 1, header=TRUE)

#Assign column names as the titles in row 1
#colnames(nonCF_MiPro_day_plot) <- nonCF_MiPro_day_plot[1, ]  
#Drop row 1
#medCFvMi_plot <- medCFvMi_plot[-1,]


# Genus order nonCF MiPro Day
x_nonCF_MiPro = tapply(nonCF_MiPro_day_plot$log2FoldChange, nonCF_MiPro_day_plot$Genus, function(x_nonCF_MiPro) max(x_nonCF_MiPro))
x_nonCF_MiPro = sort(x_nonCF_MiPro, TRUE)
nonCF_MiPro_day_plot$Genus = factor(as.character(nonCF_MiPro_day_plot$Genus), levels=names(x_nonCF_MiPro))

gen_nonCF_MiPro <- ggplot(nonCF_MiPro_day_plot, aes(x=log2FoldChange, y=Genus))+ #color=CFChange)) + 
  geom_point(size=5) + 
  scale_color_viridis(discrete=TRUE) +
  labs(x="Log2 Fold Change", y="Genus")) +
  geom_vline(xintercept = 0) +
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.title=element_text(size=12,face="bold"), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12, face="bold")) +
 geom_label(label="Decreased in MiPro", x=-15, y=24, label.size = 0.4, color = "black",  size =3.5)+
 geom_label(label="Increased in MiPro", x= 13, y=24, label.size = 0.4, color = "black", size =3.5)


#Figure S8A
gen_nonCF_MiPro

##########

#Comparing raw cultures (Day0) to Day 5 of MiPro or Day 5 of medCF
#Separate by MiPro and median-CF-MiPro, all days (use phyloseq object that has Day as characters)
load(CF_phyloseq_chday)

phyloseq_chday_nonCF <- subset_samples(CF_phyloseq_chday, Genotype == 'nonCF')
phyloseq_chday_CF <- subset_samples(CF_phyloseq_chday, Genotype == 'CF')

phyloseq_nonCF_MiPro <- subset_samples(phyloseq_chday_nonCF, Medium == 'MiPro')
phyloseq_nonCF_medCF <- subset_samples(phyloseq_chday_nonCF, Medium == 'median_CF')

phyloseq_CF_MiPro <- subset_samples(phyloseq_chday_CF, Medium == 'MiPro')
phyloseq_CF_medCF <- subset_samples(phyloseq_chday_CF, Medium == 'median_CF')

#Running this for: CF MiPro
#Preprocessing for DESEq
# Compute prevalence (how many samples each taxa appears in) of each feature, store as data.frame
prevdf_CF_MiPro = apply(X = otu_table(phyloseq_CF_MiPro),
               MARGIN = ifelse(taxa_are_rows(phyloseq_CF_MiPro), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf_CF_MiPro = data.frame(Prevalence = prevdf_CFMiPro,
                    TotalAbundance = taxa_sums(phyloseq_CF_MiPro),
                    tax_table(phyloseq_CF_MiPro))

#a plot of ASV taxa abundances 
prevdf_CF_MiPro = subset(prevdf_CFMiPro, Phylum %in% get_taxa_unique(phyloseq_CF_MiPro, "Phylum"))

ggplot(prevdf_CF_MiPro, aes(TotalAbundance, Prevalence / nsamples(phyloseq_CF_MiPro),color=Phylum)) + 
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + 
  facet_wrap(~Phylum) + theme_bw() + guides(color = FALSE, size = FALSE)

# Define prevalence threshold as 5% of total samples
prevalenceThreshold_CF_MiPro = 0.05 * nsamples(phyloseq_CF_MiPro)
prevalenceThreshold_CF_MiPro #2.7

# Remove taxa with <1% prevalence
keepTaxa_CF_MiPro = rownames(prevdf_CF_MiPro)[(prevdf_CF_MiPro$Prevalence >= prevalenceThreshold_CF_MiPro)]
phyloseq_CF_MiPro_2 = prune_taxa(keepTaxa_CF_MiPro, phyloseq_CF_MiPro)


##########################################################################################################

#significance testing at the genera level 
#merge ASV_physeq_2 at the genus level 
phyloseq_CF_MiPro_3 <- tax_glom(phyloseq_CF_MiPro_2, taxrank="Genus")

#convert to deseq object
#This will model using non-continuous Medium, but will control for patient
#Note to self: this is right. If you don't believe it, check the bioconductor vingette 'analyzing RNA-seq data with DESeq2', 'Multi-factor designs' section
ASV_deseq_genus_CF_MiPro <- phyloseq_to_deseq2(phyloseq_CF_MiPro_3, ~Patient + ch_Day_of_passage)

#workaround to deal with 0s
cts_CF_MiPro <- counts(ASV_deseq_genus_CF_MiPro)
geoMeans_CF_MiPro <- apply(cts_CF_MiPro, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
ASV_deseq_genus_CF_MiPro <- estimateSizeFactors(ASV_deseq_genus_CF_MiPro, geoMeans=geoMeans_CF_MiPro)

#deseq standard analysis
ASV_deseq_genus_CF_MiPro <- DESeq(ASV_deseq_genus_CF_MiPro)

##########################################################################################################

# pulling out our results table, we specify the object, the p-value we are going to use to filter our results 
#and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_Day_CF_MiPro <- results(ASV_deseq_genus_CF_MiPro, alpha=0.05, contrast=c("ch_Day_of_passage","Day_5", "Day_0"))

# we can get a glimpse at what this table currently holds with the summary command
summary(deseq_res_Day_CF_MiPro)

# stitch together with the ASV's taxonomic annotations for a quick look at both together
deseq_res_Day_tax_CF_MiPro <- cbind(as(deseq_res_Day_CF_MiPro, "data.frame"), as(tax_table(phyloseq_CF_MiPro_2)[row.names(deseq_res_Day_CF_MiPro), ], "matrix"))

# and sort that table by the baseMean column
sum_Day_CF_MiPro <- deseq_res_Day_tax_CF_MiPro[order(deseq_res_Day_tax_CF_MiPro$baseMean, decreasing=T), ]

write.csv(sum_Day_CF_MiPro, "sum_Day_genus_5%p_CF_MiPro.csv")

# subset this table to only include these that pass our specified significance level
sum_Day_CF_MiPro_sig <- sum_Day_CF_MiPro[which(sum_Day_CF_MiPro$padj < 0.05), ]

write.csv(sum_Day_CF_MiPro_sig, "sum_Day_5%p_0.05pval_CF_MiPro_genus_all.csv")

#Making graphs
#DESeq2 visualization methods from phyloseq tutorial
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

#read in manually annotated csv
CF_MiPro_day_plot <- read.csv("sum_Day_5%p_0.05pval_CF_MiPro_genus_all.csv", row.names = 1, header=TRUE)

#Assign column names as the titles in row 1
#colnames(nonCF_MiPro_day_plot) <- nonCF_MiPro_day_plot[1, ]  
#Drop row 1
#medCFvMi_plot <- medCFvMi_plot[-1,]


# Genus order nonCF MiPro Day
x_CF_MiPro = tapply(CF_MiPro_day_plot$log2FoldChange, CF_MiPro_day_plot$Genus, function(x_CF_MiPro) max(x_CF_MiPro))
x_CF_MiPro = sort(x_CF_MiPro, TRUE)
CF_MiPro_day_plot$Genus = factor(as.character(CF_MiPro_day_plot$Genus), levels=names(x_CF_MiPro))

gen_CF_MiPro <- ggplot(CF_MiPro_day_plot, aes(x=log2FoldChange, y=Genus))+ #color=CFChange)) + 
  geom_point(size=5) + 
  scale_color_viridis(discrete=TRUE) +
  labs(x="Log2 Fold Change", y="Genus")) +
  geom_vline(xintercept = 0) +
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.title=element_text(size=12,face="bold"), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12, face="bold")) +
 geom_label(label="Decreased in MiPro", x=-7, y=21, label.size = 0.4, color = "black",  size =3.5)+
 geom_label(label="Increased in MiPro", x= 12, y=21, label.size = 0.4, color = "black", size =3.5)



#Figure S8B
gen_CF_MiPro

```

```{r}

#Comparing raw cultures (Day0) to Day 5 of MiPro or Day 5 of medCF
#Separate by MiPro and median-CF-MiPro, all days (use phyloseq object that has Day as characters)
load(CF_phyloseq_chday)

phyloseq_chday_nonCF <- subset_samples(CF_phyloseq_chday, Genotype == 'nonCF')
phyloseq_chday_CF <- subset_samples(CF_phyloseq_chday, Genotype == 'CF')

phyloseq_nonCF_MiPro <- subset_samples(phyloseq_chday_nonCF, Medium == 'MiPro')
phyloseq_nonCF_medCF <- subset_samples(phyloseq_chday_nonCF, Medium == 'median_CF')

phyloseq_CF_MiPro <- subset_samples(phyloseq_chday_CF, Medium == 'MiPro')
phyloseq_CF_medCF <- subset_samples(phyloseq_chday_CF, Medium == 'median_CF')

#Running for medCF of nonCF samples (Day5) compared to raw
#Preprocessing for DESEq
# Compute prevalence (how many samples each taxa appears in) of each feature, store as data.frame
prevdf_nonCF_medCF = apply(X = otu_table(phyloseq_nonCF_medCF),
               MARGIN = ifelse(taxa_are_rows(phyloseq_nonCF_medCF), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf_nonCF_medCF = data.frame(Prevalence = prevdf_nonCF_medCF,
                    TotalAbundance = taxa_sums(phyloseq_nonCF_medCF),
                    tax_table(phyloseq_nonCF_medCF))

#a plot of ASV taxa abundances 
prevdf_nonCF_medCF = subset(prevdf_nonCF_medCF, Phylum %in% get_taxa_unique(phyloseq_nonCF_medCF, "Phylum"))

ggplot(prevdf_nonCF_medCF, aes(TotalAbundance, Prevalence / nsamples(phyloseq_nonCF_medCF),color=Phylum)) + 
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + 
  facet_wrap(~Phylum) + theme_bw() + guides(color = FALSE, size = FALSE)

# Define prevalence threshold as 5% of total samples
prevalenceThreshold_nonCF_medCF = 0.05 * nsamples(phyloseq_nonCF_medCF)
prevalenceThreshold_nonCF_medCF #3.9

# Remove taxa with <1% prevalence
keepTaxa_nonCF_medCF = rownames(prevdf_nonCF_medCF)[(prevdf_nonCF_medCF$Prevalence >= prevalenceThreshold_nonCF_medCF)]
phyloseq_nonCF_medCF_2 = prune_taxa(keepTaxa_nonCF_medCF, phyloseq_nonCF_medCF)
##########################################################################################################
#significance testing at the genera level 
#merge ASV_physeq_2 at the genus level 
phyloseq_nonCF_medCF_3 <- tax_glom(phyloseq_nonCF_medCF_2, taxrank="Genus")

#convert to deseq object
#This will model using non-continuous Medium, but will control for patient
#Note to self: this is right. If you don't believe it, check the bioconductor vingette 'analyzing RNA-seq data with DESeq2', 'Multi-factor designs' section
ASV_deseq_genus_nonCF_medCF <- phyloseq_to_deseq2(phyloseq_nonCF_medCF_3, ~Patient + ch_Day_of_passage)

#workaround to deal with 0s
cts_nonCF_medCF <- counts(ASV_deseq_genus_nonCF_medCF)
geoMeans_nonCF_medCF <- apply(cts_nonCF_medCF, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
ASV_deseq_genus_nonCF_medCF <- estimateSizeFactors(ASV_deseq_genus_nonCF_medCF, geoMeans=geoMeans_nonCF_medCF)

#deseq standard analysis
ASV_deseq_genus_nonCF_medCF <- DESeq(ASV_deseq_genus_nonCF_medCF)

##########################################################################################################

# pulling out our results table, we specify the object, the p-value we are going to use to filter our results 
#and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_Day_nonCF_medCF <- results(ASV_deseq_genus_nonCF_medCF, alpha=0.05, contrast=c("ch_Day_of_passage","Day_5", "Day_0"))

# we can get a glimpse at what this table currently holds with the summary command
summary(deseq_res_Day_nonCF_medCF)

# stitch together with the ASV's taxonomic annotations for a quick look at both together
deseq_res_Day_tax_nonCF_medCF <- cbind(as(deseq_res_Day_nonCF_medCF, "data.frame"), as(tax_table(phyloseq_nonCF_medCF_2)[row.names(deseq_res_Day_nonCF_medCF), ], "matrix"))

# and sort that table by the baseMean column
sum_Day_nonCF_medCF <- deseq_res_Day_tax_nonCF_medCF[order(deseq_res_Day_tax_nonCF_medCF$baseMean, decreasing=T), ]

write.csv(sum_Day_nonCF_medCF, "sum_Day_genus_5%p_nonCF_medCF.csv")

# subset this table to only include these that pass our specified significance level
sum_Day_nonCF_medCF_sig <- sum_Day_nonCF_medCF[which(sum_Day_nonCF_medCF$padj < 0.05), ]
write.csv(sum_Day_nonCF_medCF_sig, "sum_Day_5%p_0.05pval_nonCF_medCF_genus_all.csv")

#Making graphs
#DESeq2 visualization methods from phyloseq tutorial
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

#read in manually annotated csv
nonCF_medCF_day_plot <- read.csv("sum_Day_5%p_0.05pval_nonCF_medCF_genus_all.csv", row.names = 1, header=TRUE)

#Assign column names as the titles in row 1
#colnames(nonCF_MiPro_day_plot) <- nonCF_MiPro_day_plot[1, ]  
#Drop row 1
#medCFvMi_plot <- medCFvMi_plot[-1,]


# Genus order nonCF MiPro Day
x_nonCF_medCF = tapply(nonCF_medCF_day_plot$log2FoldChange, nonCF_medCF_day_plot$Genus, function(x_nonCF_medCF) max(x_nonCF_medCF))
x_nonCF_medCF = sort(x_nonCF_medCF, TRUE)
nonCF_medCF_day_plot$Genus = factor(as.character(nonCF_medCF_day_plot$Genus), levels=names(x_nonCF_medCF))

gen_nonCF_medCF <- ggplot(nonCF_medCF_day_plot, aes(x=log2FoldChange, y=Genus))+ #color=CFChange)) + 
  geom_point(size=5) + 
  scale_color_viridis(discrete=TRUE) +
  labs(x="Log2 Fold Change", y="Genus")) +
  geom_vline(xintercept = 0) +
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.title=element_text(size=12,face="bold"), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12, face="bold")) +
 geom_label(label="Decreased in medCF", x=-10, y=14, label.size = 0.4, color = "black",  size =3.5)+
 geom_label(label="Increased in medCF", x= 4, y=14, label.size = 0.4, color = "black", size =3.5)

#Figure S8C
gen_nonCF_medCF


```

```{r}


#Comparing raw cultures (Day0) to Day 5 of MiPro or Day 5 of medCF
#Separate by MiPro and median-CF-MiPro, all days (use phyloseq object that has Day as characters)
load(CF_phyloseq_chday)

phyloseq_chday_nonCF <- subset_samples(CF_phyloseq_chday, Genotype == 'nonCF')
phyloseq_chday_CF <- subset_samples(CF_phyloseq_chday, Genotype == 'CF')

phyloseq_nonCF_MiPro <- subset_samples(phyloseq_chday_nonCF, Medium == 'MiPro')
phyloseq_nonCF_medCF <- subset_samples(phyloseq_chday_nonCF, Medium == 'median_CF')

phyloseq_CF_MiPro <- subset_samples(phyloseq_chday_CF, Medium == 'MiPro')
phyloseq_CF_medCF <- subset_samples(phyloseq_chday_CF, Medium == 'median_CF')

#Running for medCF of CF samples (Day5) compared to raw
#Preprocessing for DESEq
# Compute prevalence (how many samples each taxa appears in) of each feature, store as data.frame
prevdf_CF_medCF = apply(X = otu_table(phyloseq_CF_medCF),
               MARGIN = ifelse(taxa_are_rows(phyloseq_CF_medCF), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf_CF_medCF = data.frame(Prevalence = prevdf_CF_medCF,
                    TotalAbundance = taxa_sums(phyloseq_CF_medCF),
                    tax_table(phyloseq_CF_medCF))

#a plot of ASV taxa abundances 
prevdf_CF_medCF = subset(prevdf_CF_medCF, Phylum %in% get_taxa_unique(phyloseq_CF_medCF, "Phylum"))

ggplot(prevdf_CF_medCF, aes(TotalAbundance, Prevalence / nsamples(phyloseq_CF_medCF),color=Phylum)) + 
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + 
  facet_wrap(~Phylum) + theme_bw() + guides(color = FALSE, size = FALSE)

# Define prevalence threshold as 5% of total samples
prevalenceThreshold_CF_medCF = 0.05 * nsamples(phyloseq_CF_medCF)
prevalenceThreshold_CF_medCF 
#2.7

# Remove taxa with <1% prevalence
keepTaxa_CF_medCF = rownames(prevdf_CF_medCF)[(prevdf_CF_medCF$Prevalence >= prevalenceThreshold_CF_medCF)]
phyloseq_CF_medCF_2 = prune_taxa(keepTaxa_CF_medCF, phyloseq_CF_medCF)
##########################################################################################################
#significance testing at the genera level 
#merge ASV_physeq_2 at the genus level 
phyloseq_CF_medCF_3 <- tax_glom(phyloseq_CF_medCF_2, taxrank="Genus")

#convert to deseq object
#This will model using non-continuous Medium, but will control for patient
#Note to self: this is right. If you don't believe it, check the bioconductor vingette 'analyzing RNA-seq data with DESeq2', 'Multi-factor designs' section
ASV_deseq_genus_CF_medCF <- phyloseq_to_deseq2(phyloseq_CF_medCF_3, ~Patient + ch_Day_of_passage)

#workaround to deal with 0s
cts_CF_medCF <- counts(ASV_deseq_genus_CF_medCF)
geoMeans_CF_medCF <- apply(cts_CF_medCF, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
ASV_deseq_genus_CF_medCF <- estimateSizeFactors(ASV_deseq_genus_CF_medCF, geoMeans=geoMeans_CF_medCF)

#deseq standard analysis
ASV_deseq_genus_CF_medCF <- DESeq(ASV_deseq_genus_CF_medCF)

##########################################################################################################

# pulling out our results table, we specify the object, the p-value we are going to use to filter our results 
#and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_Day_CF_medCF <- results(ASV_deseq_genus_CF_medCF, alpha=0.05, contrast=c("ch_Day_of_passage","Day_5", "Day_0"))

# we can get a glimpse at what this table currently holds with the summary command
summary(deseq_res_Day_CF_medCF)

# stitch together with the ASV's taxonomic annotations for a quick look at both together
deseq_res_Day_tax_CF_medCF <- cbind(as(deseq_res_Day_CF_medCF, "data.frame"), as(tax_table(phyloseq_CF_medCF_2)[row.names(deseq_res_Day_CF_medCF), ], "matrix"))

# and sort that table by the baseMean column
sum_Day_CF_medCF <- deseq_res_Day_tax_CF_medCF[order(deseq_res_Day_tax_CF_medCF$baseMean, decreasing=T), ]

write.csv(sum_Day_CF_medCF, "sum_Day_genus_5%p_CF_medCF.csv")

# subset this table to only include these that pass our specified significance level
sum_Day_CF_medCF_sig <- sum_Day_CF_medCF[which(sum_Day_CF_medCF$padj < 0.05), ]
write.csv(sum_Day_CF_medCF_sig, "sum_Day_5%p_0.05pval_CF_medCF_genus_all.csv")

#Making graphs
#DESeq2 visualization methods from phyloseq tutorial
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

#read in manually annotated csv
CF_medCF_day_plot <- read.csv("sum_Day_5%p_0.05pval_CF_medCF_genus_all.csv", row.names = 1, header=TRUE)

#Assign column names as the titles in row 1
#colnames(nonCF_MiPro_day_plot) <- nonCF_MiPro_day_plot[1, ]  
#Drop row 1
#medCFvMi_plot <- medCFvMi_plot[-1,]


# Genus order nonCF MiPro Day
x_CF_medCF = tapply(CF_medCF_day_plot$log2FoldChange, CF_medCF_day_plot$Genus, function(x_CF_medCF) max(x_CF_medCF))
x_CF_medCF = sort(x_CF_medCF, TRUE)
CF_medCF_day_plot$Genus = factor(as.character(CF_medCF_day_plot$Genus), levels=names(x_CF_medCF))

gen_CF_medCF <- ggplot(CF_medCF_day_plot, aes(x=log2FoldChange, y=Genus))+ #color=CFChange)) + 
  geom_point(size=5) + 
  scale_color_viridis(discrete=TRUE) +
  labs(x="Log2 Fold Change", y="Genus")) +
  geom_vline(xintercept = 0) +
  theme(axis.text=element_text(size=12), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.title=element_text(size=12,face="bold"), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=12, face="bold")) +
 geom_label(label="Decreased in medCF", x=-10, y=10, label.size = 0.4, color = "black",  size =3.5)
 #geom_label(label="Increased in medCF", x= 4, y=10, label.size = 0.4, color = "black", size =3.5)


#Figure S8D         
gen_CF_medCF

```
