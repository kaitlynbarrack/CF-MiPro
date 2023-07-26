---
title: "Figure 2C: Differentially-abundant genera between Media"
author: "Kaitlyn Barrack"
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
