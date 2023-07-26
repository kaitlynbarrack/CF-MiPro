---
title: "CFU_plots"
author: "Kaitlyn Barrack"
date: '2023-07-26'
output: html_document
---

```{r}
install.packages("readr")
library(readr)

CFU_metadata <- read_csv("CFMeta_CFU.csv") 
CFU_metadata <- CFU_metadata %>% remove_rownames %>% column_to_rownames(var="Sample_Name")

CFU_metadata.df <- as.data.frame(CFU_metadata)
colnames(CFU_metadata.df)

#Subset by Day 5 only, all genotypes

CFU_metadata.df_D5 <- subset(CFU_metadata.df, Day_of_passage == '5')

#Subset by Day 0 MiPro only, all genotypes

CFU_metadata.df_D0 <- subset(CFU_metadata.df, Day_of_passage == '0')
CFU_metadata.df_MiProD0 <- subset(CFU_metadata.df_D0, Medium == 'MiPro')

#Subset by genotype

CFU_metadata.df_CF <- subset(CFU_metadata.df, Genotype == 'CF')
CFU_metadata.df_nonCF <- subset(CFU_metadata.df, Genotype == 'nonCF')

#Subset genotypes by Day 5 only

CFU_metadata.df_CF5 <- subset(CFU_metadata.df_CF, Day_of_passage == '5')
CFU_metadata.df_nonCF5 <- subset(CFU_metadata.df_nonCF, Day_of_passage == '5')

#Subset by Day 1-5, both genotypes 
CFU_metadata.df15 <- subset(CFU_metadata.df, Day_of_passage %in% c('1','2','3','4','5'))
CFU_metadata.df_nonCF15 <- subset(CFU_metadata.df_nonCF, Day_of_passage %in% c('1','2','3','4','5'))
CFU_metadata.df_CF15 <- subset(CFU_metadata.df_CF, Day_of_passage %in% c('1','2','3','4','5'))

```

Plot CFU on y axis over time
```{r}

#organize your groups to be in order
CFU_metadata.df15$Genotype <- factor(CFU_metadata.df15$Genotype, levels=c("CF", "nonCF"))
CFU_metadata.df15$Medium <- factor(CFU_metadata.df15$Medium, levels=c("MiPro", "low_CF", "median_CF"))

CFU_metadata.df_CF15$Medium <- factor(CFU_metadata.df_CF15$Medium, levels=c("MiPro", "low_CF", "median_CF"))
CFU_metadata.df_nonCF15$Medium <- factor(CFU_metadata.df_nonCF15$Medium, levels=c("MiPro", "low_CF", "median_CF"))

CFU_metadata.df_D5$Genotype <- factor(CFU_metadata.df_D5$Genotype, levels=c("CF", "nonCF"))
CFU_metadata.df_D5$Medium <- factor(CFU_metadata.df_D5$Medium, levels=c("MiPro", "low_CF", "median_CF"))

CFU_metadata.df_MiProD0$Genotype <- factor(CFU_metadata.df_MiProD0$Genotype, levels=c("CF", "nonCF"))

library(ggplot2)

#Figure 3ABC
ggplot(CFU_metadata.df15, aes(y = Anaerobic_CFU, x = Day_of_passage, fill = Medium))+
         geom_smooth()+
         facet_grid(~Genotype)+
  scale_fill_viridis_d()+
  labs(y="log10(Anaerobic CFU/mL)", x="Day of passage")+
   theme_linedraw()+
   theme(strip.text.y = element_text(size = 15))

ggplot(CFU_metadata.df15, aes(y = Aerobic_CFU, x = Day_of_passage, fill = Medium))+
         geom_smooth()+
         facet_grid(~Genotype)+
  scale_fill_viridis_d()+
  labs(y="log10(Aerobic CFU/mL)", x="Day of passage")+
   theme_linedraw()+
   theme(strip.text.y = element_text(size = 15))

ggplot(CFU_metadata.df15, aes(y = Total_CFU, x = Day_of_passage, fill = Medium))+
         geom_smooth()+
         facet_grid(~Genotype)+
  scale_fill_viridis_d()+
  labs(y="log10(Total CFU/mL)", x="Day of passage")+
   theme_linedraw()+
   theme(strip.text.y = element_text(size = 15))

#Figure 3DEF

ggplot(CFU_metadata.df_D5, aes(y = Total_CFU, x = Day_of_passage, fill = Medium))+
         geom_smooth()+
         facet_grid(~Genotype)+
  scale_fill_viridis_d()+
  labs(y="log10(Total CFU/mL)", x="Day of passage")+
   theme_linedraw()+
   theme(strip.text.y = element_text(size = 15))

ggplot(CFU_metadata.df_D5, aes(y = Aerobic_CFU, x = Day_of_passage, fill = Medium))+
         geom_smooth()+
         facet_grid(~Genotype)+
  scale_fill_viridis_d()+
  labs(y="log10(Aerobic CFU/mL)", x="Day of passage")+
   theme_linedraw()+
   theme(strip.text.y = element_text(size = 15))

ggplot(CFU_metadata.df_D5, aes(y = Anaerobic_CFU, x = Day_of_passage, fill = Medium))+
         geom_smooth()+
         facet_grid(~Genotype)+
  scale_fill_viridis_d()+
  labs(y="log10(Anaerobic CFU/mL)", x="Day of passage")+
   theme_linedraw()+
   theme(strip.text.y = element_text(size = 15))


#Figure S1
ggplot(CFU_metadata.df_MiProD0, aes(y = Anaerobic_CFU, x = Genotype, color = Genotype))+
         geom_boxplot()+
         geom_jitter()+
         facet_grid()+
         scale_fill_viridis_d()+
  labs(y="log10(Anaerobic CFU/mL)", x="Genotype")+
  theme(strip.text.y = element_text(size = 15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.text=element_text(size=15))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=12, face="bold"))+
  #theme(legend.title = element_text(size=15))+
  #theme(legend.text = element_text(size=10))+
  theme(axis.text.x=element_blank())+ 
  scale_color_viridis_d()+
   theme_linedraw()+
   theme(strip.text.y = element_text(size = 15))

ggplot(CFU_metadata.df_MiProD0, aes(y = Aerobic_CFU, x = Genotype, color = Genotype))+
         geom_boxplot()+
         geom_jitter()+
         facet_grid()+
         scale_fill_viridis_d()+
  labs(y="log10(Aerobic CFU/mL)", x="Genotype")+
  theme(strip.text.y = element_text(size = 15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.text=element_text(size=15))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=12, face="bold"))+
  #theme(legend.title = element_text(size=15))+
  #theme(legend.text = element_text(size=10))+
  theme(axis.text.x=element_blank())+ 
  scale_color_viridis_d()+
   theme_linedraw()+
   theme(strip.text.y = element_text(size = 15))


ggplot(CFU_metadata.df_MiProD0, aes(y = Total_CFU, x = Genotype, color = Genotype))+
         geom_boxplot()+
         geom_jitter()+
         facet_grid()+
         scale_fill_viridis_d()+
  labs(y="log10(Aerobic+Anaerobic CFU/mL)", x="Genotype")+
  theme(strip.text.y = element_text(size = 15))+
  theme(strip.text.x = element_text(size = 15))+
  theme(axis.text=element_text(size=15))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.title=element_text(size=12, face="bold"))+
  #theme(legend.title = element_text(size=15))+
  #theme(legend.text = element_text(size=10))+
  theme(axis.text.x=element_blank())+ 
  scale_color_viridis_d()+
   theme_linedraw()+
   theme(strip.text.y = element_text(size = 15))


```

```{r}
#Linear modeling

summary(lm(Total_CFU ~ Medium + Patient + Genotype + Day_of_passage, data= CFU_metadata.df15))
summary(lm(Aerobic_CFU ~ Medium + Patient + Genotype + Day_of_passage, data= CFU_metadata.df15))
summary(lm(Anaerobic_CFU ~ Medium + Patient + Genotype + Day_of_passage, data= CFU_metadata.df15))

#by genotype:CF - used for Fig 3
summary(lm(Total_CFU ~ Medium + Patient + Day_of_passage, data= CFU_metadata.df_CF15))
summary(lm(Aerobic_CFU ~ Medium + Patient + Day_of_passage, data= CFU_metadata.df_CF15))
summary(lm(Anaerobic_CFU ~ Medium + Patient + Day_of_passage, data= CFU_metadata.df_CF15))

summary(lmer(Total_CFU ~ Medium + Day_of_passage + (1 | Patient), data=CFU_metadata.df_CF15))
summary(lmer(Aerobic_CFU ~ Medium + Day_of_passage + (1 | Patient), data=CFU_metadata.df_CF15))
summary(lmer(Anaerobic_CFU ~ Medium + Day_of_passage + (1 | Patient), data=CFU_metadata.df_CF15))

#by genotype:nonCF - used for Fig 3
summary(lm(Total_CFU ~ Medium + Patient + Day_of_passage, data= CFU_metadata.df_nonCF))
summary(lm(Aerobic_CFU ~ Medium + Patient + Day_of_passage, data= CFU_metadata.df_nonCF))
summary(lm(Anaerobic_CFU ~ Medium + Patient + Day_of_passage, data= CFU_metadata.df_nonCF))

summary(lmer(Total_CFU ~ Medium + Day_of_passage + (1 | Patient), data=CFU_metadata.df_nonCF15))
summary(lmer(Aerobic_CFU ~ Medium + Day_of_passage + (1 | Patient), data=CFU_metadata.df_nonCF15))
summary(lmer(Anaerobic_CFU ~ Medium + Day_of_passage + (1 | Patient), data=CFU_metadata.df_nonCF15))

#by genotype and Day 5: CF
summary(lmer(Total_CFU ~ Medium + (1|Patient), data= CFU_metadata.df_CF5))
summary(lmer(Aerobic_CFU ~ Medium + (1|Patient), data= CFU_metadata.df_CF5))
summary(lmer(Anaerobic_CFU ~ Medium + (1|Patient), data= CFU_metadata.df_CF5))

#by genotype and Day 5: nonCF
summary(lmer(Total_CFU ~ Medium + (1|Patient), data= CFU_metadata.df_nonCF5))
summary(lmer(Aerobic_CFU ~ Medium + (1|Patient), data= CFU_metadata.df_nonCF5))
summary(lmer(Anaerobic_CFU ~ Medium + (1|Patient), data= CFU_metadata.df_nonCF5))

#by Day 0 MiPro (inoculum)
summary(lm(Total_CFU ~ Genotype, data= CFU_metadata.df_MiProD0))
summary(lm(Anaerobic_CFU ~ Genotype, data= CFU_metadata.df_MiProD0))
summary(lm(Aerobic_CFU ~ Genotype, data= CFU_metadata.df_MiProD0))

```

