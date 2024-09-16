#Family level plasticity from trade-off experiment
#Authors: Sarah Britton and Goggy Davidowitz

library(cowplot)
library(gridExtra)
library(ggplot2) #for making plots
library(extrafont) #for plot fonts
library(ggpubr) #for combining ggplots
library(tidyr) 
library(emmeans) #for post hoc tests
library(lmerTest) # for extracting p vals from mixed models
library(lme4) # mixed models
library(car)
library(knitr)
library(patchwork)
library(lmerTest) # for extracting p vals from mixed models
library(dplyr) 

#Read in data and convert variables as needed
tradeoff_data <- read.csv("data/tradeoff_data.csv")
tradeoff_data<-tradeoff_data |>
  mutate(Full_Treatment=as.factor(Full_Treatment)) |>
  mutate(Photo_Treatment=as.factor(Photo_Treatment)) |>
  mutate(Diet_Treatment=as.factor(Diet_Treatment)) |>
  mutate(Family=as.factor(Family)) |>
  mutate(Sex=as.factor(Sex))

immune_data<-read.csv("data/immune_data.csv")
immune_data<-immune_data |>
  mutate(Concentration=as.numeric(Concentration)) |>
  mutate(Concentration_log = log10(Concentration)) |>
  mutate(Total_log = log10(Total_Melanin)) 

KOH_data <- read.csv("data/KOH_data.csv")
KOH_data <- KOH_data |>
  mutate(muscle_percentage=muscle_percent*100) |>
  filter(exclude != "Yes")

wing_data <- read.csv("data/wing_data.csv")

#Merge data sets
tradeoff_data_immune <- merge(tradeoff_data, immune_data, by="Larval_ID")
tradeoff_data_adult_intermediate <- merge(tradeoff_data, KOH_data, by="Larval_ID")
tradeoff_data_adult <- merge(tradeoff_data_adult_intermediate, wing_data, by="Larval_ID")


#Split up data
low_diet<-tradeoff_data|>
  filter(Diet_Treatment == 'Low') |>
  filter(!is.na(Family))

high_diet <- tradeoff_data|>
  filter(Diet_Treatment == 'High') |>
  filter(!is.na(Family))

short_photo <- tradeoff_data|>
  filter(Photo_Treatment == 'Short') |>
  filter(!is.na(Family))

long_photo <- tradeoff_data|>
  filter(Photo_Treatment == 'Long') |>
  filter(!is.na(Family))


#colors for plots
diet_colors <- c("#FF7F0E","#1B9E77")

#labels for plots
treatment_order <- c("LongPhotoLowTyr", "LongPhotoHighTyr", "ShortPhotoLowTyr", "ShortPhotoHighTyr") 
my_labels <- c("Long Photo\nLow Tyr", "Long Photo\nHigh Tyr", "Short Photo\nLow Tyr", "Short Photo\nHigh Tyr")

treatment_order2 <- c("LongPhotoHighTyr", "ShortPhotoHighTyr", "LongPhotoLowTyr",  "ShortPhotoLowTyr") 
my_labels2 <- c("Long Photo\nHigh Tyr", "Short Photo\nHigh Tyr", "Long Photo\nLow Tyr", "Short Photo\nLow Tyr")

family_labels<-c("1", "2", "3", "4", "5", "6")

#By diet
high_diet_family <- ggplot(high_diet, aes(x=Photo_Treatment, y=Percent_G)) + 
  stat_summary(aes(group=Family,color=Family), fun=mean, shape="square", size=0.7) +
  stat_summary(aes(group=Family,color=Family), geom="errorbar",  fun.data= "mean_se",  width = 0.1, size = 0.5) + 
  stat_summary(aes(group=Family,color=Family), geom = "line", fun = mean, size=1) +
  scale_color_discrete(labels=family_labels) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman")) + 
  xlab("Photoperiod") + ylab("Percent melanin") + ylim(5,75) + ggtitle("High tyrosine diet") +
  guides(color = guide_legend(nrow = 1))
high_diet_family

low_diet_family <- ggplot(low_diet, aes(x=Photo_Treatment, y=Percent_G)) + 
  stat_summary(aes(group=Family,color=Family), fun=mean, shape="square", size=0.7) +
  stat_summary(aes(group=Family,color=Family), geom="errorbar",  fun.data= "mean_se",  width = 0.1, size = 0.5) + 
  stat_summary(aes(group=Family,color=Family), geom = "line", fun = mean, size=1) +
  scale_color_discrete(labels=family_labels) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman")) + 
  xlab("Photoperiod") + ylab("Percent melanin") + ylim(5,75) + ggtitle("Low tyrosine diet") +
  guides(color = guide_legend(nrow = 1))
low_diet_family

high_diet_mod <- lm(Percent_G ~ Photo_Treatment + Family + Photo_Treatment*Family, data=high_diet)
Anova(high_diet_mod, Type=2)
qqnorm(residuals(high_diet_mod)) 

low_diet_mod <- lm(Percent_G ~ Photo_Treatment + Family + Photo_Treatment*Family, data=low_diet)
Anova(low_diet_mod, Type=2)
qqnorm(residuals(low_diet_mod)) 

ggarrange(high_diet_family, low_diet_family, 
          font.label = list(size=14, family="Times New Roman"),labels=c("A: High Tyr Diet", "B: Low Tyr Diet"),
          ncol=2, hjust=-1, align="hv")

#By photo
short_photo_family <- ggplot(short_photo, aes(x=Diet_Treatment, y=Percent_G)) + 
  stat_summary(aes(group=Family,color=Family), fun=mean, shape="square", size=0.7) +
  stat_summary(aes(group=Family,color=Family), geom="errorbar",  fun.data= "mean_se",  width = 0.1, size = 0.5) + 
  stat_summary(aes(group=Family,color=Family), geom = "line", fun = mean, size=1) +
  scale_color_discrete(labels=family_labels) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman")) + 
  xlab("Diet") + ylab("Percent melanin") +ylim(50,80) + ggtitle("Short photoperiod") +
  guides(color = guide_legend(nrow = 1))
short_photo_family

long_photo_family <- ggplot(long_photo, aes(x=Diet_Treatment, y=Percent_G)) + 
  stat_summary(aes(group=Family,color=Family), fun=mean, shape="square", size=0.7) +
  stat_summary(aes(group=Family,color=Family), geom="errorbar",  fun.data= "mean_se",  width = 0.1, size = 0.5) + 
  stat_summary(aes(group=Family,color=Family), geom = "line", fun = mean, size=1) +
  scale_color_discrete(labels=family_labels) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman")) + 
  xlab("Diet") + ylab("Percent melanin") + ylim(0,15) + ggtitle("Long photoperiod") +
  guides(color = guide_legend(nrow = 1))
long_photo_family


short_photo_mod <- lm(Percent_G ~ Diet_Treatment*Family, data=short_photo)
Anova(short_photo_mod, Type=2)
qqnorm(residuals(high_diet_mod)) 

long_photo_mod <- lm(Percent_G ~ Diet_Treatment*Family, data=long_photo)
Anova(long_photo_mod, Type=2)
qqnorm(residuals(high_diet_mod))

ggarrange(short_photo_family, long_photo_family, 
          font.label = list(size=14, family="Times New Roman"),labels=c("C: Long Photoperiod", "D: Short Photoperiod"),
          ncol=2,hjust=-1, align="hv")


(high_diet_family + low_diet_family)/(short_photo_family + long_photo_family) + 
  plot_layout(guides="collect", axis_titles ="collect") & theme(legend.position = "bottom", legend.box = "horizantal")
+ plot_annotation(tag_levels = "A") 

combined_plot <- (high_diet_family + low_diet_family) /
  (short_photo_family + long_photo_family) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect", axis_titles = "collect") &
  theme(legend.position = "bottom") 

combined_plot

full_family_mod <- lm(Percent_G ~  Photo_Treatment*Family*Diet_Treatment, data=tradeoff_data)
Anova(full_family_mod, Type=2)

full_family_mod2 <- lmer(Percent_G ~  Photo_Treatment + Diet_Treatment + (1|Family) + 
                           (1|Family:Photo_Treatment) + (1|Family:Diet_Treatment) +
                           (1|Family:Photo_Treatment:Diet_Treatment), data=tradeoff_data)
summary(full_family_mod2)
Anova(full_family_mod2)


