# The effect of diet on melanin pigmentation in Hyles lineata
#Authors: Sarah Britton and Goggy Davidowitz

# Libraries
library(ggplot2) #for making plots
library(car) #for Anova function
library(multcomp) #for Tukey tests
library(extrafont) #for plot fonts
library(ggpubr) #for combining ggplots
library(emmeans) #for post hoc tests
library(lmerTest) # for extracting p vals from mixed models
library(lme4) # mixed models
library(cowplot)
library(dplyr) #for piping and grouping


#Read in data, convert variables and combine data sets as needed
diet_data <- read.csv(file="data/diet_data.csv")
diet_data$diet_treatment <- as.factor(diet_data$diet_treatment) #convert to factor
extraction_data<-read.csv(file="data/Extraction_data.csv")
extraction_data_full<- merge(diet_data, extraction_data, by=c("ID", "diet_treatment")) 
extraction_data_full$extraction_run <- as.factor(extraction_data_full$extraction_run) #convert to factor
extraction_data_full$concentration_photo_mass<-extraction_data_full$total_extracted/extraction_data_full$photo_mass
extraction_data_full$concentration_blocking_mass<-extraction_data_full$total_extracted/extraction_data_full$blocking_mass

#Check distribution of  response variables
hist(diet_data$avg_percent_noB)
hist(diet_data$darkness_noB)
hist(extraction_data_full$concentration_photo_mass)


#These lists will be used for making plots for ImageJ data
diet_colors = c("orange", "orange3", "green2", "green4", "yellow") #assigns colors to each diet treatment for figures
diet_order = c('LOW', 'MEDIUM', 'HIGH', 'LOW-2', 'HIGH-2') #puts diet in an order that makes sense for figures
diet_labels=c('Low-P/\nLow-T', 'Low-P/\nMed-T', 'Low-P/\nHigh-T', 'High-P/\nLow-T', 'High-P/\nHigh-T')

#These lists will be used for making plots for extraction data
diet_order_extract = c('LOW','HIGH') #puts diet in an order that makes sense for figures
diet_labels_extract=c('Low-P/\nLow-T','Low-P/\nHigh-T')
diet_colors_extract = c("orange", "green2") #assigns colors to each diet treatment for figures

####Statistics####
#Create table of descriptive stats (mean and sd) for all 4 treatment groups
diet_data %>%
  group_by(diet_treatment) %>%
  summarize(area_mean = mean(avg_percent, na.rm = TRUE),
            area_sd = sd(avg_percent, na.rm = TRUE),
            area_median = median(avg_percent, na.rm = TRUE),
            darkness_mean = mean(darkness, na.rm = TRUE),
            darkness_sd = sd(darkness, na.rm = TRUE),
            darkness_median = median(darkness, na.rm = TRUE),
            size_mean = mean(photo_mass, na.rm = TRUE),
            size_sd = sd(photo_mass, na.rm = TRUE),
            size_median = median(photo_mass, na.rm = TRUE),
            fourth_mean = mean(X4th_food_change, na.rm = TRUE),
            fourth_sd = sd(X4th_food_change, na.rm = TRUE),
            fourth_median = median(X4th_food_change, na.rm = TRUE),
            fifth_mean = mean(X5th_food_change, na.rm = TRUE),
            fifth_sd = sd(X5th_food_change, na.rm = TRUE),
            fifth_median = median(X5th_food_change, na.rm = TRUE)) %>%
  as.data.frame

#Melanin Proportion
proportion_mod <- lm(avg_percent ~ diet_treatment, data=diet_data)
summary(proportion_mod)
emmeans(proportion_mod, specs="diet_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(proportion_mod)) #normality of residuals
plot(fitted(proportion_mod), residuals(proportion_mod)) #homoskedasticity 

#Darkness
darkness_mod <- lm(darkness ~ diet_treatment, data=diet_data)
summary(darkness_mod)
emmeans(darkness_mod, specs="diet_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(darkness_mod)) #normality of residuals
plot(fitted(darkness_mod), residuals(darkness_mod)) #homoskedasticity 

#Size
size_mod<-lm(photo_mass ~ diet_treatment, data = diet_data)
summary(size_mod)
emmeans(size_mod, specs="diet_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(size_mod)) #normality of residuals
plot(fitted(size_mod),residuals(size_mod)) #homoskedasticity

size_regress<-lm(avg_percent ~ photo_mass, data = diet_data)
summary(size_regress)

size_regress2<-lm(darkness ~ photo_mass, data = diet_data)
summary(size_regress2)

#Feeding
Feeding_Test_4<-lm(X4th_food_change ~ diet_treatment, data = diet_data)
summary(Feeding_Test_4)
emmeans(Feeding_Test_4, specs="diet_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Feeding_Test_4))
plot(fitted(Feeding_Test_4),residuals(Feeding_Test_4)) #homoskedasticity

Feeding_Test_5<-lm(X5th_food_change ~ diet_treatment, data = diet_data)
summary(Feeding_Test_5)
emmeans(Feeding_Test_5, specs="diet_treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Feeding_Test_5))
plot(fitted(Feeding_Test_5),residuals(Feeding_Test_5)) #homoskedasticity

#Extraction
t.test(extraction_data_full$concentration_photo_mass~extraction_data_full$diet_treatment)

extraction_test <- lm(concentration_photo_mass ~ avg_percent*darkness, data=extraction_data_full)
summary(extraction_test)

####Figures for Paper####
#Figure 2
area_graph_paper<-ggplot(diet_data, aes(x=factor(diet_treatment, level = diet_order),  avg_percent)) + #x and y to use throughout, reorders x var
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.17, alpha=0.5)+
  stat_summary(fun ="mean", shape=18, size=1.2) + #add mean value points to graph (as diamonds)
  annotate(geom="text", x=c(1,2,3,4,5), y=c(87,89,93,91,90),label=c("a", "ab", "b", "b", "b")) + #adds sig letter labels
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + #text font and size, remove legend
  scale_x_discrete(labels = diet_labels) +
  xlab("Diet Treatment") + ylab("Percent melanic area") #axis labels
area_graph_paper

darkness_graph_paper<-ggplot(diet_data, aes(x=factor(diet_treatment, level = diet_order),  darkness)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.17,alpha=0.5) +
  stat_summary(fun="mean", shape=18, size=1.2) +
  annotate(geom="text", x=c(1,2,3,4,5), y=c(26.3,27,26.2,26.7,27.5 ),label=c("a", "ab", "a", "b", "b")) +
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  scale_x_discrete(labels = diet_labels) +
  xlab("Diet Treatment") + ylab("Darkness")
darkness_graph_paper

extraction_paper <- ggplot(extraction_data_full, aes(x=factor(diet_treatment, level=diet_order_extract), concentration_photo_mass)) +
  geom_boxplot (outlier.shape = NA) +
  geom_jitter(width=0.15, alpha=0.5) +
  stat_summary(fun="mean",shape=18, size=1.2) +
  geom_signif(comparison=list(c("LOW","HIGH")), map_signif_level = TRUE, textsize=7, vjust=0.5)+
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  scale_x_discrete(labels = diet_labels_extract) +
  xlab("Diet Treatment") + ylab("Melanin concentration\n(ug melanin/g body mass)")
extraction_paper

ggarrange(area_graph_paper, darkness_graph_paper, extraction_paper, 
          font.label = list(size=14, family="Times New Roman"),labels=c("A", "B", "C"),
          nrow=1,
          hjust=-8, align="hv")

#Figure 3
size_graph_paper<-ggplot(diet_data, aes(x=factor(diet_treatment, level = diet_order),  photo_mass)) +
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.17, alpha=0.5)+
  stat_summary(fun.y="mean",shape=18, size=1.2) +
  scale_color_manual(values=diet_colors) +  
  annotate(geom="text", x=c(1,2,3,4,5), y=c(4.7,3.4,5.8,3.6,4.1),label=c("a", "c", "a", "bc", "ab")) + 
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  scale_x_discrete(labels = diet_labels) +
  xlab("Diet Treatment") + ylab("Mass (g)")
size_graph_paper

feeding_4th_paper <- ggplot(diet_data, aes(x=factor(diet_treatment, level = diet_order), y=X4th_food_change)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.13, alpha=0.5) +
  stat_summary(fun.y="mean", shape=18, size=1.2) +
  annotate(geom="text", x=c(1,2,3,4,5), y=c(430, 400, 590, 370, 410),label=c("ab", "c", "a", "c", "bc")) +
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  scale_x_discrete(labels = diet_labels) +
  xlab("Diet Treatment") + ylab("4th instar: Change in food (mg)")
feeding_4th_paper
  
feeding_5th_paper <-ggplot(diet_data, aes(x=factor(diet_treatment, level = diet_order), y=X5th_food_change)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width=0.13, alpha=0.5) +
  stat_summary(fun.y="mean", shape=18, size=1.2) +
  annotate(geom="text", x=c(1,2,3,4,5), y=c(1760,1380,2050,1220,1190),label=c("a", "b", "a", "b", "b")) +
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  scale_x_discrete(labels = diet_labels) +
  xlab("Diet Treatment") + ylab("5th instar: Change in food (mg)")
feeding_5th_paper

ggarrange(size_graph_paper, feeding_4th_paper, feeding_5th_paper, 
          font.label = list(size=14, family="Times New Roman"),labels=c("A", "B", "C"),
          nrow=1,
          hjust=-7.5, align="hv")

#Figure S1
size_regress<-ggplot(diet_data, aes(photo_mass, avg_percent))+ geom_point(alpha=0.7)+ 
  geom_smooth(method=lm, color="Black") + 
  annotate(geom="text", x=4.5, y=96, label="r"^2~"=0.04 , p=0.002")+
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Larval mass (g)") + ylab("Proportion melanic area (%)")+
  ylim(50, 100) 
size_regress

size_regress_2<-ggplot(diet_data, aes(photo_mass, darkness))+ geom_point(alpha=0.7)+ 
  geom_smooth(method=lm, color="Black") + 
  annotate(geom="text", x=4.5, y=29, label="r"^2~"=0.10, p<0.001")+
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Larval mass (g)") + ylab("Darkness")+
  ylim(15,30) 
size_regress_2 

ggarrange(size_regress, size_regress_2,
          font.label = list(size=14, family="Times New Roman"),labels=c("A", "B"),
          nrow=1,
          hjust=-7.5, align="hv")

####Figures for presentations####

#Melanin Percent
area_graph_slides<-ggplot(diet_data, aes(x=factor(diet_treatment, level = diet_order),  avg_percent)) + 
  geom_boxplot(aes(color=diet_treatment), outlier.shape=NA) +  
  geom_jitter(aes(color=diet_treatment), width=0.17, alpha=0.5)+
  stat_summary(fun.y="mean", aes(color=diet_treatment), shape=18, size=1) + 
  scale_color_manual(values=diet_colors) + 
  annotate(geom="text", x=c(1,2,3,4,5), y=c(86,88,92,90,89),label=c("a", "ab", "b", "b", "b")) + 
  scale_x_discrete(labels = diet_labels) +
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Diet Treatment") + ylab("Percent melanic area") 
area_graph_slides

#Darkness
darkness_graph_slides<-ggplot(diet_data, aes(x=factor(diet_treatment, level = diet_order),  darkness)) + 
  geom_boxplot(aes(color=diet_treatment), outlier.shape=NA) +  
  geom_jitter(aes(color=diet_treatment), width=0.17, alpha=0.5)+  stat_summary(fun.y="mean", aes(color=diet_treatment), shape=18, size=1) +
  scale_color_manual(values=diet_colors) + 
  annotate(geom="text", x=c(1,2,3,4,5), y=c(26,26.8,26,26.7,27.5 ),label=c("a", "ab", "a", "b", "b")) +
  scale_x_discrete(labels = diet_labels) +
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Diet Treatment") + ylab("Darkness")
darkness_graph_slides

#Correlation between melanin measures
ggplot(diet_data, aes(avg_percent, darkness))+ geom_point(aes(color=diet_treatment))+ 
  geom_smooth(aes(color=diet_treatment), method=lm, se=FALSE) + 
  scale_color_manual(values=diet_colors) +
  theme_classic(base_size = 20) + xlab("Percent Area") + ylab("Darkness")

#Size
size_graph_slides<-ggplot(diet_data, aes(x=factor(diet_treatment, level = diet_order),  photo_mass)) +
  geom_boxplot(aes(color=diet_treatment), outlier.shape=NA)+
  geom_jitter(aes(color=diet_treatment), width=0.17, alpha=0.7)+
  stat_summary(fun.y="mean", aes(color=diet_treatment), shape=18, size=1) +
  scale_color_manual(values=diet_colors) +  
  annotate(geom="text", x=c(1,2,3,4,5), y=c(4.6,3.3,5.7,3.5,4),label=c("a", "c", "a", "bc", "ab")) + 
  #annotate(geom="text", x=4.5, y=5.3, label="n=258, p<0.001")+
  scale_x_discrete(labels = diet_labels) +
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Diet Treatment") + ylab("Mass (g)")
size_graph_slides


#Feeding
feeding_graph_slides<-ggplot(diet_data, aes(x=factor(diet_treatment, level = diet_order))) +
  geom_boxplot(aes(y=X4th_food_change, color=diet_treatment), outlier.shape = NA) + 
  geom_jitter(aes(y=X4th_food_change, color=diet_treatment), width=0.17, alpha=0.7)+
  stat_summary(fun.y="mean", aes(y=X4th_food_change, color=diet_treatment), shape=18, size=1) +
  geom_boxplot(aes(y=X5th_food_change, color=diet_treatment), outlier.shape=NA) +
  geom_jitter(aes(y=X5th_food_change, color=diet_treatment), width=0.17, alpha=0.7)+
  stat_summary(fun.y="mean", aes(y=X5th_food_change, color=diet_treatment), shape=18, size=1) +
  scale_color_manual(values=diet_colors) + 
  annotate(geom="text", x=c(1,2,3,4,5), y=c(1750,1380,2050,1220,1190),label=c("a", "b", "a", "b", "b")) +
  annotate(geom="text", x=c(1,2,3,4,5), y=c(500, 430, 650, 400, 490),label=c("cd", "e", "c", "e", "de")) +
  #annotate(geom="text", x=4.5, y=2000,label="5th: n=108, p<0.001") +
  #annotate(geom="text", x=4.5, y=1800,label="4th: n=108, p<0.001") +
  scale_x_discrete(labels = diet_labels) +
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Diet Treatment") + ylab("Change in food (mg)")
feeding_graph_slides

#Extraction
ggplot(extraction_data_full, aes(x=factor(diet_treatment, level=diet_order_extract), concentration_photo_mass)) +
  geom_boxplot (aes(color=diet_treatment), outlier.shape=NA) +
  geom_jitter(aes(color=diet_treatment), alpha=0.8) +
  stat_summary(fun.y="mean", aes(color=diet_treatment), shape=18, size=1) +
  scale_color_manual(values=diet_colors_extract) +
  geom_signif(comparison=list(c("LOW","HIGH")), map_signif_level = TRUE)+
  scale_x_discrete(labels = diet_labels_extract) +
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Diet Treatment") + ylab("Melanin concentration (ug melanin/g mass)")

ggplot(extraction_data_full, aes(avg_percent, concentration_photo_mass))+ geom_point(size=1) + 
  geom_smooth(method=lm, se=TRUE, fill="gray80") + 
  annotate(geom="text", x=60, y=150, label="r"^2~"=0.04 , p=0.03")+
  theme_classic(base_size = 16) + theme(legend.position="none", text=element_text(family="Times New Roman")) +
  ylab("Melanin Concentration (ug melanin/g body mass)") + xlab("Percent Area Melanin")

ggplot(extraction_data_full, aes(concentration_photo_mass, melanin_metric_noB))+ geom_point(size=1)+ 
  geom_smooth(method=lm, se=TRUE, fill="gray80") + 
  #annotate(geom="text", x=20, y=150, label="r"^2~"=0.04 , p=0.03")+
  theme_classic(base_size = 16) + theme(legend.position="none", text=element_text(family="Times New Roman")) +
  ylab("Melanin") + xlab("Concentration")



