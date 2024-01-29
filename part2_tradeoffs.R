#Melanin pigmentation trade-offs in Hyles lineata
#Authors: Sarah Britton and Goggy Davidowitz

# Libraries
library(ggplot2) #for making plots
library(extrafont) #for plot fonts
library(ggpubr) #for combining ggplots
library(ggthemr) 
library(tidyr) 
library(emmeans) #for post hoc tests
library(lmerTest) # for extracting p vals from mixed models
library(lme4) # mixed models
library(dplyr) #load last to avoid conflicts between packages

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


#Merge data sets
tradeoff_data_immune <- merge(tradeoff_data, immune_data, by="Larval_ID")
head(tradeoff_data_immune)

tradeoff_data_muscle <- merge(tradeoff_data, KOH_data, by="Larval_ID")
head(tradeoff_data_muscle)

#Split into two data sets for each diet treatment, needed below
lowTyr_data <- tradeoff_data_immune |>
  filter(Diet_Treatment == "Low")

highTyr_data <- tradeoff_data_immune |>
  filter(Diet_Treatment == "High")

#### Part 1- Immunity #### 

#Visualize  Data
ggplot(aes(x=Avg_Percent_Melanin), data=tradeoff_data_immune) + geom_histogram() 
qqnorm(tradeoff_data_immune$Avg_Percent_Melanin)
#qq plot looks funky! bimodal? due to no intermediate phenotypes?

ggplot(aes(x=Darkness), data=tradeoff_data_immune) + geom_histogram()
qqnorm(tradeoff_data_immune$Darkness)
#looks great!

ggplot(aes(x=Concentration_log), data=tradeoff_data_immune) + geom_histogram()
qqnorm(tradeoff_data_immune$Concentration)
qqnorm(tradeoff_data_immune$Concentration_log)
#qq plot looked funky! Much better after log transformation

ggplot(aes(x=Total_log), data=tradeoff_data_immune) + geom_histogram()
qqnorm(tradeoff_data_immune$Total_log)
#looks good

ggplot(aes(x=Food_Change), data=tradeoff_data_immune) + geom_histogram()
qqnorm(tradeoff_data_immune$Food_Change)
#looks great!

ggplot(aes(x=Mass), data=tradeoff_data_immune) + geom_histogram()
qqnorm(tradeoff_data_immune$Mass)
#looks great!

#Models
Percent_Melanin_mod <- lm(Avg_Percent_Melanin ~ Full_Treatment + Mass, data=tradeoff_data_immune)
summary(Percent_Melanin_mod)
emmeans(Percent_Melanin_mod, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Percent_Melanin_mod)) #normality of residuals

Darkness_mod <- lm(Darkness ~ Photo_Treatment + Mass, data=tradeoff_data_immune)
summary(Darkness_mod)
emmeans(Darkness_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Darkness_mod)) #normality of residuals

Concentration_mod<-lm(Concentration_log ~ Full_Treatment, data=tradeoff_data_immune)
summary(Concentration_mod)
emmeans(Concentration_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Concentration_mod)) #normality of residuals


Total_mod<-lm(Total_log ~ Full_Treatment, data=tradeoff_data_immune)
summary(Total_mod)
emmeans(Total_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Concentration_mod)) #normality of residuals

Interaction_mod <-lm(Concentration_log ~ Photo_Treatment + Diet_Treatment + Photo_Treatment*Diet_Treatment, data=tradeoff_data_immune)
summary(Interaction_mod)
#No interaction between diet and photo

mass_mod <- lm(Concentration_log ~ Mass, data=tradeoff_data_immune)
summary(mass_mod)
Conc_mod2 <- lm(Concentration_log ~ Avg_Percent_Melanin, data=highTyr_data)
summary(Conc_mod2)

Conc_mod3 <- lm(Concentration_log ~ Avg_Percent_Melanin, data=lowTyr_data)
summary(Conc_mod3)

Food_mod <- lm(Food_Change ~ Full_Treatment, data=tradeoff_data_immune)
summary(Food_mod)
emmeans(Food_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Food_mod)) #normality of residuals

Mass_mod <- lm(Mass ~ Full_Treatment, data=tradeoff_data_immune)
summary(Mass_mod)
emmeans(Mass_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Mass_mod)) #normality of residuals

#Figures

treatment_colors <- c("#D55E00","#0072B2", "#E69F00", "#56B4E9") #assigns colors to each treatment for figures
treatment_colors2 <- c("#D55E00","#E69F00","#0072B2", "#56B4E9")

treatment_order <- c("LongPhotoHighTyr", "ShortPhotoHighTyr","LongPhotoLowTyr", "ShortPhotoLowTyr") #puts diet in an order that makes sense for figures
my_labels <- c("Low Melanin\nHigh Tyr", "High Melanin\nHigh Tyr", "Low Melanin\nLow Tyr", "High Melanin\nLow Tyr")
my_labels2 <- c("Low Melanin\nHigh Tyr", "Low Melanin\nLow Tyr", "High Melanin\nHigh Tyr", "Low Melanin\nLow Tyr", "High Melanin\nLow Tyr")

melanin_plot<-ggplot(tradeoff_data_immune, aes(Full_Treatment, y=Avg_Percent_Melanin)) + #x and y to use throughout
  geom_boxplot(aes(color=Full_Treatment), size=1) +  #boxplot with treatments as different colors, resize boxes
  stat_summary(fun.y="mean", color=treatment_colors2, shape=18, size=1) + #add mean value points to graph (as diamonds)
  scale_color_manual(values=treatment_colors2) + #chooses colors for reach diet treatment (from above)
  annotate(geom="text", x=c(1,2,3,4), y=c(43,44,92,93),label=c("a", "a", "b", "b")) + #adds sig letter labels
  annotate(geom="text", x=3.5, y=25, label="F=615.3, p<0.001")+
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Treatment") + ylab("Percent Area") + #axis labels
  scale_x_discrete(labels= my_labels2)

melanin_plot

ggplot(tradeoff_data_immune, aes(x=Photo_Treatment, y=Avg_Percent_Melanin)) + geom_boxplot(size=0.5, width=0.6) +
  stat_summary(fun.y="mean", shape=18, size=1) +
  annotate(geom="text", x=c(1,2), y=c(44, 93),label=c("a", "b")) + 
  annotate(geom="text", x=1, y=80, label="F(2,158)=1,220, p<0.001") +
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Photoperiod") + ylab("Percent Area")

ggplot(tradeoff_data_immune, aes(x=Photo_Treatment, y=Darkness)) + geom_boxplot(size=0.5, width=0.6) +
  stat_summary(fun.y="mean", shape=18, size=1) +
  annotate(geom="text", x=c(1,2), y=c(19, 20),label=c("a", "b")) + 
  annotate(geom="text", x=1, y=23, label="F(2,158)=70.79, p<0.001") +
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Photoperiod") + ylab("Darkness")



ggplot(tradeoff_data_immune, aes(x=Photo_Treatment, Avg_Percent_Melanin)) + #x and y to use throughout
  geom_boxplot(aes(color=Photo_Treatment), size=1) +  #boxplot with treatments as different colors, resize boxes
  stat_summary(fun.y="mean", color=Photo_Treatment, shape=18, size=1) + #add mean value points to graph (as diamonds)
  scale_color_manual(values=treatment_colors) + #chooses colors for reach diet treatment (from above)
  annotate(geom="text", x=c(1,2,3,4), y=c(43,44,92,93),label=c("a", "a", "b", "b")) + #adds sig letter labels
  annotate(geom="text", x=3.5, y=25, label="F=615.3, p<0.001")+
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Photoperiod") + ylab("Percent Area") + #axis labels
  #scale_x_discrete(labels= my_labels)
  
  
  melanin_plot

darkness_plot<-ggplot(tradeoff_data_immune, aes(x=Full_Treatment, y=Darkness)) + #x and y to use throughout
  geom_boxplot(aes(color=Full_Treatment), size=1) +  #boxplot with treatments as different colors, resize boxes
  stat_summary(fun.y="mean", color=treatment_colors, shape=18, size=1) + #add mean value points to graph (as diamonds)
  scale_color_manual(values=treatment_colors) + #chooses colors for reach diet treatment (from above)
  annotate(geom="text", x=c(1,2,3,4), y=c(19,19,20,20),label=c("a", "a", "b", "b")) + #adds sig letter labels
  annotate(geom="text", x=3.5, y=23, label="F=35.39.3, p<0.001")+
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Treatment") + ylab("Darkness") + #axis labels
  scale_x_discrete(labels= my_labels)

darkness_plot

conc_plot<-ggplot(tradeoff_data_immune, aes(x=factor(Full_Treatment, level = treatment_order), y=Concentration_log)) + #x and y to use throughout
  geom_boxplot(aes(color=Full_Treatment), size=0.7, width=0.8) +  #boxplot with treatments as different colors, resize boxes
  stat_summary(fun.y="mean", color=treatment_colors, shape=18, size=1) + #add mean value points to graph (as diamonds)
  scale_color_manual(values=treatment_colors2) + #chooses colors for reach diet treatment (from above)
  annotate(geom="text", x=c(1,2,3,4), y=c(2.6,2.3,2.2,1.9),label=c("a", "b", "bc", "c")) + #adds sig letter labels
  annotate(geom="text", x=3.5, y=2.7, label="F(3,138)=15.58, p<0.001")+
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Treatment") + ylab("log Melanin Concentration (ug/g wire)") + 
  scale_x_discrete(labels= my_labels)


conc_plot


total_plot<-ggplot(tradeoff_data_immune, aes(x=factor(Full_Treatment, level = treatment_order), y=Total_log)) + #x and y to use throughout
  geom_boxplot(aes(color=Full_Treatment), size=0.7, width=0.8) +  #boxplot with treatments as different colors, resize boxes
  stat_summary(fun.y="mean", color=treatment_colors, shape=18, size=1) + #add mean value points to graph (as diamonds)
  scale_color_manual(values=treatment_colors2) + #chooses colors for reach diet treatment (from above)
  annotate(geom="text", x=c(1,2,3,4), y=c(2,1.8,1.7,1.3),label=c("a", "b", "bc", "c")) + #adds sig letter labels
  annotate(geom="text", x=3.5, y=2.5, label="F(3,138)=16.18, p<0.001")+
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Treatment") + ylab("Total melanin (ug)") + 
  scale_x_discrete(labels= my_labels)


total_plot

regression_plot<-ggplot(tradeoff_data_immune, aes(x=Avg_Percent_Melanin, y=Concentration_log, color=Full_Treatment)) +
  geom_point()+ geom_smooth(method=lm, se=FALSE) + 
  scale_color_manual(values=treatment_colors) +
  #annotate(geom="text", x=4.5, y=92, label="r"^2~"=0.11 , p=0.007")+
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Percent Melanin") + ylab("log Melanin Concentration (ug/g wire)") 


regression_plot


food_plot<-ggplot(tradeoff_data_immune, aes(x=Full_Treatment, y=Food_Change)) + #x and y to use throughout
  geom_boxplot(aes(color=Full_Treatment), size=1) +  #boxplot with treatments as different colors, resize boxes
  stat_summary(fun.y="mean", color=treatment_colors2, shape=18, size=1)+ #add mean value points to graph (as diamonds)
  scale_color_manual(values=treatment_colors2) + #chooses colors for reach diet treatment (from above)
  annotate(geom="text", x=c(1,2,3,4), y=c(6.2,6.9,5.3,5.8),label=c("ab", "a", "b", "b")) + #adds sig letter labels
  annotate(geom="text", x=3.5, y=6.5, label="F=4.49, p=0.005")+
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Treatment") + ylab("Consumption") + #axis labels
  scale_x_discrete(labels= my_labels2)

food_plot

mass_plot<-ggplot(tradeoff_data_immune, aes(x=Full_Treatment, y=Mass)) + #x and y to use throughout
  geom_boxplot(aes(color=Full_Treatment), size=1) +  #boxplot with treatments as different colors, resize boxes
  stat_summary(fun.y="mean", color=treatment_colors2, shape=18, size=1) + #add mean value points to graph (as diamonds)
  scale_color_manual(values=treatment_colors2) + #chooses colors for reach diet treatment (from above)
  annotate(geom="text", x=c(1,2,3,4), y=c(6.7,6.6,5.7,5.2),label=c("a", "a", "b", "b")) + #adds sig letter labels
  annotate(geom="text", x=3.5, y=6.5, label="F=17.58, p<0.001")+
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Treatment") + ylab("Mass") + #axis labels
  scale_x_discrete(labels= my_labels2)

mass_plot


ggarrange(melanin_plot, darkness_plot, 
          font.label = list(family="Times New Roman"),labels=c("A", "B"), 
          ncol = 2, hjust=-6.5, align="hv")

ggarrange(food_plot, mass_plot, 
          font.label = list(family="Times New Roman"),labels=c("A", "B"), 
          ncol = 2, hjust=-6.5, align="hv")

boxplot(Avg_Percent_Melanin ~ Photo_Treatment, data=tradeoff_data_full)

ggplot(highTyr_data, aes(x=Avg_Percent_Melanin, y=Concentration_log)) +
  geom_point(aes(color=Photo_Treatment)) + geom_smooth(method=lm, color="black", se=FALSE) + 
  scale_color_manual(values=c("#D55E00","#0072B2")) +
  annotate(geom="text", x=50, y=0.5, label="F(1,67) = 6.30, p = 0.01")+
  annotate(geom="text", x=50, y=0.25, label="y = -0.003x + 1.6")+
  ylim(0, 2.4)+
  theme_classic(base_size = 18) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Percent Melanin") + ylab("log Melanin Concentration (ug/g wire)") 

ggplot(lowTyr_data, aes(x=Avg_Percent_Melanin, y=Concentration_log)) +
  geom_point(aes(color=Photo_Treatment)) + geom_smooth(method=lm, color="black", linetype=2, se=FALSE) + 
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  annotate(geom="text", x=50, y=0.5, label="F(1,71) = 0.79, p = 0.38")+
  annotate(geom="text", x=50, y=0.25, label="y = -0.0007x + 1.3")+
  ylim(0,2.4) + 
  theme_classic(base_size = 18) + theme(legend.position="none", text=element_text(family="Times New Roman")) +
  xlab("Percent Melanin") + ylab("log Melanin Concentration (ug/g wire)") 

mod1 <- lm(Concentration_log ~ Avg_Percent_Melanin, data=highTyr_data)
summary(mod1)

mod2 <- lm(Concentration_log ~ Avg_Percent_Melanin, data=lowTyr_data)
summary(mod2)

#### Part 2- Muscle Mass ####

#Visualize  Data
ggplot(aes(x=Avg_Percent_Melanin), data=tradeoff_data_muscle) + geom_histogram() 
qqnorm(tradeoff_data_muscle$Avg_Percent_Melanin)
#qq plot looks funky! bimodal? due to no intermediate phenotypes?


ggplot(aes(x=Muscle_per_full_mass), data=tradeoff_data_muscle) + geom_histogram()
qqnorm(tradeoff_data_muscle$Muscle_per_full_mass)
#looks good!


ggplot(aes(x=Full_Mass), data=tradeoff_data_muscle) + geom_histogram()
qqnorm(tradeoff_data_muscle$Full_Mass)
#looks great!

ggplot(aes(x=Devo_Time), data=tradeoff_data_muscle) + geom_histogram()
qqnorm(tradeoff_data_muscle$Devo_Time)
#looks great! remember to use count data

#Models
Percent_Melanin_mod2 <- lm(Avg_Percent_Melanin ~ Full_Treatment, data=tradeoff_data_muscle)
summary(Percent_Melanin_mod2)
emmeans(Percent_Melanin_mod2, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Percent_Melanin_mod)) #normality of residuals

Muscle_mod <- lm(Muscle_per_full_mass ~ Full_Treatment, data=tradeoff_data_muscle)
summary(Muscle_mod)
emmeans(Muscle_mod, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Muscle_mod)) #normality of residuals

Mass_mod2 <- lm(Full_Mass ~ Full_Treatment, data=tradeoff_data_muscle)
summary(Mass_mod2)
emmeans(Mass_mod2, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Mass_mod2)) #normality of residuals

Devo_mod <- glm(Devo_Time ~ Full_Treatment, data=tradeoff_data_muscle, family="quasipoisson")
summary(Devo_mod)
emmeans(Devo_mod, specs="Full_Treatment") |> pairs(adjust="tukey")


