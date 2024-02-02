# The effect of diet on melanin pigmentation in Hyles lineata
#Authors: Sarah Britton and Goggy Davidowitz

# Libraries
library(ggplot2) #for making plots
library(car) #for Anova function
library(multcomp) #for Tukey tests
library(extrafont) #for plot fonts
library(dplyr) #for piping and grouping
library(ggpubr) #for combining ggplots
library(emmeans) #for post hoc tests
library(lmerTest) # for extracting p vals from mixed models
library(lme4) # mixed models

#Read in Data
#convert variables and combine data sets as needed
diet_data <- read.csv(file="data/Diet_Data.csv")
diet_data$Diet_Treatment <- as.factor(diet_data$Diet_Treatment) #convert to factor

extraction_data<-read.csv(file="data/Extraction_data.csv")
extraction_data_full<- merge(diet_data, extraction_data, by=c("ID", "Diet_Treatment")) 

extraction_data_full$Extraction_Run <- as.factor(extraction_data_full$Extraction_Run) #convert to factor

extraction_data_full$Concentration_PhotoMass<-extraction_data_full$Total_Extracted/extraction_data_full$Photo_Mass
extraction_data_full$Concentration_BlockingMass<-extraction_data_full$Total_Extracted/extraction_data_full$Blocking_Mass

feeding_data <- read.csv(file="data/Feeding_Data.csv")
feeding_data$Diet_Treatment <- as.factor(feeding_data$Diet_Treatment) #convert to factor

feeding_data<- merge(diet_data, feeding_data, by=c("ID", "Diet_Treatment")) 


#Check distribution of  response variables
hist(diet_data$Melanin_Percent)
hist(diet_data$Darkness)
hist(extraction_data_full$Concentration_PhotoMass)


#These lists will be used for making plots for ImageJ data
diet_colors = c("orange", "orange3", "green2", "green4", "yellow") #assigns colors to each diet treatment for figures
diet_order = c('LOW', 'MEDIUM','HIGH','LOW-2', 'HIGH-2') #puts diet in an order that makes sense for figures

#These lists will be used for making plots for extraction data
diet_order_extract = c('LOW','HIGH') #puts diet in an order that makes sense for figures
diet_colors_extract = c("orange", "green2") #assigns colors to each diet treatment for figures

#subset data into diet treatments for within treatment analysis
diet_data_low <- subset(diet_data, Diet_Treatment=='LOW')
diet_data_high <- subset(diet_data, Diet_Treatment=='HIGH')
diet_data_med <- subset(diet_data, Diet_Treatment=='MEDIUM')
diet_data_low2 <- subset(diet_data, Diet_Treatment=='LOW-2')
diet_data_high2 <- subset(diet_data, Diet_Treatment=='HIGH-2')

extraction_data_low <- subset(extraction_data_full, Diet_Treatment=='LOW')
extraction_data_high <- subset(extraction_data_full, Diet_Treatment=='HIGH')


#### Image J Data #####

#Statistics
#Create table of descriptive stats (mean and sd) for all 4 treatment groups
diet_data %>%
  group_by(Diet_Treatment) %>%
  summarize(area_mean = mean(Melanin_Percent, na.rm = TRUE),
            area_sd = sd(Melanin_Percent, na.rm = TRUE),
            area_median = median(Melanin_Percent, na.rm = TRUE),
            darkness_mean = mean(Darkness, na.rm = TRUE),
            darkness_sd = sd(Darkness, na.rm = TRUE),
            darkness_median = median(Darkness, na.rm = TRUE),
            size_mean = mean(Photo_Mass, na.rm = TRUE),
            size_sd = sd(Photo_Mass, na.rm = TRUE),
            size_median = median(Photo_Mass, na.rm = TRUE)) %>%
  as.data.frame

#Melanin Percent
Percent_Melanin_mod <- lm(Melanin_Percent ~ Diet_Treatment, data=diet_data)
summary(Percent_Melanin_mod)
emmeans(Percent_Melanin_mod, specs="Diet_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Percent_Melanin_mod)) #normality of residuals
plot(fitted(Percent_Melanin_mod), residuals(Percent_Melanin_mod)) #homoskedasticity 


#Darkness
Darkness_mod <- lm(Darkness ~ Diet_Treatment, data=diet_data)
summary(Darkness_mod)
emmeans(Darkness_mod, specs="Diet_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Darkness_mod)) #normality of residuals
plot(fitted(Darkness_mod), residuals(Darkness_mod)) #homoskedasticity 

#Figures 

#Melanin Percent
area_graph<-ggplot(diet_data, aes(x=factor(Diet_Treatment, level = diet_order),  Melanin_Percent)) + #x and y to use throughout, reorders x var
  geom_boxplot(aes(color=Diet_Treatment), size=1) +  #boxplot with treatments as differnt colors, resize boxes
  stat_summary(fun.y="mean", aes(color=Diet_Treatment), shape=18, size=1) + #add mean value points to graph (as diamonds)
  scale_color_manual(values=diet_colors) + #chooses colors for reach diet treatment (from above)
  annotate(geom="text", x=c(1,2,3,4,5), y=c(86,88,92,90,89),label=c("a", "ab", "b", "b", "b")) + #adds sig letter labels
  annotate(geom="text", x=4.5, y=95, label="n=258, p=0.002")+
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Diet Treatment") + ylab("Percent Area") #axis labels

area_graph

#Darkness
darkness_graph<-ggplot(diet_data, aes(x=factor(Diet_Treatment, level = diet_order),  Darkness)) + 
  geom_boxplot(aes(color=Diet_Treatment), size=1) + 
  stat_summary(fun.y="mean", aes(color=Diet_Treatment), shape=18, size=1) +
  scale_color_manual(values=diet_colors) + 
  annotate(geom="text", x=c(1,2,3,4,5), y=c(26,26.8,26,26.7,27.5 ),label=c("a", "ab", "a", "b", "b")) +
  annotate(geom="text", x=4.5, y=28.2, label="n=258, p<0.001")+
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Diet Treatment") + ylab("Darkness")

darkness_graph

#combine ImageJ results
pigmentation_combo_graph<-ggarrange(area_graph, darkness_graph, 
                       font.label = list(family="Times New Roman"), labels=c("A", "B"),
                       ncol = 2, hjust=-6.5, align="hv")
pigmentation_combo_graph
ggsave("pigmentation_combo_plot.jpeg")

#Correlation between melanin measures
ggplot(diet_data, aes(Melanin_Percent, Darkness))+ geom_point(aes(color=Diet_Treatment))+ 
  geom_smooth(aes(color=Diet_Treatment), method=lm, se=FALSE) + 
  scale_color_manual(values=diet_colors) +
  theme_classic(base_size = 20) + xlab("Percent Area") + ylab("Darkness")

#### Size ####

#Statistics 

Size_mod<-lm(Photo_Mass ~ Diet_Treatment, data = diet_data)
summary(Size_mod)
emmeans(Size_mod, specs="Diet_Treatment") |> pairs(adjust="tukey")
Photo_Size_Tukey
qqnorm(residuals(Size_mod)) #normality of residuals
plot(fitted(Size_mod),residuals(Size_mod)) #homoskedasticity

Size_Regress<-lm(Melanin_Percent ~ Photo_Mass, data = diet_data)
summary(Size_Regress)

Size_Regress_2<-lm(Darkness ~ Photo_Mass, data = diet_data)
summary(Size_Regress_2)

#Split up by treatment
Size_Regress_low <-lm(Melanin_Percent~ Photo_Mass, data=diet_data_low)
summary(Size_Regress_low)

Size_Regress_high <-lm(Melanin_Percent ~ Photo_Mass, data=diet_data_high)
summary(Size_Regress_high)

Size_Regress_low2 <-lm(Melanin_Percent ~ Photo_Mass, data=diet_data_low2)
summary(Size_Regress_low2)

Size_Regress_high2 <-lm(Melanin_Percent ~ Photo_Mass, data=diet_data_high2)
summary(Size_Regress_high2)

Size_Regress_med <-lm(Melanin_Percent~ Photo_Mass, data=diet_data_med)
summary(Size_Regress_med)

#Figures
size_graph<-ggplot(diet_data, aes(x=factor(Diet_Treatment, level = diet_order),  Photo_Mass)) +
  geom_boxplot(aes(color=Diet_Treatment), size=1)+
  stat_summary(fun.y="mean", aes(color=Diet_Treatment), shape=18, size=1) +
  scale_color_manual(values=diet_colors) +  
  annotate(geom="text", x=c(1,2,3,4,5), y=c(4.6,3.3,5.7,3.5,4),label=c("a", "c", "a", "bc", "ab")) + 
  annotate(geom="text", x=4.5, y=5.3, label="n=258, p<0.001")+
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Diet Treatment") + ylab("Mass (g)")

size_graph

ggplot(diet_data, aes(x=factor(Diet_Treatment, level = diet_order),  Devo_Total_Time)) +
  geom_boxplot(aes(color=Diet_Treatment), size=1)+
  stat_summary(fun.y="mean", aes(color=Diet_Treatment), shape=18, size=1) +
  scale_color_manual(values=diet_colors) +  
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Diet Treatment") + ylab("Development Time (days)")

size_regress<-ggplot(diet_data, aes(Photo_Mass, Melanin_Percent))+ geom_point()+ 
  geom_smooth(method=lm, se=FALSE) + 
  annotate(geom="text", x=4.5, y=92, label="r"^2~"=0.04 , p=0.005")+
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Mass (g)") + ylab("Percent Area")

size_regress

size_regress_2<-ggplot(diet_data, aes(Photo_Mass, Darkness))+ geom_point()+ 
  geom_smooth(method=lm, se=FALSE) + 
  annotate(geom="text", x=4.5, y=27, label="r"^2~"=0.09 , p<0.001")+
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Mass (g)") + ylab("Darkness")

size_regress_2

size_combo_graph<-ggarrange(size_regress, size_regress_2, 
                                    font.label = list(family="Times New Roman"), labels=c("A", "B"),
                                    ncol = 2, hjust=-6.5, align="hv")
size_combo_graph
ggsave("pigmentation_combo_plot.jpeg")

#### Feeding ####

#Statistics

feeding_data %>%
  group_by(Diet_Treatment) %>%
  summarize(fourth_mean = mean(fourth_food_change, na.rm = TRUE),
            fourth_sd = sd(fourth_food_change, na.rm = TRUE),
            fourth_median = median(fourth_food_change, na.rm = TRUE),
            fifth_mean = mean(fifth_food_change, na.rm = TRUE),
            fifth_sd = sd(fifth_food_change, na.rm = TRUE),
            fifth_median = median(fifth_food_change, na.rm = TRUE)) %>%
  as.data.frame

Feeding_Test_4<-lm(fourth_food_change ~ Diet_Treatment, data = feeding_data)
summary(Feeding_Test_4)
emmeans(Feeding_Test_4, specs="Diet_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Feeding_Test_4))
plot(fitted(Feeding_Test_4),residuals(Feeding_Test_4)) #homoskedasticity

Feeding_Test_5<-lm(fifth_food_change ~ Diet_Treatment, data = feeding_data)
summary(Feeding_Test_5)
emmeans(Feeding_Test_5, specs="Diet_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Feeding_Test_5))
plot(fitted(Feeding_Test_5),residuals(Feeding_Test_5)) #homoskedasticity

feeding_graph<-ggplot(feeding_data, aes(x=factor(Diet_Treatment, level = diet_order))) +
  geom_boxplot(aes(y=fourth_food_change, color=Diet_Treatment), size=1) + 
  stat_summary(fun.y="mean", aes(y=fourth_food_change, color=Diet_Treatment), shape=18, size=1) +
  geom_boxplot(aes(y=fifth_food_change, color=Diet_Treatment), size=1) +
  stat_summary(fun.y="mean", aes(y=fifth_food_change, color=Diet_Treatment), shape=18, size=1) +
  scale_color_manual(values=diet_colors) + 
  annotate(geom="text", x=c(1,2,3,4,5), y=c(1750,1380,2050,1220,1190),label=c("a", "b", "a", "b", "b")) +
  annotate(geom="text", x=c(1,2,3,4,5), y=c(500, 430, 650, 400, 490),label=c("cd", "e", "c", "e", "de")) +
  annotate(geom="text", x=4.5, y=2000,label="5th: n=108, p<0.001") +
  annotate(geom="text", x=4.5, y=1800,label="4th: n=108, p<0.001") +
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Diet Treatment") + ylab("Change in food (mg)")

feeding_graph


#Make combo figures
ggarrange(size_graph, feeding_graph, 
          font.label = list(family="Times New Roman"),labels=c("A", "B"), 
          ncol = 2, hjust=-6.5, align="hv")

ggarrange(size_graph, feeding_graph, size_regress, size_regress_2,
          labels=c("A", "B", "C", "D"), ncol = 2, nrow=2)



#### Extraction ####

#Statistics- Extraction Results
t.test(extraction_data_full$Concentration_PhotoMass~extraction_data_full$Diet_Treatment)


Conc_Percent <-lm(Melanin_Percent ~ Concentration_PhotoMass+Diet_Treatment+Concentration_PhotoMass:Diet_Treatment, data=extraction_data_full)
summary(Conc_Percent)

Conc_Percent_low <-lm(Melanin_Percent ~ Concentration_PhotoMass, data=extraction_data_low)
summary(Conc_Percent_low)

Conc_Percent_high <-lm(Melanin_Percent ~ Concentration_PhotoMass, data=extraction_data_high)
summary(Conc_Percent_high)

Conc_Dark <-lm(Darkness ~ Concentration_PhotoMass+Diet_Treatment+Concentration_PhotoMass:Diet_Treatment, data=extraction_data_full)
summary(Conc_Dark)

Conc_Dark_low <-lm(Darkness ~ Concentration_PhotoMass, data=extraction_data_low)
summary(Conc_Dark_low)

Conc_Dark_high <-lm(Darkness ~ Concentration_PhotoMass, data=extraction_data_high)
summary(Conc_Dark_high)



#Figures- Extraction Results 
ggplot(extraction_data_full, aes(x=factor(Diet_Treatment, level=diet_order_extract), Concentration_PhotoMass)) +
  geom_boxplot (aes(color=Diet_Treatment), size=1)+
  stat_summary(fun.y="mean", aes(color=Diet_Treatment), shape=18, size=1) +
  scale_color_manual(values=diet_colors_extract) +
  geom_signif(comparison=list(c("LOW","HIGH")), map_signif_level = TRUE)+
  annotate(geom="text", x=2.4, y=150, label="t(72)=2.87,")+
  annotate(geom="text", x=2.4, y=140, label="p=0.005")+
  theme_classic(base_size = 16) + theme(legend.position="none",text=element_text(family="Times New Roman")) +
  xlab("Diet Treatment") + ylab("Melanin Concentration (ug melanin/g body mass)")

ggplot(extraction_data_full, aes(Concentration_PhotoMass, Melanin_Percent))+ geom_point(aes(color=Diet_Treatment), size=1)+ 
  geom_smooth(aes(color=Diet_Treatment), method=lm, se=TRUE, fill="gray80") + 
  scale_color_manual(values=diet_colors_extract) +
  theme_classic(base_size = 16) + theme(legend.position="none", text=element_text(family="Times New Roman")) +
  xlab("Melanin Concentration (ug melanin/g body mass)") + ylab("Percent Area Melanin")

ggplot(extraction_data_full, aes(Concentration_PhotoMass, Darkness))+ geom_point(aes(color=Diet_Treatment), size=1)+ 
  geom_smooth(aes(color=Diet_Treatment), method=lm, se=TRUE, fill="gray80") + 
  scale_color_manual(values=diet_colors_extract) +
  theme_classic(base_size = 16) + theme(legend.position="none", text=element_text(family="Times New Roman")) +
  xlab("Melanin Concentration (ug melanin/g body mass)") + ylab("Darkness")


#Citations
packages_in_use <- c(sessionInfo()$basePkgs, names(sessionInfo()$loadedOnly))
the_citations_list <- lapply(X=packages_in_use, FUN=citation)
the_citations_list 