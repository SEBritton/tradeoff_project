#Melanin pigmentation trade-offs in Hyles lineata
#Authors: Sarah Britton and Goggy Davidowitz

intstall.packages("parameters")
install.packages("parameters")
install.packages("MuMIn")

# Libraries
library(cowplot)
library(gridExtra)
library(ggplot2) #for making plots
library(extrafont) #for plot fonts
library(ggpubr) #for combining ggplots
library(tidyr) 
library(emmeans) #for post hoc tests
library(lmerTest) # for extracting p vals from mixed models
library(lme4) # mixed models
library(piecewiseSEM)
library(multcompView)
library(car)
library(knitr)
library(patchwork)
library(lmerTest)
library(dplyr) #load last to avoid conflicts between packages
library(MuMIn)


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


#Split up by diet treatment
lowdiet_immune<-tradeoff_data_immune |>
  filter(Diet_Treatment == 'Low')

highdiet_immune <- tradeoff_data_immune |>
  filter(Diet_Treatment == 'High')

lowdiet_adult<-tradeoff_data_adult |>
  filter(Diet_Treatment == 'Low')

highdiet_adult <- tradeoff_data_adult |>
  filter(Diet_Treatment == 'High')

#colors for plots
diet_colors <- c("#FF7F0E","#1B9E77")

#labels for plots
treatment_order <- c("LongPhotoLowTyr", "LongPhotoHighTyr", "ShortPhotoLowTyr", "ShortPhotoHighTyr") 
my_labels <- c("Long Photo\nLow Tyr", "Long Photo\nHigh Tyr", "Short Photo\nLow Tyr", "Short Photo\nHigh Tyr")

treatment_order2 <- c("LongPhotoHighTyr", "ShortPhotoHighTyr", "LongPhotoLowTyr",  "ShortPhotoLowTyr") 
my_labels2 <- c("Long Photo\nHigh Tyr", "Short Photo\nHigh Tyr", "Long Photo\nLow Tyr", "Short Photo\nLow Tyr")

family_labels<-c("1", "2", "3", "4", "5", "6")

####Treatment Effects####

#Visualize data
ggplot(aes(x=Avg_Percent_Melanin), data=tradeoff_data) + geom_histogram() 
ggplot(aes(x=Darkness), data=tradeoff_data) + geom_histogram()
ggplot(aes(x=Mass), data=tradeoff_data_immune) + geom_histogram()
ggplot(aes(x=full_mass), data=tradeoff_data_adult) + geom_histogram()
ggplot(aes(x=Devo_5th), data=tradeoff_data_adult) + geom_histogram()
ggplot(aes(x=Devo_Total), data=tradeoff_data_adult) + geom_histogram(binwidth = 1)

#Statistics
tradeoff_data %>%
  group_by(Photo_Treatment) %>%
  summarize(area_mean = mean(Avg_Percent, na.rm = TRUE),
            area_sd = sd(Avg_Percent, na.rm = TRUE),
            area_median = median(Avg_Percent, na.rm = TRUE),
            darkness_mean = mean(Darkness, na.rm = TRUE),
            darkness_sd = sd(Darkness, na.rm = TRUE),
            darkness_median = median(Darkness, na.rm = TRUE)) %>%
  as.data.frame

tradeoff_data_adult %>%
  group_by(Photo_Treatment) %>%
  summarize(mass_mean = mean(full_mass, na.rm = TRUE),
            mass_sd = sd(full_mass, na.rm = TRUE)) %>%
  as.data.frame

tradeoff_data_immune %>%
  group_by(Photo_Treatment) %>%
  summarize(mass_mean = mean(Mass, na.rm = TRUE),
            mass_sd = sd(Mass, na.rm = TRUE)) %>%
  as.data.frame


#Pigmentation
Percent_Melanin_mod <- lmer(Percent_G ~ Photo_Treatment + Diet_Treatment + Photo_Treatment:Diet_Treatment + (1|Family), data=tradeoff_data)
summary(Percent_Melanin_mod)
emmeans(Percent_Melanin_mod, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Percent_Melanin_mod)) 
plot(Percent_Melanin_mod)

Darkness_mod <- lmer(Darkness ~ Photo_Treatment + Diet_Treatment + Photo_Treatment:Diet_Treatment + (1|Family), data=tradeoff_data)
summary(Darkness_mod)
emmeans(Darkness_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Darkness_mod)) 
plot(Darkness_mod)

#Size
Food_mod<-lmer(Food_Change ~ Photo_Treatment + Diet_Treatment + Photo_Treatment*Diet_Treatment + (1|Family), data=tradeoff_data_immune)
summary(Food_mod)
emmeans(Food_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Food_mod)) 
plot(Food_mod)

var.test(Food_Change ~ Diet_Treatment, data=tradeoff_data_immune)
#no difference in variance

Mass_mod <- lmer(Mass ~ Photo_Treatment + Diet_Treatment + Photo_Treatment*Diet_Treatment + (1|Family), data=tradeoff_data_immune)
summary(Mass_mod)
emmeans(Mass_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Mass_mod)) 
plot(Mass_mod)

var.test(Mass ~ Diet_Treatment, data=tradeoff_data_immune)
#no difference in variance


Adult_mass_mod <- lmer(full_mass ~ Photo_Treatment + Diet_Treatment + Photo_Treatment*Diet_Treatment +(1|Family), data=tradeoff_data_adult)
summary(Adult_mass_mod)
emmeans(Devo_mod, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Adult_mass_mod)) 
plot(Adult_mass_mod)

var.test(full_mass ~ Diet_Treatment, data=tradeoff_data_adult)
#no difference in variance

Devo_mod <- glm(Devo_Total ~ Photo_Treatment + Diet_Treatment + Photo_Treatment*Diet_Treatment, data=tradeoff_data_adult, family="quasipoisson")
summary(Devo_mod)
emmeans(Devo_mod, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Devo_mod)) 
plot(Devo_mod)

melanin_vs_size <- lm(Mass~Darkness, data = tradeoff_data_immune)
summary(melanin_vs_size)
qqnorm(residuals(melanin_vs_size)) 
plot(melanin_vs_size)

#Figures
melanin_plot<-ggplot(tradeoff_data, aes(x=factor(Full_Treatment, level = treatment_order), y=Avg_Percent)) + #x and y to use throughout
  geom_boxplot( size=0.7, outlier.shape=NA) + #boxplot with colors for each treatment
  geom_jitter(width=0.25, alpha=0.5) + #show raw data points
  stat_summary(fun.y="mean", shape=18, size=1.2) + #add mean value points to graph (as diamonds)
  annotate(geom="text", x=c(1,2,3,4), y=c(43,50,94,94),label=c("a", "b", "c", "c")) + #adds sig letter labels
  #scale_color_manual(values=diet_colors)+
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Treatment") + ylab("Percent melanic\narea") + ylim(0,100)+
  scale_x_discrete(labels= my_labels) 
melanin_plot

darkness_plot<-ggplot(tradeoff_data_immune, aes(x=factor(Full_Treatment, level = treatment_order), y=Darkness)) +
  geom_boxplot(size=0.7, outlier.shape=NA) + 
  geom_jitter(width=0.25, alpha=0.5) + #show raw data points
  stat_summary(fun.y="mean", shape=18, size=1.2) + 
  annotate(geom="text", x=c(1,2,3,4), y=c(18.5,20,18.5,20),label=c("a", "a", "b", "b")) + 
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Darkness") + 
  scale_x_discrete(labels= my_labels)
darkness_plot

food_plot<-ggplot(tradeoff_data_immune, aes(x=factor(Full_Treatment, level = treatment_order), y=Food_Change)) + #x and y to use throughout
  geom_boxplot(size=0.7, outlier.shape=NA) +
  geom_jitter(width=0.25, alpha=0.5) + 
  stat_summary(fun.y="mean", shape=18, size=1) + #add mean value points to graph (as diamonds)
  annotate(geom="text", x=c(1,2,3,4), y=c(7, 6.3,5.9,5.4),label=c("a", "ab", "b", "b")) + #adds sig letter labels
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Treatment") + ylab("Food consumption\n(g/ 24hr)") + ylim(0,7.5)+
  scale_x_discrete(labels= my_labels)
food_plot

mass_plot<-ggplot(tradeoff_data_immune, aes(x=factor(Full_Treatment, level = treatment_order), y=Mass)) + 
  geom_boxplot(size=0.7, outlier.shape=NA) +
  geom_jitter(width=0.25, alpha=0.5) +   
  stat_summary(fun.y="mean", shape=18, size=1) +
  annotate(geom="text", x=c(1,2,3,4), y=c(6.5,6.8,5.2, 5.8),label=c("a", "a", "b", "b")) + 
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Larval mass (g)") + ylim(0,7) +
  scale_x_discrete(labels= my_labels)
mass_plot

adult_mass_plot<-ggplot(tradeoff_data_adult, aes(x=factor(Full_Treatment, level = treatment_order), y=full_mass)) + 
  geom_boxplot(size=0.7, outlier.shape=NA) +
  geom_jitter(width=0.25, alpha=0.5) + 
  stat_summary(fun.y="mean", shape=18, size=1) + 
  annotate(geom="text", x=c(1,2,3,4), y=c(0.92,0.86,0.62,0.66),label=c("a", "a", "b", "b")) + 
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Adult mass (g)") + ylim(0,1) + 
  scale_x_discrete(labels= my_labels)
adult_mass_plot

#Size graph
ggarrange(melanin_plot, food_plot, mass_plot, adult_mass_plot, 
          font.label = list(size=14, family="Times New Roman"),labels=c("A", "B", "C", "D"),
          nrow=2,ncol=2,hjust=-8, align="hv")

devo_plot<-ggplot(tradeoff_data_adult, aes(x=factor(Full_Treatment, level = treatment_order), y=Devo_Total)) + 
  geom_boxplot(aes(color=Full_Treatment), size=0.7, outlier.shape=NA) +
  geom_jitter(aes(color=Full_Treatment), width=0.25, alpha=0.5) + 
  stat_summary(fun.y="mean", color=treatment_colors, shape=18, size=1) + 
  scale_color_manual(values=treatment_colors2) + 
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Development Time") + 
  scale_x_discrete(labels= my_labels)
devo_plot


#### Part 1: Immunity #### 
#Visualize data
ggplot(aes(x=Concentration_log), data=tradeoff_data_immune) + geom_histogram()
ggplot(aes(x=Total_log), data=tradeoff_data_immune) + geom_histogram()
ggplot(aes(x=Food_Change), data=tradeoff_data_immune) + geom_histogram()

#Statistics
tradeoff_data_immune %>%
  group_by(Full_Treatment) %>%
  summarize(conc_mean = mean(Concentration_log, na.rm = TRUE),
            conc_sd = sd(Concentration_log, na.rm = TRUE),
            conc_median = median(Concentration_log, na.rm = TRUE),
            total_mean = mean(Total_log, na.rm = TRUE),
            total_sd = sd(Total_log, na.rm = TRUE),
            total_median = median(Total_log, na.rm = TRUE)) %>%
  as.data.frame
            
Concentration_mod<-lmer(Concentration_log ~  Percent_G + Diet_Treatment + Mass + (1|Family), data=tradeoff_data_immune)
summary(Concentration_mod)
emmeans(Concentration_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Concentration_mod)) 
plot(Concentration_mod)

Total_mod<-lm(Total_log ~ Photo_Treatment + Diet_Treatment + Photo_Treatment*Diet_Treatment, data=tradeoff_data_immune)
summary(Total_mod)
emmeans(Total_mod, specs="Full_Treatment") %>% pairs(adjust="tukey")
qqnorm(residuals(Concentration_mod)) 

Concentration_low <- lmer(Concentration_log ~ Percent_G + (1|Family), data=lowdiet_immune)
summary(Concentration_low)
MuMIn::r.squaredGLMM(Concentration_low)

Concentration_high<-lmer(Concentration_log ~  Percent_G + (1|Family), data=highdiet_immune)
summary(Concentration_high)
MuMIn::r.squaredGLMM(Concentration_high)


#Figures
conc_plot<-ggplot(tradeoff_data_immune, aes(x=factor(Full_Treatment, level = treatment_order2), y=Concentration_log)) + 
  geom_boxplot(size=0.7, outlier.shape=NA) +
  geom_jitter(width=0.25, alpha=0.5) + 
  stat_summary(fun.y="mean", shape=18, size=1) + 
  annotate(geom="text", x=c(1,2,3,4), y=c(2.5,2.3,2.2,1.8),label=c("a", "b", "bc", "c")) + 
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("log Melanin concentration\n(ug melanin/g wire)") + 
  scale_x_discrete(labels= my_labels2)
conc_plot

conc_plot2<-ggplot(tradeoff_data_immune, aes(x=Percent_G, y=Concentration_log)) + 
  geom_point(aes(color=Diet_Treatment)) +
  geom_smooth(method="lm", aes(linetype=Diet_Treatment, color=Diet_Treatment)) + 
  scale_color_manual(values=diet_colors) + 
  theme_classic(base_size = 15)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Percent melanic area") + ylab("log Melanin concentration\n(ug melanin/g wire)")  + xlim(0,85)
conc_plot2


#### Part 2: Muscle Mass ####
#Visualize  Data
ggplot(aes(x=muscle_percent), data=tradeoff_data_adult) + geom_histogram()
ggplot(aes(x=muscle_resids), data=tradeoff_data_adult) + geom_histogram()

#Statistics
Muscle_mod <- lm(muscle_percent ~ Full_Treatment + full_mass, data=tradeoff_data_adult)
summary(Muscle_mod)
emmeans(Muscle_mod, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(Muscle_mod)) #normality of residuals

#Residuals analysis
residual_mod <- lm(log(muscle_mass)~log(full_mass), data=tradeoff_data_adult) 
resid <- residual_mod$residuals

tradeoff_data_adult <- tradeoff_data_adult |>
  mutate(muscle_resids=residual_mod$residuals)

residual_mod_2 <- lmer(muscle_resids ~ Percent_G + Diet_Treatment + Percent_G*Diet_Treatment + (1|Family), data=tradeoff_data_adult)
summary(residual_mod_2)
emmeans(residual_mod_2, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(residual_mod_2)) 
plot(residual_mod_2)

#Figures
#residuals
resid_plot<-ggplot(tradeoff_data_adult, aes(x=factor(Full_Treatment, level = treatment_order2), y=muscle_resids)) + 
  geom_boxplot(size=0.7, outlier.shape=NA) +
  geom_jitter(width=0.25, alpha=0.5) + 
  stat_summary(fun.y="mean", shape=18, size=1) + 
  annotate(geom="text", x=c(1,2,3,4), y=c(0.23,0.23,0.40,0.19),label=c("a", "a", "ab", "b")) + 
  theme_classic(base_size = 16)+ theme(legend.position="none",text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Muscle mass residuals") +
  scale_x_discrete(labels= my_labels2) 
resid_plot

resid_plot2<-ggplot(tradeoff_data_adult, aes(x=Percent_G, y=muscle_resids)) + 
  geom_point(aes(color=Diet_Treatment)) +
  geom_smooth(method="lm", aes(color=Diet_Treatment), linetype="dashed", show.legend = F) + 
  scale_color_manual(values=diet_colors) + 
  theme_classic(base_size = 15)+ theme(text=element_text(family="Times New Roman")) + 
  xlab("Percent melanic area") + ylab("Muscle mass residuals") +  xlim(0,85) + labs(color="Diet Treatment")
resid_plot2

#raw muscle mass
muscle_plot<-ggplot(tradeoff_data_adult, aes(x=factor(Full_Treatment, level = treatment_order), y=muscle_percent)) + 
  geom_boxplot(aes(color=Full_Treatment), size=0.7, outlier.shape=NA) +
  geom_jitter(aes(color=Full_Treatment), width=0.25, alpha=0.5) + 
  stat_summary(fun.y="mean", color=treatment_colors, shape=18, size=1) + 
  scale_color_manual(values=treatment_colors2) + 
  #annotate(geom="text", x=c(1,2,3,4), y=c(17.4,19,18,16.3),label=c("a", "b", "a", "a")) + #
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Muscle (g)") + 
  scale_x_discrete(labels= my_labels)
muscle_plot


#### Part 3: Wings ####
#Visualize data
ggplot(aes(x=fore_darkness), data=tradeoff_data_adult) + geom_histogram()
ggplot(aes(x=hind_darkness), data=tradeoff_data_adult) + geom_histogram()
ggplot(aes(x=log(hind_percent)), data=tradeoff_data_adult) + geom_histogram() #right skewed

#Models
forewing_mod <- lmer(fore_darkness~Percent_G + Diet_Treatment + Percent_G*Diet_Treatment +(1|Family), data=tradeoff_data_adult)
summary(forewing_mod)
emmeans(forewing_mod, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(forewing_mod)) 
plot(forewing_mod)

hindwing_dark_mod <- lmer(hind_darkness~Melanin_Metric+ Diet_Treatment + Melanin_Metric*Diet_Treatment +(1|Family), data=tradeoff_data_adult)
summary(hindwing_dark_mod)
emmeans(hindwing_dark_mod, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(hindwing_dark_mod))
plot(hindwing_dark_mod)

hindwing_dark_low <- lmer(hind_darkness~Percent_G + (1|Family), data=lowdiet_adult)
summary(hindwing_dark_low)
MuMIn::r.squaredGLMM(hindwing_dark_low)

hindwing_dark_high <- lmer(hind_darkness~Percent_G + (1|Family), data=highdiet_adult)
summary(hindwing_dark_high)
MuMIn::r.squaredGLMM(hindwing_dark_high)


hindwing_area_mod <- lmer(hind_percent~ Percent_G + Diet_Treatment + Percent_G*Diet_Treatment +(1|Family), data=tradeoff_data_adult)
summary(hindwing_area_mod)
emmeans(hindwing_area_mod, specs="Full_Treatment") |> pairs(adjust="tukey")
qqnorm(residuals(hindwing_area_mod)) 
plot(hindwing_area_mod)

wing_tradeoffs1 <- lm(hind_darkness~hind_percent, data=lowhdiet_adult)
summary(wing_tradeoffs1)

wing_tradeoffs2 <- lm(hind_darkness~fore_darkness, data=lowdiet_adult)
summary(wing_tradeoffs2)

wing_tradeoffs3<- lm(hind_percent~muscle_resids, data=tradeoff_data_adult)
summary(wing_tradeoffs3)

ggplot(tradeoff_data_adult, aes(x=hind_darkness, y=fore_darkness)) + 
  geom_point() +
  geom_smooth(method="lm") + 
  theme_classic(base_size = 15)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Hind %") + ylab("Forewing darkness") 

#Figures
forewing_plot<-ggplot(tradeoff_data_adult, aes(x=factor(Full_Treatment, level = treatment_order2), y=fore_darkness)) + 
  geom_boxplot(size=0.7, outlier.shape=NA) +
  geom_jitter(width=0.25, alpha=0.5) +
  stat_summary(fun.y="mean", shape=18, size=1) + 
  annotate(geom="text", x=c(1,2,3,4), y=c(232,230,233,229),label=c("a", "ab", "bc", "c")) + 
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Forewing darkness") + 
  scale_x_discrete(labels= my_labels2)
forewing_plot

forewing_plot2<-ggplot(tradeoff_data_adult, aes(x=Percent_G, y=fore_darkness)) + 
  geom_point(aes(color=Diet_Treatment)) +
  geom_smooth(method="lm", linetype="dashed", aes(color=Diet_Treatment)) + 
  scale_color_manual(values=diet_colors) + 
  theme_classic(base_size = 15)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Percent melanic area") + ylab("Forewing darkness") + xlim(0,85)
forewing_plot2

hindwing_dark_plot<-ggplot(tradeoff_data_adult, aes(x=factor(Full_Treatment, level = treatment_order2), y=hind_darkness)) + 
  geom_boxplot(size=0.7, outlier.shape=NA) +
  geom_jitter(width=0.25, alpha=0.5) +
  stat_summary(fun.y="mean", shape=18, size=1) + 
  annotate(geom="text", x=c(1,2,3,4), y=c(14.9,12.8,14.8,13.9),label=c("a", "b", "ab", "b")) + 
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) +
  xlab("Treatment") + ylab("Hindwing darkness") + 
  scale_x_discrete(labels= my_labels2)
hindwing_dark_plot

hindwing_dark_plot2<-ggplot(tradeoff_data_adult, aes(x=Percent_G, y=hind_darkness)) + 
  geom_point(aes(color=Diet_Treatment)) +
  geom_smooth(method="lm", aes(linetype=Diet_Treatment, color=Diet_Treatment)) + 
  scale_color_manual(values=diet_colors) + 
  theme_classic(base_size = 15)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Percent melanic area") + ylab("Hindwing darkness")  + xlim(0,85)
hindwing_dark_plot2

hindwing_area_plot<-ggplot(tradeoff_data_adult, aes(x=factor(Full_Treatment, level = treatment_order2), y=hind_percent)) + 
  geom_boxplot(size=0.7, outlier.shape=NA) +
  geom_jitter(width=0.25, alpha=0.5) +
  stat_summary(fun.y="mean", shape=18, size=1) + 
  annotate(geom="text", x=c(1,2,3,4), y=c(49,45,48,41),label=c("a", "a", "b", "b")) + 
  theme_classic(base_size = 16)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Treatment") + ylab("Hindwing proportion\nmelanic area (%)") + 
  scale_x_discrete(labels= my_labels2)
hindwing_area_plot

hindwing_area_plot2<-ggplot(tradeoff_data_adult, aes(x=Percent_G, y=hind_percent)) + 
  geom_point(aes(color=Diet_Treatment)) +
  geom_smooth(method="lm", linetype="dashed", aes(color=Diet_Treatment)) + 
  scale_color_manual(values=diet_colors) + 
  theme_classic(base_size = 15)+ theme(legend.position="none", text=element_text(family="Times New Roman")) + 
  xlab("Percent melanic area") + ylab("Hindwing percent\nmelanic area")  + xlim(0,85)
hindwing_area_plot2



figure_5<-ggdraw() +
  draw_plot(conc_plot2, x = 0.08, y = .5, width=0.35, height=0.46) +
  draw_plot(resid_plot2, x = 0.45, y = .5, width=0.50, height=0.46) +
  draw_plot(forewing_plot2, x = 0, y = 0, width=0.33, height=0.46) +
  draw_plot(hindwing_area_plot2, x = 0.33, y = 0,width=0.33, height=0.46) +
  draw_plot(hindwing_dark_plot2,  x = 0.66, y = 0, width=0.33, height=0.46) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15, family="Times New Roman",
                  x = c(0.06, 0.45, 0, 0.33, 0.66), y = c(1, 1, 0.49, 0.49, 0.49))
figure_5
             

figure_s2<-ggdraw() +
  draw_plot(conc_plot, x = 0.13, y = .5, width=0.34, height=0.46) +
  draw_plot(resid_plot, x = 0.55, y = .5, width=0.34, height=0.46) +
  draw_plot(forewing_plot, x = 0, y = 0, width=0.33, height=0.46) +
  draw_plot(hindwing_area_plot, x = 0.33, y = 0,width=0.33, height=0.46) +
  draw_plot(hindwing_dark_plot, x = 0.66, y = 0, width=0.33, height=0.46) +
  draw_plot_label(label = c("A", "B", "C", "D", "E"), size = 15, family="Times New Roman",
                  x = c(0.1, 0.55, 0, 0.31, 0.66), y = c(1, 1, 0.49, 0.49, 0.49))
figure_s2
          

ggplot(tradeoff_data_immune, aes(x=Avg_Percent_Melanin, y=Concentration_log)) + 
  geom_point(aes(color=Full_Treatment)) +
  geom_smooth(aes(color=Full_Treatment), method=lm, se=TRUE, fill="gray80") + 
  scale_color_manual(values=treatment_colors2) +
  theme_classic(base_size = 16)+ theme(text=element_text(family="Times New Roman")) + #text font and size, remove legend
  xlab("Melanin") + ylab("Forewing darkness") 



#### Family differences ####

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

#by diet
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

#by photo
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




ggplot(tradeoff_data, aes(x=Photo_Treatment, y=Percent_G)) + 
  stat_summary(aes(group=Diet_Treatment,color=Diet_Treatment), fun=mean, shape="square", size=0.7) +
  stat_summary(aes(group=Diet_Treatment,color=Diet_Treatment), geom="errorbar",  fun.data= "mean_se",  width = 0.1, size = 0.5) + 
  stat_summary(aes(group=Diet_Treatment,color=Diet_Treatment), geom = "line", fun = mean, size=1) +
  theme_classic(base_size = 20)+ theme(text=element_text(family="Times New Roman")) + 
  xlab("Diet") + ylab("Percent melanin") 


