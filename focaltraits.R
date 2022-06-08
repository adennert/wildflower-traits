####ALLISON DENNERT, Wildflower Traits Analyses
####"Experimental addition of marine-derived nutrients affects wildflower traits 
####in a coastal meta-ecosystem" by A. Dennert, E. Elle, J. Reynolds

#### Last updated; June 8th, 2022

#This R code analyzes/visualizes 1) isotopes and % nitrogen 2) total leaf area,
#and 3) whole plant floral display size, 4) and seed set  in Douglas' Aster (ASSU), 
#Yarrow (ACMI), silverweed (POAN), and common red paintbrush (CAMI) across four 
#treatments: salmon, algae, both, control. 25 randomized blocks were used for 
#all species, with all 4 treatments within a block. Up to 3 individuals of each 
#species were measured per treatment plot. Individuals are nested within 
#species, nested within plot, nested within block, except for the isotope data
#where one individual of each species was sampled

####0. PACKAGE LOADING####
library(dplyr)
library(lme4)
library(visreg)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(plotrix)
library(sjPlot)
library(MuMIn)
library(bootpredictlme4)
library(tidyr)
library(car)
library(MASS)
library(glmmTMB)
library(DHARMa)
library(emmeans)
library(cowplot)

####00. READ/CLEAN THE DATA, SET UP CODE####
f <- read.csv("Raw Data/focaltraits.csv", 
              strip.white = TRUE)
#contains N = 2195

# ensure that each variable type is correct; change variable type if need be,
# rename variables for clarity, and set control group first (as intercept)
str(f)
f$year <- as.factor(f$year)
f$block <- as.factor(f$block) 
f$indiv <- as.factor(f$indiv)
f <- f %>% rename(treatment = plot)
levels(f[,'treatment']) 
levels(f[,'plot.id'])
f$treatment = factor(f$treatment,levels(f$treatment)[c(4,1,2,3)])
levels(f[,'treatment']) 
f$display.size.4 <- as.integer(f$display.size.4) 
f$seed.set.2 <- as.integer(f$seed.set.2) 
f$seed.set.3 <- as.integer(f$seed.set.3) 
str(f)

# calculate leaf area, total leaf area, total display area, investment,
# and total seed set
f <- f %>% dplyr::mutate(leaf.area = leaf.length*leaf.width) %>% 
  dplyr::mutate(total.leaf.area = leaf.area*total.leaves)

#replace NAs with zeros for ease of display area calculation only
f <- f %>% mutate(replace.display.1 = ifelse(is.na(display.size.1), 0, display.size.1))
f <- f %>% mutate(replace.display.2 = ifelse(is.na(display.size.2), 0, display.size.2))
f <- f %>% mutate(replace.display.3 = ifelse(is.na(display.size.3), 0, display.size.3))
f <- f %>% mutate(replace.display.4 = ifelse(is.na(display.size.4), 0, display.size.4))
f <- f %>% mutate(replace.display.5 = ifelse(is.na(display.size.5), 0, display.size.5))
f <- f %>% mutate(replace.display.6 = ifelse(is.na(display.size.6), 0, display.size.6))
f <- f %>% mutate(replace.display.7 = ifelse(is.na(display.size.7), 0, display.size.7))
f <- f %>% mutate(replace.display.8 = ifelse(is.na(display.size.8), 0, display.size.8))

#calculate area of each flower and add them together to calculate total display
#and remove zero values as those were formerly NAs with no data
f <- f %>% 
  dplyr::mutate(total.display.area = (pi*(replace.display.1/2)^2) + 
                                            (pi*(replace.display.2/2)^2) +
                                            (pi*(replace.display.3/2)^2) +
                                            (pi*(replace.display.4/2)^2) +
                                            (pi*(replace.display.5/2)^2) +
                                            (pi*(replace.display.6/2)^2) +
                                            (pi*(replace.display.7/2)^2) +
                                            (pi*(replace.display.8/2)^2))

# replace NA seed values with zeros for ease of total seed calculation only
f <- f %>% mutate(replace.seed.1 = ifelse(is.na(seed.set.1), 0, seed.set.1))
f <- f %>% mutate(replace.seed.2 = ifelse(is.na(seed.set.2), 0, seed.set.2))
f <- f %>% mutate(replace.seed.3 = ifelse(is.na(seed.set.3), 0, seed.set.3))
f <- f %>% mutate(replace.seed.4 = ifelse(is.na(seed.set.4), 0, seed.set.4))
f <- f %>% mutate(replace.seed.5 = ifelse(is.na(seed.set.5), 0, seed.set.5))
f <- f %>% mutate(replace.seed.6 = ifelse(is.na(seed.set.6), 0, seed.set.6))
f <- f %>% mutate(replace.seed.7 = ifelse(is.na(seed.set.7), 0, seed.set.7))
f <- f %>% mutate(replace.seed.8 = ifelse(is.na(seed.set.8), 0, seed.set.8))
f <- f %>% mutate(replace.seed.9 = ifelse(is.na(seed.set.9), 0, seed.set.9))
f <- f %>% mutate(replace.seed.10 = ifelse(is.na(seed.set.10), 0, seed.set.10))
f <- f %>% mutate(replace.seed.11 = ifelse(is.na(seed.set.11), 0, seed.set.11))
f <- f %>% mutate(replace.seed.12 = ifelse(is.na(seed.set.12), 0, seed.set.12))
f <- f %>% mutate(replace.seed.13 = ifelse(is.na(seed.set.13), 0, seed.set.13))
f <- f %>% mutate(replace.seed.14 = ifelse(is.na(seed.set.14), 0, seed.set.14))

#add seeds together to calculate total seeds
f <- f %>% dplyr::mutate(total.seeds = replace.seed.1 + replace.seed.2 + replace.seed.3 + 
                                    replace.seed.4 + replace.seed.5 + replace.seed.6 + 
                                    replace.seed.7 + replace.seed.8 + replace.seed.9 + 
                                    replace.seed.10 + replace.seed.11 + replace.seed.12 + 
                                    replace.seed.13 + replace.seed.14) 

#round total seed set to *UP* to nearest whole number
f$total.seeds <- ceiling(f$total.seeds)

# read and clean isotope data from 2017, 2018, 2019
iso <- read.csv("Raw Data/2017leafisotopes.csv", stringsAsFactors=TRUE)
#N = 629
head(iso)

#create new variables for C:N ratios and %N
iso <- iso %>% dplyr::mutate(C.N = C.ug/N.ug) %>% 
  dplyr::mutate(percent.N = ((N.ug*0.001)/sample.mg)*100)

str(iso)
iso$year <- as.factor(iso$year)
iso$block <- as.factor(iso$block)
levels(iso[,'treatment']) 
iso$treatment = factor(iso$treatment,levels(iso$treatment)[c(4,1,2,3)])
levels(iso[,'treatment']) 

#organize labels, line separation, color palletes
#dplyr::recode(treatment, 
              #"D" = "Control", "B" = "Fucus", 
              #"A" = "Salmon", "C" = "Salmon& Fucus"))
#levels(pl$treatment) <- gsub(" ", "\n", levels(pl$treatment))
#levels(pl$treatment) <- gsub("Salmon&", "Salmon &", levels(pl$treatment))
#levels(pl$treatment)

labels.sp <- c(ACMI = "Yarrow", ASSU = "Douglas' Aster", 
               CAMI = "Paintbrush", POAN = "Silverweed")
plantPalette <- c("#000000", "#9999CC", "#FF3333", "#CCCC33")

#make a summary function for means and SE in plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m - (1.96*std.error(x))
  ymax <- m + (1.96*std.error(x))
  return(c(y=m,ymin=ymin,ymax=ymax))
}


####1. ISOTOPES ####
####PLOTTING ####
ggplot(iso, aes(x = treatment, y = d15N, colour =  species)) +
  geom_jitter(position=position_jitter(0.2), alpha = 0.25, na.rm = TRUE, 
              size = 2) +
  facet_wrap(.~species, ncol = 4, labeller = 
               labeller(species = labels.sp)) +
  stat_summary(fun.data=data_summary, color="black", size = 0.7) +
  scale_colour_manual(values = plantPalette) +
  #expand_limits(y = 0) +
  theme_classic(24) +
  labs(x = "Treatment", y = expression(paste(delta^{15} ~ N))) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(hjust = 1), 
        legend.position="none") 

ggplot(iso, aes(x = treatment, y = percent.N, colour =  species)) +
  geom_jitter(position=position_jitter(0.2), alpha = 0.25, na.rm = TRUE, 
              size = 3) +
  facet_wrap(.~species, ncol = 4, labeller = 
               labeller(species = labels.sp)) +
  stat_summary(fun.data=data_summary, color="black", size = 0.75) +
  scale_colour_manual(values = plantPalette) +
  #expand_limits(y = 0) +
  theme_classic(32) +
  labs(x = "Treatment", y = "Leaf % Nitrogen") +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = -1.5),
        legend.position="none", axis.title.x = element_text(margin = margin(t = 35)))


ggsave("Figures/N_pubfig.png", height=9, width=16)


####MODELLING####
##MODEL THE DATA
hist(iso$d15N, breaks = 500) 

d15N.model <- lmer(d15N ~ treatment * species * year 
                   - treatment:species:year 
                   + (1|block), 
                   data = iso, 
                   na.action = "na.omit", REML = TRUE)

summary(d15N.model)

#Check Model: Use DHARMa package to evaluate fit
##GLMM Model Diagnostics: Simulate residuals (DHARMa Package)##
res <- simulateResiduals(fittedModel = d15N.model, n = 500)
plot(res)
#sig deviation very likely due to large sample size, line is also straight;
#see DHARMa vignette for info and justifcation

#The lines in the resid vs predict should be horizontal at 0.25/0.5/0.75.
#Warning: with DHARMa diagnosing glmmTMB, package author states the following:
#"With strong random effects, this can sometimes create diagonal patterns from 
#bottom left to top right in the res ~ pred plot. All other tests and plots 
#should work as desired"

#Check Individual Predictors
par(mfrow = c(1,1))
plotResiduals(res$scaledResiduals, iso$treatment)
plotResiduals(res$scaledResiduals, iso$species)
plotResiduals(res$scaledResiduals, iso$year)
#box plots look reasonable despite negative Levene test; red likely due to sample
#size again, see DHARMa vignette for justification

# Goodness of fit test
testUniformity(simulationOutput = res) 
# Deviation issue most likely reflective of large sample size (DHARMa vignette),
# dispersion looks good, line very straight given sample size

#Check for outliers
testOutliers(simulationOutput = res) #number of outliers due to number of res 
#simulations; when I experimentally decreased the number of simulations from 500 
#to 400 the outlier issue went away; see DHARMa vignette

# Check Dispersion
testDispersion(res) # No over/underdispersion 
#(Small p-value (<0.05) indicates dispersion problem)

#Full coeff plot for referencee
plot_model(d15N.model, type = "eff",  sort = TRUE,
           vline.color = "grey") + theme_classic()

#Plot least squares means for each level of each categorical variable
emiso <- emmeans(d15N.model, ~treatment*year*species)
emiso
plot(emiso, comparison = TRUE)
#calculate Dunnett-adjusted p-values using Dunnett comparison to control for each
#species and year
contrast(emiso, method = "dunnett", adjust = "dunnett", by = c("species", "year"))


####PUBLIATION FIGURE####

#Plot least squares mean for each level of categorical variable, but 
#average across all three years
emisoyear <- emmeans(d15N.model, ~treatment*species)
emisosummary <- as.data.frame(emisoyear)
plot(emisoyear, comparison = TRUE)


ggplot(NULL, aes(x = treatment, y = d15N, colour =  species)) +
  geom_jitter(data = iso, position=position_jitter(0.2), alpha = 0.25, na.rm = TRUE, 
              size = 3) +
  facet_wrap(.~species, ncol = 4, scales="free_y", labeller = 
               labeller(species = labels.sp)) +
  geom_pointrange(inherit.aes = FALSE, data = emisosummary, aes(x = treatment, 
                                                                y = emmean, ymin = lower.CL, ymax = upper.CL), color = "black", 
                  size = 0.75) +
  scale_colour_manual(values = plantPalette) +
  theme_classic(32) +
  labs(x = "Treatment", y = expression(paste(Leaf ~ delta^{15} ~ N))) +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = -1.5),
        legend.position="none", axis.title.x = element_text(margin = margin(t = 35)))

ggsave("Figures/iso_pubfig.png", height=9, width=16)






####2. PERCENT N ####
####PLOTTING ####

ggplot(iso, aes(x = treatment, y = percent.N, colour =  species)) +
  geom_jitter(position=position_jitter(0.2), alpha = 0.25, na.rm = TRUE, 
              size = 3) +
  facet_wrap(.~species, ncol = 4, labeller = 
               labeller(species = labels.sp)) +
  stat_summary(fun.data=data_summary, color="black", size = 0.75) +
  scale_colour_manual(values = plantPalette) +
  #expand_limits(y = 0) +
  theme_classic(32) +
  labs(x = "Treatment", y = "Leaf % Nitrogen") +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = -1.5),
        legend.position="none", axis.title.x = element_text(margin = margin(t = 35)))


ggsave("Figures/N_pubfig.png", height=9, width=16)

ggplot(subset(iso, species %in% c("ACMI")), aes(x = d15N, y = percent.N, colour =  treatment)) +
  geom_jitter(position=position_jitter(0.2), alpha = 0.25, na.rm = TRUE, 
              size = 3) +
  facet_wrap(.~species, ncol = 2, labeller = 
               labeller(species = labels.sp),  scales = "free_x") +
  #stat_summary(fun.data=data_summary, color="black", size = 0.75) +
  scale_colour_manual(values = plantPalette) +
  #expand_limits(y = 0) +
  geom_smooth(method = "lm") +
  theme_classic(32) +
  labs(x = "Nitrogen-15", y = "Leaf % Nitrogen") +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = -1.5),
        axis.title.x = element_text(margin = margin(t = 35)))



####MODELLING####
##MODEL THE DATA
hist(iso$percent.N, breaks = 500) 

percent.N.model <- lmer(percent.N ~ treatment * species * year 
                        - treatment:species:year 
                        + (1|block), 
                        data = iso, 
                        na.action = "na.omit", REML = TRUE)

summary(percent.N.model)

#Check Model: Use DHARMa package to evaluate fit
##GLMM Model Diagnostics: Simulate residuals (DHARMa Package)##
res <- simulateResiduals(fittedModel = percent.N.model, n = 500)
plot(res)
#sig deviation likely due to sample size, line is also straight;
#see DHARMa vignette for info

#The lines in the resid vs predict should be horizontal at 0.25/0.5/0.75.
#Warning: with DHARMa diagnosing glmmTMB, package author states the following:
#"With strong random effects, this can sometimes create diagonal patterns from 
#bottom left to top right in the res ~ pred plot. All other tests and plots 
#should work as desired"

#Check Individual Predictors
par(mfrow = c(1,1))
plotResiduals(res$scaledResiduals, iso$treatment)
plotResiduals(res$scaledResiduals, iso$species)
plotResiduals(res$scaledResiduals, iso$year)
#box plots look reasonable despite negative Levene test; red likely due to sample
#size again, see DHARMa vignette for justification

# Goodness of fit test
testUniformity(simulationOutput = res) 
# Deviation issue most likely reflective of large sample size (DHARMa vignette),
# dispersion looks good, line very straight given sample size

#Check for outliers
testOutliers(simulationOutput = res) #number of outliers due to number of res 
#simulations; when I experimentally decreased the number of simulations from 500 
#to 400 the outlier issue went away; see DHARMa vignette

# Check Dispersion
testDispersion(res) # No over/underdispersion 
#(Small p-value (<0.05) indicates dispersion problem)

#Full coeff plot for referencee
plot_model(percent.N.model, type = "eff",  sort = TRUE,
           vline.color = "grey") + theme_classic()

#Plot least squares means for each level of each categorical variable
emN <- emmeans(percent.N.model, ~treatment*year*species)
emN
plot(emN, comparison = TRUE)
#calculate Dunnett-adjusted p-values using Dunnett comparison to control for each
#species and year
contrast(emN, method = "dunnett", adjust = "dunnett", by = c("species", "year"))

emmeans(percent.N.model, ~species)

####PUBLIATION FIGURE####

#Plot least squares mean for each level of categorical variable, but 
#average across all three years
emNyear <- emmeans(percent.N.model, ~treatment*species)
emNsummary <- as.data.frame(emNyear)
plot(emNyear, comparison = TRUE)


ggplot(NULL, aes(x = treatment, y = percent.N, colour =  species)) +
  geom_jitter(data = iso, position=position_jitter(0.2), alpha = 0.25, na.rm = TRUE, 
              size = 3) +
  facet_wrap(.~species, ncol = 4, scales="free_y", labeller = 
               labeller(species = labels.sp)) +
  geom_pointrange(inherit.aes = FALSE, data = emNsummary, aes(x = treatment, 
                                                              y = emmean, ymin = lower.CL, ymax = upper.CL), color = "black", 
                  size = 0.75) +
  scale_colour_manual(values = plantPalette) +
  theme_classic(32) +
  labs(x = "Treatment", y = "Leaf % Nitrogen") +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = -1.5),
        legend.position="none", axis.title.x = element_text(margin = margin(t = 35)))

ggsave("Figures/percentN_pubfig.png", height=9, width=16)





####3. TOTAL LEAF AREA####
####PLOTTING####
##PLOT THE RAW DATA
ggplot(f, aes(x = treatment, y = log(total.leaf.area), colour =  species)) +
  geom_jitter(position=position_jitter(0.2), alpha = 0.25, na.rm = TRUE, 
              size = 2) +
  facet_wrap(.~species, ncol = 4, scales="free_y", labeller = 
               labeller(species = labels.sp)) +
  stat_summary(fun.data=data_summary, color="black", size = 0.7) +
  scale_colour_manual(values = plantPalette) +
  #expand_limits(y = 0) +
  theme_classic(28) +
  labs(x = "Treatment", y = expression(paste('log(Total leaf area (mm'^2,')'))) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, 
                                                                   hjust = 1), 
        legend.position="none") 
ggsave("Figures/totalleafarea_rawdata.png", height=9, width=16)


####MODELLING####
##MODEL THE DATA
hist(f$total.leaf.area, breaks = 500) #need to remove impossibly high outlier values
filtered.f.leaf <- f %>%dplyr::filter(total.leaf.area < 100000)
hist(filtered.f.leaf$total.leaf.area, breaks = 500)
#removes 15 of 2195 values

total.leaf.area.model <- glmmTMB(total.leaf.area ~ treatment * species * year 
                                    - treatment:species:year 
                                    + (1|block/plot.id), 
                                    data = filtered.f.leaf, 
                                    family = Gamma(link = "log"), ziformula = ~0, 
                                    dispformula = ~year*treatment*species,
                                    na.action = "na.omit",  
                                    control = glmmTMBControl(optimizer = nlminb))

summary(total.leaf.area.model)

#Check Model: Use DHARMa package to evaluate fit
##GLMM Model Diagnostics: Simulate residuals (DHARMa Package)##
res <- simulateResiduals(fittedModel = total.leaf.area.model, n = 500)
plot(res)
#sig deviation likely due to sample size along, line is also straight;
#see DHARMa vignette for info

#The red lines in the resid vs predict should be horizontal at 0.25/0.5/0.75.
#Warning: with DHARMa diagnosing glmmTMB, package author states the following:
#"With strong random effects, this can sometimes create diagonal patterns from 
#bottom left to top right in the res ~ pred plot. All other tests and plots 
#should work as desired"

#Check Individual Predictors
f2 <- tidyr::drop_na(filtered.f.leaf, total.leaf.area)
par(mfrow = c(1,1))
plotResiduals(res$scaledResiduals, f2$treatment)
plotResiduals(res$scaledResiduals, f2$species)
plotResiduals(res$scaledResiduals, f2$year)
#box plots look reasonable despite negative Levene test; red likely due to sample
#size again, see DHARMa vignette for justification

# Goodness of fit test
testUniformity(simulationOutput = res) 
# Deviation issue most likely reflective of large sample size (DHARMa vignette),
# dispersion looks good, line very straight given sample size

#Check for outliers
testOutliers(simulationOutput = res) #number of outliers likely die to sample size
#will not remove any simply to improve fit, other than prior removal of 
#impossibly high values

# Check Dispersion
testDispersion(res) # No over/underdispersion 
#(Small p-value (<0.05) indicates dispersion problem)

#Full coeff plot for referencee
plot_model(total.leaf.area.model, type = "eff",  sort = TRUE,
           vline.color = "grey") + theme_classic()
           
#Plot least squares means for each level of each categorical variable
emleaf <- emmeans(total.leaf.area.model, ~treatment*year*species)
emleaf
plot(emleaf, comparison = TRUE)
#calculate Dunnett-adjusted p-values using Dunnett comparison to control for each
#species and year
contrast(emleaf, method = "dunnett", adjust = "dunnett", by = c("species", "year"))



##coeff with NS interactions removed
#plot_model(total.leaf.area.model, type = "est", sort.est = FALSE, vline.color = "grey",
#           terms = c("treatmentA", "treatmentB", "treatmentC", "speciesASSU", 
#                     "speciesCAMI","speciesPOAN", "year2018", "year2019",
#                     "treatmentC:speciesASSU", "treatmentA:speciesASSU",
#                     "speciesPOAN:year2018", "treatmentA:speciesPOAN",
#                     "speciesPOAN:year2019"),
#           group.terms = c(1,2,1,2,2,1,2,1,1,1,1,1,1),
#           order.terms = c(1,2,3,4,5,6,7,8,9,10,11,12,13),
#           colors = c("black", "grey"), title = 
#             "Coeff Estimates: Total Leaf Area Model", dot.size = 5,
#           line.size = 2) +
#  theme_classic(30) 
#ggsave("Figures/modeloutputs_totalleafarea.png",  height=9, width=16)

####PUBLICATION FIGURE####

#Plot least squares mean for each level of categorical variable, but 
#average across all three years
emleafyear <- emmeans(total.leaf.area.model, ~treatment*species)
emleafsummary <- as.data.frame(emleafyear)
plot(emleafyear)

ggplot(NULL, aes(x = treatment, y = log(total.leaf.area), colour =  species)) +
  geom_jitter(data = filtered.f.leaf, position=position_jitter(0.2), alpha = 0.2, na.rm = TRUE, 
              size = 3) +
  facet_wrap(.~species, ncol = 4, scales="free_y", labeller = 
               labeller(species = labels.sp)) +
  geom_pointrange(inherit.aes = FALSE, data = emleafsummary, aes(x = treatment, 
              y = emmean, ymin = lower.CL, ymax = upper.CL), color = "black", 
              size = 0.75) +
  scale_colour_manual(values = plantPalette) +
  theme_classic(32) +
  labs(x = "Treatment", y = expression(paste('log(Leaf area estimate (mm'^2,'))'))) +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = -1.5),
        legend.position="none", axis.title.x = element_text(margin = margin(t = 35)))
  
ggsave("Figures/totalleafarea_pubfig.png", height=9, width=16)



####4.TOTAL DISPLAY AREA####
####PLOTTING####
##PLOT THE RAW DATA
ggplot(filtered.display.f, aes(x = treatment, y = (total.display.area), colour =  species)) +
  geom_jitter(position=position_jitter(0.2), alpha = 0.25, na.rm = TRUE, 
              size = 2) +
  facet_wrap(.~species, ncol = 4, scales="free_y", labeller = 
               labeller(species = labels.sp)) +
  stat_summary(fun.data=data_summary, color="black", size = 0.7) +
  scale_colour_manual(values = plantPalette) +
  #expand_limits(y = 0) +
  theme_classic(28) +
  labs(x = "Treatment", y = expression(paste('log(Total Display area (mm'^2,'))'))) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, 
                                                                   hjust = 1), 
        legend.position="none") 
ggsave("Figures/totaldisplayarea_rawdata.png", height=9, width=16)


####MODELLING####
##MODEL THE DATA
hist(f$total.display.area, breaks = 500)

filtered.display.f <- f %>% 
  dplyr::filter(total.display.area>0) #excludes NA and zero values as they are 
#meaningless due to prior transformation, N = 1346 display measurements
hist(filtered.display.f$total.display.area, breaks = 500)

#remove outlier values above 10000, and impossible value of 3
filtered.display.f <- filtered.display.f %>% 
  dplyr::filter(total.display.area<8000) %>% 
  dplyr::filter(total.display.area > 8) #removes erroneous impossible values
  
#removes 8 of 1346 impossible values, now N = 1338

total.display.area.model <- glmmTMB(total.display.area ~ treatment * species * year 
                                    - treatment:species:year 
                                    + (1|block/plot.id), 
                            data = filtered.display.f, 
                            family = Gamma(link = "log"), ziformula = ~0, 
                            dispformula = ~treatment*species*year,
                            na.action = "na.omit",  
                            control = glmmTMBControl(optimizer = nlminb))

summary(total.display.area.model)
r.squaredGLMM(total.display.area.model)


#Check Model: Use DHARMa package to evaluate fit
##GLMM Model Diagnostics: Simulate residuals (DHARMa Package)##
res <- simulateResiduals(fittedModel = total.display.area.model, n = 500)
plot(res)
#deviation obs vs. expected is significant, but likely only de to sample size
#according to the dharma package vignette

#The red lines in the resid vs predict should be horizontal at 0.25/0.5/0.75.
#Warning: with DHARMa diagnosing glmmTMB, package author states the following:
#"With strong random effects, this can sometimes create diagonal patterns from 
#bottom left to top right in the res ~ pred plot. All other tests and plots 
#should work as desired"

#Check Individual Predictors
f2 <- tidyr::drop_na(filtered.display.f, total.display.area)
par(mfrow = c(1,1))
par(mfrow = c(1,1))
plotResiduals(res$scaledResiduals, f2$treatment)
plotResiduals(res$scaledResiduals, f2$species)
plotResiduals(res$scaledResiduals, f2$year)

# Goodness of fit test
testUniformity(simulationOutput = res) 
# Likely only an artifact large sample size, see DHARMa vignette

#Check for outliers
testOutliers(simulationOutput = res) #Outliers present--need to deal with this

# Check Dispersion
testDispersion(res) # No over/underdispersion 
#(Small p-value (<0.05) indicates dispersion problem)

#Full coeff plot
plot_model(total.display.area.model, type = "est", sort.est = TRUE, vline.color = "grey") +
  theme_classic()

#Plot least squares means for each level of each categorical variable
emdisplay <- emmeans(total.display.area.model, ~treatment*year*species)
emdisplay
plot(emdisplay, comparison = TRUE)
#calculate Dunnett-adjusted p-values using Dunnett comparison to control for each
#species and year
contrast(emdisplay, method = "dunnett", adjust = "dunnett", by = c("species", "year"))


##coeff with NS interactions removed
#plot_model(total.display.area.model, type = "est", sort.est = FALSE, vline.color = "grey",
#           terms = c("treatmentA", "treatmentB", "treatmentC", "speciesASSU", 
#                     "speciesCAMI","speciesPOAN", "year2018", "year2019",
#                     "treatmentB:year2018", "treatmentA:speciesASSU",
#                     "speciesCAMI:year2019", "speciesPOAN:year2018",
#                     "speciesASSU:year2018", "speciesPOAN:year2019",
#                     "speciesASSU:year2019"),
#           group.terms = c(2,1,2,1,1,1,1,1,1,1,1,1,1,1,1),
#           order.terms = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),
#           colors = c("black", "grey"), title = 
#             "Coeff Estimates: Display Area", dot.size = 5,
#           line.size = 2) +
#  theme_classic(30) 
#ggsave("Figures/modeloutputs_totaldisplayarea.png",  height=9, width=16)

####PUBLICATION FIGURE####

#Plot least squares mean for each level of categorical variable, but 
#average across all three years
emdisplayyear <- emmeans(total.display.area.model, ~treatment*species)
emdisplaysummary <- as.data.frame(emdisplayyear)
plot(emdisplayyear, comparison = TRUE)

ggplot(NULL, aes(x = treatment, y = log(total.display.area), colour =  species)) +
  geom_jitter(data = filtered.display.f, position=position_jitter(0.2), alpha = 0.2, na.rm = TRUE, 
              size = 3) +
  facet_wrap(.~species, ncol = 4, scales="free_y", labeller = 
               labeller(species = labels.sp)) +
  geom_pointrange(inherit.aes = FALSE, data = emdisplaysummary, aes(x = treatment, 
                                                                 y = emmean, ymin = lower.CL, ymax = upper.CL), color = "black", 
                  size = 0.75) +
  scale_colour_manual(values = plantPalette) +
  theme_classic(32) +
  labs(x = "Treatment", y = expression(paste('log(Display area estimate (mm'^2,'))'))) +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = -1.5),
        legend.position="none", axis.title.x = element_text(margin = margin(t = 35)))
ggsave("Figures/totaldisplayarea_pubfig.png", height=9, width=16)




####5. SEED SET####
####PLOTTING####
##PLOT THE RAW DATA
filtered.seeds.f <- f %>% dplyr::filter(total.seeds < 6000) %>% 
  dplyr::filter(total.seeds > 0)

ggplot(filtered.seeds.f, aes(x = treatment, y = log(total.seeds), colour =  species)) +
  geom_jitter(position=position_jitter(0.2), alpha = 0.25, na.rm = TRUE, 
              size = 2) +
  facet_wrap(.~species, ncol = 4, scales="free_y", labeller = 
               labeller(species = labels.sp)) +
  stat_summary(fun.data=data_summary, color="black", size = 0.7) +
  scale_colour_manual(values = plantPalette) +
  #expand_limits(y = 0) +
  theme_classic(28) +
  labs(x = "Treatment", y = expression(paste('Total seeds'))) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, 
                                                                   hjust = 1), 
        legend.position="none") 
ggsave("Figures/totalseeds_rawdata.png", height=9, width=16)

####MODELLING####
##MODEL THE DATA
hist(f$total.seeds, breaks = 500)

filtered.seeds.f <- f %>% dplyr::filter(total.seeds < 6000) %>% 
  dplyr::filter(total.seeds > 0)
hist(filtered.seeds.f$total.seeds, breaks = 500)

total.seed.model <- glmmTMB(total.seeds ~ treatment * species * year 
                                 - treatment:species:year 
                                 + (1|block/plot.id), 
                                 data = filtered.seeds.f, 
                                 family = nbinom2(link = "log"), 
                                 ziformula = ~year,
                                 dispformula = ~species*year,
                                 na.action = "na.omit", 
                                 control = glmmTMBControl(optimizer = nlminb))

summary(total.seed.model)


#Check Model: Use DHARMa package to evaluate fit
##GLMM Model Diagnostics: Simulate residuals (DHARMa Package)##
res <- simulateResiduals(fittedModel = total.seed.model, n = 500)
plot(res)

#The red lines in the resid vs predict should be horizontal at 0.25/0.5/0.75.
#Warning: with DHARMa diagnosing glmmTMB, package author states the following:
#"With strong random effects, this can sometimes create diagonal patterns from 
#bottom left to top right in the res ~ pred plot. All other tests and plots 
#should work as desired"

#Check Individual Predictors
par(mfrow = c(1,1))
plotResiduals(res$scaledResiduals, filtered.seeds.f$treatment)
plotResiduals(res$scaledResiduals, filtered.seeds.f$species)
plotResiduals(res$scaledResiduals, filtered.seeds.f$year)

# Goodness of fit test
testUniformity(simulationOutput = res) 
# No Uniformity Issue, deviation due to large sample size (see DHARMa vignette)

#Check for outliers
testOutliers(simulationOutput = res) #Outliers present, likely an artifact of high
#sample size, see DHARMa vingette for info

# Check Dispersion
testDispersion(res) # No over/underdispersion 
#(Small p-value (<0.05) indicates dispersion problem)

#Full coeff plot
plot_model(total.seed.model, type = "est", sort.est = TRUE, vline.color = "grey") +
  theme_classic()

#Plot least squares means for each level of each categorical variable
emseed <- emmeans(total.seed.model, ~treatment*year*species)
emseed
plot(emseed, comparison = TRUE)
#calculate Dunnett-adjusted p-values using Dunnett comparison to control for each
#species and year
contrast(emseed, method = "dunnett", adjust = "dunnett", by = c("species", "year"))


####PUBLICATION FIGURE####

#Plot least squares mean for each level of categorical variable, but 
#average across all three years
emseedyear <- emmeans(total.seed.model, ~treatment*species)
emseedsummary <- as.data.frame(emseedyear)
plot(emseedyear, comparison = TRUE)

ggplot(NULL, aes(x = treatment, y = log(total.seeds), colour =  species)) +
  geom_jitter(data = filtered.seeds.f, position=position_jitter(0.2), alpha = 0.2, na.rm = TRUE, 
              size = 3) +
  facet_wrap(.~species, ncol = 4, scales="free_y", labeller = 
               labeller(species = labels.sp)) +
  geom_pointrange(inherit.aes = FALSE, data = emseedsummary, aes(x = treatment, 
                                                                    y = emmean, ymin = lower.CL, ymax = upper.CL), color = "black", 
                  size = 0.75) +
  scale_colour_manual(values = plantPalette) +
  theme_classic(32) +
  labs(x = "Treatment", y = "Estimated seed set") +
  theme(axis.text.x = element_text(hjust = 0,
                                   vjust = -1.5),
        legend.position="none", axis.title.x = element_text(margin = margin(t = 35)))


ggsave("Figures/seedset_pubfig.png", height=9, width=16)

