# NSF FX data analysis for the Waiakea biomass data

# There are a number of ways we could analyze biomass growth rates
# for these trees. Each of these will be separate for each species.
# - Compute absolute and relative growth rates for each time period
# - Compute absolute and relative growth rates for each year
# - Compute absolute and relative growth rates from t0 to each tf
# For each of these, we want to analyze growth as a function of 
# treatment, which answers our initial scientific question:
# Is the species limited by N and/or P at each treatment level?
# We also want to account for the fact that growth rate varies with
# tree size, and possibly with year.

# How to handle odd datapoints? By odd datapoints, I mean ones that 
# result in negative growth or are otherwise anomalous compared to the
# previous/expected growth rate of the tree.
# There are a reasons why datapoints could be odd:
# - Measurement or data input error
# - Point loss of biomass (branches etc.) due to physical damage
# - Persistent loss of biomass due to pest damage or other illness
# The way we deal with odd datapoints would be different for each.
# - Measurement or data input error should correct themselves over time, 
# although negative relative growth rates will be tricky.
# - Point loss of biomass should be a "restart" of growth measurements.
# - Persistent loss of biomass due to pest damage or other illness should
# be removed from the analysis, as it does not address the question.

# I'm editing the "Use_growth_0X" columns in
# "HI_W_FX_Size_FoliarCNIsotope_Data.csv" to reflect these odd datapoints
# based on biomass trends.
# 1 means not odd in any way.
# 0.5 means use July 2019 if not using January 2019, but not if using January 2019.
# 0 means sufficiently ill or damaged that it shouldn't be used.
# The "Use_growth_0X" columns in 
# "HI_W_FX_Height_Diam_FoliarCNIsotope_Data.csv" were based on diameter,
# height, and canopy dimension data independently.

rm(list=ls())
library(nlme)
library(emmeans)

setwd("/Users/duncanmenge/Documents/Academia/Grants/NSF/2014_Strategies/Data/Master datafiles/")

dat <- read.csv("HI_W_FX_Size_FoliarCNIsotope_Data.csv")[1:108,1:111]

# Dates
# _01: January 2016
# _02: July 2016
# _03: June 2017
# _04: June 2018
# _05: January 2019
# _06: July 2019

# Option to print out how many are in each category (uncomment)

#print(dat[dat$Species=="GLSE" & dat$Treatment=="LN",c(33,44,59,76,87)])
#print(dat[dat$Species=="GLSE" & dat$Treatment=="MN",c(33,44,59,76,87)])
#print(dat[dat$Species=="GLSE" & dat$Treatment=="HN",c(33,44,59,76,87)])
#print(dat[dat$Species=="GLSE" & dat$Treatment=="PHN",c(33,44,59,76,87)])

#print(dat[dat$Species=="CAEQ" & dat$Treatment=="LN",c(33,44,59,76,87)])
#print(dat[dat$Species=="CAEQ" & dat$Treatment=="MN",c(33,44,59,76,87)])
#print(dat[dat$Species=="CAEQ" & dat$Treatment=="HN",c(33,44,59,76,87)])
#print(dat[dat$Species=="CAEQ" & dat$Treatment=="PHN",c(33,44,59,76,87)])

#print(dat[dat$Species=="PSCA" & dat$Treatment=="LN",c(33,44,59,76,87)])
#print(dat[dat$Species=="PSCA" & dat$Treatment=="MN",c(33,44,59,76,87)])
#print(dat[dat$Species=="PSCA" & dat$Treatment=="HN",c(33,44,59,76,87)])
#print(dat[dat$Species=="PSCA" & dat$Treatment=="PHN",c(33,44,59,76,87)])

#############################################################
################ Calculate growth rates #####################
#############################################################

# "AGR" is absolute growth rate and "RGR" is relative growth rate. 
# Each is calculated for aboveground biomass, 
# for the following time periods: 

# Annual/semiannual increments:

# AGB

# January 2016 to July 2016
dat$AGR_AGB_2016_01_2016 <- (dat$AGB_est_kg_02 - dat$AGB_est_kg_01)/((dat$Days_02 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_01_2016 <- (log(dat$AGB_est_kg_02) - log(dat$AGB_est_kg_01))/((dat$Days_02 - dat$Days_01)/365.25)

# July 2016 to June 2017
dat$AGR_AGB_2016_2017 <- (dat$AGB_est_kg_03 - dat$AGB_est_kg_02)/((dat$Days_03 - dat$Days_02)/365.25)
dat$RGR_AGB_2016_2017 <- (log(dat$AGB_est_kg_03) - log(dat$AGB_est_kg_02))/((dat$Days_03 - dat$Days_02)/365.25)

# June 2017 to June 2018
dat$AGR_AGB_2017_2018 <- (dat$AGB_est_kg_04 - dat$AGB_est_kg_03)/((dat$Days_04 - dat$Days_03)/365.25)
dat$RGR_AGB_2017_2018 <- (log(dat$AGB_est_kg_04) - log(dat$AGB_est_kg_03))/((dat$Days_04 - dat$Days_03)/365.25)

# June 2018 to July 2019
dat$AGR_AGB_2018_2019 <- (dat$AGB_est_kg_06 - dat$AGB_est_kg_04)/((dat$Days_06 - dat$Days_04)/365.25)
dat$RGR_AGB_2018_2019 <- (log(dat$AGB_est_kg_06) - log(dat$AGB_est_kg_04))/((dat$Days_06 - dat$Days_04)/365.25)

# June 2018 to January 2019
dat$AGR_AGB_2018_2019_01 <- (dat$AGB_est_kg_05 - dat$AGB_est_kg_04)/((dat$Days_05 - dat$Days_04)/365.25)
dat$RGR_AGB_2018_2019_01 <- (log(dat$AGB_est_kg_05) - log(dat$AGB_est_kg_04))/((dat$Days_05 - dat$Days_04)/365.25)

# January 2019 to July 2019
dat$AGR_AGB_2019_01_2019 <- (dat$AGB_est_kg_06 - dat$AGB_est_kg_05)/((dat$Days_06 - dat$Days_05)/365.25)
dat$RGR_AGB_2019_01_2019 <- (log(dat$AGB_est_kg_06) - log(dat$AGB_est_kg_05))/((dat$Days_06 - dat$Days_05)/365.25)

# 2016 to tf increments:

# AGB

# January 2016 to June 2018
dat$AGR_AGB_2016_2018 <- (dat$AGB_est_kg_04 - dat$AGB_est_kg_01)/((dat$Days_04 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_2018 <- (log(dat$AGB_est_kg_04) - log(dat$AGB_est_kg_01))/((dat$Days_04 - dat$Days_01)/365.25)

# January 2016 to January 2019
dat$AGR_AGB_2016_2019_01 <- (dat$AGB_est_kg_05 - dat$AGB_est_kg_01)/((dat$Days_05 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_2019_01 <- (log(dat$AGB_est_kg_05) - log(dat$AGB_est_kg_01))/((dat$Days_05 - dat$Days_01)/365.25)

# January 2016 to July 2019
dat$AGR_AGB_2016_2019 <- (dat$AGB_est_kg_06 - dat$AGB_est_kg_01)/((dat$Days_06 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_2019 <- (log(dat$AGB_est_kg_06) - log(dat$AGB_est_kg_01))/((dat$Days_06 - dat$Days_01)/365.25)

###############################################################
########################### Stats #############################
###############################################################

# Only use trees/measurements with Use_growth_0X >= 0.5
# Response variables:
# PSCA and GLSE and CAEQ ...
# AGB for ...
# 2017: 
# - AGR 2016-2017
# - RGR 2016-2017
# 2018: 
# - AGR 2016-2017, 2017-2018
# - RGR 2016-2017, 2017-2018
# - AGR 2016-2018
# - RGR 2016-2018
# 2019: 
# - AGR 2016-2017, 2017-2018, 2018-2019
# - RGR 2016-2017, 2017-2018, 2018-2019
# - AGR 2016-2019
# - RGR 2016-2019
#
# Driver: Treatment
# Covariates: For multi-year, AGB (fixed), tree (random)

########
# PSCA #
########

# Change this to set a different base
# MN as base
#dat[dat$Treatment == "MN",]$Treatment <- "A_MN"
#dat[dat$Treatment == "HN",]$Treatment <- "N_HN"
# Back to HN as base
#dat[dat$Treatment == "A_MN",]$Treatment <- "MN"
#dat[dat$Treatment == "N_HN",]$Treatment <- "HN"
# LN as base
#dat[dat$Treatment == "LN",]$Treatment <- "A_LN"
#dat[dat$Treatment == "MN",]$Treatment <- "B_MN"
# Back to HN as base
#dat[dat$Treatment == "A_LN",]$Treatment <- "LN"
#dat[dat$Treatment == "B_MN",]$Treatment <- "MN"

### 
# PSCA Aboveground Biomass
###

##### RGR

### 2016-2017

summary(lm_PSCA_AGB_RGR_2016_2017 <- lm(RGR_AGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_PSCA_AGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# b b a b

### 2016-2018

# Cumulative RGR
summary(lm_PSCA_AGB_RGR_2016_2018 <- lm(RGR_AGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_PSCA_AGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# b b a c

# Annual RGR
RGR_PSCA_AGB_678 <- 
	c(
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018
	)
Treatment_PSCA_678 <- 
	c(
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment
	)
AGB_PSCA_678 <- 
	c(
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03
	)
Tree_PSCA_678 <- 
	as.factor(c(
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID
	))
summary(lme_PSCA_AGB_RGR_201678 <- lme(RGR_PSCA_AGB_678 ~ 
	Treatment_PSCA_678 + AGB_PSCA_678,
	random=~1 | Tree_PSCA_678))
anova(lme_PSCA_AGB_RGR_201678)
emmeans(lme_PSCA_AGB_RGR_201678, list(pairwise ~ Treatment_PSCA_678), adjust = "tukey")
# ab b a b

### 2016-2019

# Cumulative RGR
summary(lm_PSCA_AGB_RGR_2016_2019 <- lm(RGR_AGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_PSCA_AGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# b b a b

# Annual RGR
RGR_PSCA_AGB_6789 <- 
	c(
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$RGR_AGB_2018_2019
	)
Treatment_PSCA_6789 <- 
	c(
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$Treatment
	)
AGB_PSCA_6789 <- 
	c(
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGB_est_kg_04
	)
Tree_PSCA_6789 <- 
	as.factor(c(
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID,
	dat[dat$Species=="PSCA" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$PID
	))
summary(lme_PSCA_AGB_RGR_2016789 <- lme(RGR_PSCA_AGB_6789 ~ 
	Treatment_PSCA_6789 + AGB_PSCA_6789,
	random=~1 | Tree_PSCA_6789))
anova(lme_PSCA_AGB_RGR_2016789)
emmeans(lme_PSCA_AGB_RGR_2016789, list(pairwise ~ Treatment_PSCA_6789), adjust = "tukey")
# ab ab a b

##### AGR

### 2016-2017

summary(lm_PSCA_AGB_AGR_2016_2017 <- lm(log(AGR_AGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_PSCA_AGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_PSCA_AGB_AGR_2016_2018 <- lm(log(AGR_AGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_PSCA_AGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSCA_AGB_678 <- 
	c(
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018
	)
summary(lme_PSCA_AGB_AGR_201678 <- lme(log(AGR_PSCA_AGB_678) ~ 
	Treatment_PSCA_678 + AGB_PSCA_678,
	random=~1 | Tree_PSCA_678))
anova(lme_PSCA_AGB_AGR_201678)
emmeans(lme_PSCA_AGB_AGR_201678, list(pairwise ~ Treatment_PSCA_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_PSCA_AGB_AGR_2016_2019 <- lm(log(AGR_AGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_PSCA_AGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_PSCA_AGB_6789 <- 
	c(
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018,
	dat[dat$Species=="PSCA" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGR_AGB_2018_2019
	)
summary(lme_PSCA_AGB_AGR_2016789 <- lme(log(AGR_PSCA_AGB_6789) ~ 
	Treatment_PSCA_6789 + AGB_PSCA_6789,
	random=~1 | Tree_PSCA_6789))
anova(lme_PSCA_AGB_AGR_2016789)
emmeans(lme_PSCA_AGB_AGR_2016789, list(pairwise ~ Treatment_PSCA_6789), adjust = "tukey")
# a a a a



########
# GLSE #
########

# Change this to set a different base
# MN as base
#dat[dat$Treatment == "MN",]$Treatment <- "A_MN"
#dat[dat$Treatment == "HN",]$Treatment <- "N_HN"
# Back to HN as base
#dat[dat$Treatment == "A_MN",]$Treatment <- "MN"
#dat[dat$Treatment == "N_HN",]$Treatment <- "HN"
# LN as base
#dat[dat$Treatment == "LN",]$Treatment <- "A_LN"
#dat[dat$Treatment == "MN",]$Treatment <- "B_MN"
# Back to HN as base
#dat[dat$Treatment == "A_LN",]$Treatment <- "LN"
#dat[dat$Treatment == "B_MN",]$Treatment <- "MN"

### 
# GLSE Aboveground Biomass
###

##### RGR

### 2016-2017

summary(lm_GLSE_AGB_RGR_2016_2017 <- lm(RGR_AGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_GLSE_AGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_GLSE_AGB_RGR_2016_2018 <- lm(RGR_AGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_GLSE_AGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_GLSE_AGB_678 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018
	)
Treatment_GLSE_678 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment
	)
AGB_GLSE_678 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03
	)
Tree_GLSE_678 <- 
	as.factor(c(
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID
	))
summary(lme_GLSE_AGB_RGR_201678 <- lme(RGR_GLSE_AGB_678 ~ 
	Treatment_GLSE_678 + AGB_GLSE_678,
	random=~1 | Tree_GLSE_678))
anova(lme_GLSE_AGB_RGR_201678)
emmeans(lme_GLSE_AGB_RGR_201678, list(pairwise ~ Treatment_GLSE_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_GLSE_AGB_RGR_2016_2019 <- lm(RGR_AGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_GLSE_AGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_GLSE_AGB_6789 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$RGR_AGB_2018_2019
	)
Treatment_GLSE_6789 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$Treatment
	)
AGB_GLSE_6789 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGB_est_kg_04
	)
Tree_GLSE_6789 <- 
	as.factor(c(
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID,
	dat[dat$Species=="GLSE" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$PID
	))
summary(lme_GLSE_AGB_RGR_2016789 <- lme(RGR_GLSE_AGB_6789 ~ 
	Treatment_GLSE_6789 + AGB_GLSE_6789,
	random=~1 | Tree_GLSE_6789))
anova(lme_GLSE_AGB_RGR_2016789)
emmeans(lme_GLSE_AGB_RGR_2016789, list(pairwise ~ Treatment_GLSE_6789), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_GLSE_AGB_AGR_2016_2017 <- lm(log(AGR_AGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_GLSE_AGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_GLSE_AGB_AGR_2016_2018 <- lm(log(AGR_AGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_GLSE_AGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_GLSE_AGB_678 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018
	)
summary(lme_GLSE_AGB_AGR_201678 <- lme(log(AGR_GLSE_AGB_678) ~ 
	Treatment_GLSE_678 + AGB_GLSE_678,
	random=~1 | Tree_GLSE_678))
anova(lme_GLSE_AGB_AGR_201678)
emmeans(lme_GLSE_AGB_AGR_201678, list(pairwise ~ Treatment_GLSE_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_GLSE_AGB_AGR_2016_2019 <- lm(log(AGR_AGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_GLSE_AGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_GLSE_AGB_6789 <- 
	c(
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018,
	dat[dat$Species=="GLSE" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGR_AGB_2018_2019
	)
summary(lme_GLSE_AGB_AGR_2016789 <- lme(log(AGR_GLSE_AGB_6789) ~ 
	Treatment_GLSE_6789 + AGB_GLSE_6789,
	random=~1 | Tree_GLSE_6789))
anova(lme_GLSE_AGB_AGR_2016789)
emmeans(lme_GLSE_AGB_AGR_2016789, list(pairwise ~ Treatment_GLSE_6789), adjust = "tukey")
# a a a a



########
# CAEQ #
########

# Change this to set a different base
# MN as base
#dat[dat$Treatment == "MN",]$Treatment <- "A_MN"
#dat[dat$Treatment == "HN",]$Treatment <- "N_HN"
# Back to HN as base
#dat[dat$Treatment == "A_MN",]$Treatment <- "MN"
#dat[dat$Treatment == "N_HN",]$Treatment <- "HN"
# LN as base
#dat[dat$Treatment == "LN",]$Treatment <- "A_LN"
#dat[dat$Treatment == "MN",]$Treatment <- "B_MN"
# Back to HN as base
#dat[dat$Treatment == "A_LN",]$Treatment <- "LN"
#dat[dat$Treatment == "B_MN",]$Treatment <- "MN"

### 
# CAEQ Aboveground Biomass
###

##### RGR

### 2016-2017

summary(lm_CAEQ_AGB_RGR_2016_2017 <- lm(RGR_AGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_CAEQ_AGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_CAEQ_AGB_RGR_2016_2018 <- lm(RGR_AGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_CAEQ_AGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_CAEQ_AGB_678 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018
	)
Treatment_CAEQ_678 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment
	)
AGB_CAEQ_678 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03
	)
Tree_CAEQ_678 <- 
	as.factor(c(
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID
	))
summary(lme_CAEQ_AGB_RGR_201678 <- lme(RGR_CAEQ_AGB_678 ~ 
	Treatment_CAEQ_678 + AGB_CAEQ_678,
	random=~1 | Tree_CAEQ_678))
anova(lme_CAEQ_AGB_RGR_201678)
emmeans(lme_CAEQ_AGB_RGR_201678, list(pairwise ~ Treatment_CAEQ_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_CAEQ_AGB_RGR_2016_2019 <- lm(RGR_AGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_CAEQ_AGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_CAEQ_AGB_6789 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$RGR_AGB_2018_2019
	)
Treatment_CAEQ_6789 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$Treatment
	)
AGB_CAEQ_6789 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGB_est_kg_04
	)
Tree_CAEQ_6789 <- 
	as.factor(c(
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID,
	dat[dat$Species=="CAEQ" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$PID
	))
summary(lme_CAEQ_AGB_RGR_2016789 <- lme(RGR_CAEQ_AGB_6789 ~ 
	Treatment_CAEQ_6789 + AGB_CAEQ_6789,
	random=~1 | Tree_CAEQ_6789))
anova(lme_CAEQ_AGB_RGR_2016789)
emmeans(lme_CAEQ_AGB_RGR_2016789, list(pairwise ~ Treatment_CAEQ_6789), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_CAEQ_AGB_AGR_2016_2017 <- lm(log(AGR_AGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_CAEQ_AGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_CAEQ_AGB_AGR_2016_2018 <- lm(log(AGR_AGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_CAEQ_AGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_CAEQ_AGB_678 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018
	)
summary(lme_CAEQ_AGB_AGR_201678 <- lme(log(AGR_CAEQ_AGB_678) ~ 
	Treatment_CAEQ_678 + AGB_CAEQ_678,
	random=~1 | Tree_CAEQ_678))
anova(lme_CAEQ_AGB_AGR_201678)
emmeans(lme_CAEQ_AGB_AGR_201678, list(pairwise ~ Treatment_CAEQ_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_CAEQ_AGB_AGR_2016_2019 <- lm(log(AGR_AGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_CAEQ_AGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_CAEQ_AGB_6789 <- 
	c(
	dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018,
	dat[dat$Species=="CAEQ" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGR_AGB_2018_2019
	)
summary(lme_CAEQ_AGB_AGR_2016789 <- lme(log(AGR_CAEQ_AGB_6789) ~ 
	Treatment_CAEQ_6789 + AGB_CAEQ_6789,
	random=~1 | Tree_CAEQ_6789))
anova(lme_CAEQ_AGB_AGR_2016789)
emmeans(lme_CAEQ_AGB_AGR_2016789, list(pairwise ~ Treatment_CAEQ_6789), adjust = "tukey")
# a a a a




#############################################################################
#############################################################################
#############################################################################