# NSF FX data analysis for the Volcano biomass data

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
# "HI_V_FX_Size_FoliarCNIsotope_Data.csv" to reflect these odd datapoints
# based on biomass trends.
# 1 means not odd in any way.
# 0.5 means use July 2019 if not using January 2019, but not if using January 2019.
# 0 means sufficiently ill or damaged that it shouldn't be used.

rm(list=ls())
library(nlme)
library(emmeans)

setwd("/Users/duncanmenge/Documents/Academia/Grants/NSF/2014_Strategies/Data/Master datafiles/")

dat <- read.csv("HI_V_FX_Size_FoliarCNIsotope_Data.csv")[1:96,1:111]

# Dates
# _01: January 2016
# _02: July 2016
# _03: June 2017
# _04: June 2018
# _05: January 2019
# _06: July 2019

# Option to print out how many are in each category (uncomment)

#print(dat[dat$Species=="ACKO" & dat$Treatment=="LN",c(22,33,44,59,76,87)])
#print(dat[dat$Species=="ACKO" & dat$Treatment=="MN",c(22,33,44,59,76,87)])
#print(dat[dat$Species=="ACKO" & dat$Treatment=="HN",c(22,33,44,59,76,87)])
#print(dat[dat$Species=="ACKO" & dat$Treatment=="PHN",c(22,33,44,59,76,87)])

#print(dat[dat$Species=="MOFA" & dat$Treatment=="LN",c(22,33,44,59,76,87)])
#print(dat[dat$Species=="MOFA" & dat$Treatment=="MN",c(22,33,44,59,76,87)])
#print(dat[dat$Species=="MOFA" & dat$Treatment=="HN",c(22,33,44,59,76,87)])
#print(dat[dat$Species=="MOFA" & dat$Treatment=="PHN",c(22,33,44,59,76,87)])

#print(dat[dat$Species=="DOVI" & dat$Treatment=="LN",c(22,33,44,59,76,87)])
#print(dat[dat$Species=="DOVI" & dat$Treatment=="MN",c(22,33,44,59,76,87)])
#print(dat[dat$Species=="DOVI" & dat$Treatment=="HN",c(22,33,44,59,76,87)])
#print(dat[dat$Species=="DOVI" & dat$Treatment=="PHN",c(22,33,44,59,76,87)])

#############################################################
################ Calculate growth rates #####################
#############################################################

# "AGR" is absolute growth rate and "RGR" is relative growth rate. 
# Each is calculated for aboveground biomass, 
# for the following time periods: 

# Since some of the trees were replanted in July 2016 or January 2017,
# I need to account for that.

# Annual/semiannual increments:

# AGB

# July 2016 to June 2017
dat$AGR_AGB_2016_2017 <- (dat$AGB_est_kg_03 - dat$AGB_est_kg_02)/((dat$Days_03 - dat$Days_02)/365.25)
dat$RGR_AGB_2016_2017 <- (log(dat$AGB_est_kg_03) - log(dat$AGB_est_kg_02))/((dat$Days_03 - dat$Days_02)/365.25)

replant <- which(dat$Use_growth_02==0 & dat$Use_growth_01==1)
dat[replant,]$AGR_AGB_2016_2017 <-
	(dat[replant,]$AGB_est_kg_03 - dat[replant,]$AGB_est_kg_01)/
	((dat[replant,]$Days_03 - dat[replant,]$Days_01)/365.25)
dat[replant,]$RGR_AGB_2016_2017 <-
	(log(dat[replant,]$AGB_est_kg_03) - log(dat[replant,]$AGB_est_kg_01))/
	((dat[replant,]$Days_03 - dat[replant,]$Days_01)/365.25)

# June 2017 to June 2018
dat$AGR_AGB_2017_2018 <- (dat$AGB_est_kg_04 - dat$AGB_est_kg_03)/((dat$Days_04 - dat$Days_03)/365.25)
dat$RGR_AGB_2017_2018 <- (log(dat$AGB_est_kg_04) - log(dat$AGB_est_kg_03))/((dat$Days_04 - dat$Days_03)/365.25)

# June 2018 to July 2019
dat$AGR_AGB_2018_2019 <- (dat$AGB_est_kg_06 - dat$AGB_est_kg_04)/((dat$Days_06 - dat$Days_04)/365.25)
dat$RGR_AGB_2018_2019 <- (log(dat$AGB_est_kg_06) - log(dat$AGB_est_kg_04))/((dat$Days_06 - dat$Days_04)/365.25)

# 2016 to tf increments:

# AGB

# January 2016 to June 2018
dat$AGR_AGB_2016_2018 <- (dat$AGB_est_kg_04 - dat$AGB_est_kg_01)/((dat$Days_04 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_2018 <- (log(dat$AGB_est_kg_04) - log(dat$AGB_est_kg_01))/((dat$Days_04 - dat$Days_01)/365.25)

# Replant defined the opposite way here
replant <- which(dat$Use_growth_01==0 & dat$Use_growth_02==1)
dat[replant,]$AGR_AGB_2016_2018 <- 
	(dat[replant,]$AGB_est_kg_04 - dat[replant,]$AGB_est_kg_02)/
	((dat[replant,]$Days_04 - dat[replant,]$Days_02)/365.25)
dat[replant,]$RGR_AGB_2016_2018 <- 
	(log(dat[replant,]$AGB_est_kg_04) - log(dat[replant,]$AGB_est_kg_02))/
	((dat[replant,]$Days_04 - dat[replant,]$Days_02)/365.25)

# January 2016 to July 2019
dat$AGR_AGB_2016_2019 <- (dat$AGB_est_kg_06 - dat$AGB_est_kg_01)/((dat$Days_06 - dat$Days_01)/365.25)
dat$RGR_AGB_2016_2019 <- (log(dat$AGB_est_kg_06) - log(dat$AGB_est_kg_01))/((dat$Days_06 - dat$Days_01)/365.25)

dat[replant,]$AGR_AGB_2016_2019 <- 
	(dat[replant,]$AGB_est_kg_06 - dat[replant,]$AGB_est_kg_02)/
	((dat[replant,]$Days_06 - dat[replant,]$Days_02)/365.25)
dat[replant,]$RGR_AGB_2016_2019 <- 
	(log(dat[replant,]$AGB_est_kg_06) - log(dat[replant,]$AGB_est_kg_02))/
	((dat[replant,]$Days_06 - dat[replant,]$Days_02)/365.25)

###############################################################
########################### Stats #############################
###############################################################

# Only use trees/measurements with Use_growth_0X >= 0.5
# Response variables:
# DOVI and ACKO and MOFA ...
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
# DOVI #
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
# DOVI Aboveground Biomass
###

##### RGR

### 2016-2017

summary(lm_DOVI_AGB_RGR_2016_2017 <- lm(RGR_AGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_DOVI_AGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_DOVI_AGB_RGR_2016_2018 <- lm(RGR_AGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_DOVI_AGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_DOVI_AGB_678 <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018
	)
Treatment_DOVI_678 <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment
	)
AGB_DOVI_678 <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03
	)
AGB_DOVI_678[is.na(AGB_DOVI_678)] <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,][is.na(AGB_DOVI_678),]$AGB_est_kg_02
	)
Tree_DOVI_678 <- 
	as.factor(c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID
	))
summary(lme_DOVI_AGB_RGR_201678 <- lme(RGR_DOVI_AGB_678 ~ 
	Treatment_DOVI_678 + AGB_DOVI_678,
	random=~1 | Tree_DOVI_678))
anova(lme_DOVI_AGB_RGR_201678)
emmeans(lme_DOVI_AGB_RGR_201678, list(pairwise ~ Treatment_DOVI_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_DOVI_AGB_RGR_2016_2019 <- lm(RGR_AGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_DOVI_AGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_DOVI_AGB_6789 <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$RGR_AGB_2018_2019
	)
Treatment_DOVI_6789 <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$Treatment
	)
AGB_DOVI_6789 <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGB_est_kg_04
	)
AGB_DOVI_6789[is.na(AGB_DOVI_6789)] <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,][is.na(AGB_DOVI_6789),]$AGB_est_kg_02
	)
Tree_DOVI_6789 <- 
	as.factor(c(
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID,
	dat[dat$Species=="DOVI" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$PID
	))
summary(lme_DOVI_AGB_RGR_2016789 <- lme(RGR_DOVI_AGB_6789 ~ 
	Treatment_DOVI_6789 + AGB_DOVI_6789,
	random=~1 | Tree_DOVI_6789))
anova(lme_DOVI_AGB_RGR_2016789)
emmeans(lme_DOVI_AGB_RGR_2016789, list(pairwise ~ Treatment_DOVI_6789), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_DOVI_AGB_AGR_2016_2017 <- lm(log(AGR_AGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_DOVI_AGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_DOVI_AGB_AGR_2016_2018 <- lm(log(AGR_AGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_DOVI_AGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_DOVI_AGB_678 <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018
	)
summary(lme_DOVI_AGB_AGR_201678 <- lme(log(AGR_DOVI_AGB_678) ~ 
	Treatment_DOVI_678 + AGB_DOVI_678,
	random=~1 | Tree_DOVI_678))
anova(lme_DOVI_AGB_AGR_201678)
emmeans(lme_DOVI_AGB_AGR_201678, list(pairwise ~ Treatment_DOVI_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_DOVI_AGB_AGR_2016_2019 <- lm(log(AGR_AGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_DOVI_AGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_DOVI_AGB_6789 <- 
	c(
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018,
	dat[dat$Species=="DOVI" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGR_AGB_2018_2019
	)
summary(lme_DOVI_AGB_AGR_2016789 <- lme(log(AGR_DOVI_AGB_6789) ~ 
	Treatment_DOVI_6789 + AGB_DOVI_6789,
	random=~1 | Tree_DOVI_6789))
anova(lme_DOVI_AGB_AGR_2016789)
emmeans(lme_DOVI_AGB_AGR_2016789, list(pairwise ~ Treatment_DOVI_6789), adjust = "tukey")
# a a a a



########
# ACKO #
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
# ACKO Aboveground Biomass
###

##### RGR

### 2016-2017

summary(lm_ACKO_AGB_RGR_2016_2017 <- lm(RGR_AGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_ACKO_AGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative RGR
summary(lm_ACKO_AGB_RGR_2016_2018 <- lm(RGR_AGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_ACKO_AGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ACKO_AGB_678 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018
	)
Treatment_ACKO_678 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment
	)
AGB_ACKO_678 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03
	)
AGB_ACKO_678[is.na(AGB_ACKO_678)] <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,][is.na(AGB_ACKO_678),]$AGB_est_kg_02
	)
Tree_ACKO_678 <- 
	as.factor(c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID
	))
summary(lme_ACKO_AGB_RGR_201678 <- lme(RGR_ACKO_AGB_678 ~ 
	Treatment_ACKO_678 + AGB_ACKO_678,
	random=~1 | Tree_ACKO_678))
anova(lme_ACKO_AGB_RGR_201678)
emmeans(lme_ACKO_AGB_RGR_201678, list(pairwise ~ Treatment_ACKO_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative RGR
summary(lm_ACKO_AGB_RGR_2016_2019 <- lm(RGR_AGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_ACKO_AGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual RGR
RGR_ACKO_AGB_6789 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$RGR_AGB_2018_2019
	)
Treatment_ACKO_6789 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$Treatment
	)
AGB_ACKO_6789 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGB_est_kg_04
	)
AGB_ACKO_6789[is.na(AGB_ACKO_6789)] <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,][is.na(AGB_ACKO_6789),]$AGB_est_kg_02
	)
Tree_ACKO_6789 <- 
	as.factor(c(
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID,
	dat[dat$Species=="ACKO" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$PID
	))
summary(lme_ACKO_AGB_RGR_2016789 <- lme(RGR_ACKO_AGB_6789 ~ 
	Treatment_ACKO_6789 + AGB_ACKO_6789,
	random=~1 | Tree_ACKO_6789))
anova(lme_ACKO_AGB_RGR_2016789)
emmeans(lme_ACKO_AGB_RGR_2016789, list(pairwise ~ Treatment_ACKO_6789), adjust = "tukey")
# a a a a

##### AGR

### 2016-2017

summary(lm_ACKO_AGB_AGR_2016_2017 <- lm(log(AGR_AGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_ACKO_AGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_ACKO_AGB_AGR_2016_2018 <- lm(log(AGR_AGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_ACKO_AGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ACKO_AGB_678 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018
	)
summary(lme_ACKO_AGB_AGR_201678 <- lme(log(AGR_ACKO_AGB_678) ~ 
	Treatment_ACKO_678 + AGB_ACKO_678,
	random=~1 | Tree_ACKO_678))
anova(lme_ACKO_AGB_AGR_201678)
emmeans(lme_ACKO_AGB_AGR_201678, list(pairwise ~ Treatment_ACKO_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_ACKO_AGB_AGR_2016_2019 <- lm(log(AGR_AGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_ACKO_AGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_ACKO_AGB_6789 <- 
	c(
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018,
	dat[dat$Species=="ACKO" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGR_AGB_2018_2019
	)
summary(lme_ACKO_AGB_AGR_2016789 <- lme(log(AGR_ACKO_AGB_6789) ~ 
	Treatment_ACKO_6789 + AGB_ACKO_6789,
	random=~1 | Tree_ACKO_6789))
anova(lme_ACKO_AGB_AGR_2016789)
emmeans(lme_ACKO_AGB_AGR_2016789, list(pairwise ~ Treatment_ACKO_6789), adjust = "tukey")
# a a a a



########
# MOFA #
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
# MOFA Aboveground Biomass
###

##### RGR

### 2016-2017

summary(lm_MOFA_AGB_RGR_2016_2017 <- lm(RGR_AGB_2016_2017 ~ 
	Treatment,
	data=dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_MOFA_AGB_RGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# ab ab b a

### 2016-2018

# Cumulative RGR
summary(lm_MOFA_AGB_RGR_2016_2018 <- lm(RGR_AGB_2016_2018 ~ 
	Treatment,
	data=dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_MOFA_AGB_RGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a ab b a

# Annual RGR
RGR_MOFA_AGB_678 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018
	)
Treatment_MOFA_678 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment
	)
AGB_MOFA_678 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03
	)
AGB_MOFA_678[is.na(AGB_MOFA_678)] <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,][is.na(AGB_MOFA_678),]$AGB_est_kg_02
	)
Tree_MOFA_678 <- 
	as.factor(c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID
	))
summary(lme_MOFA_AGB_RGR_201678 <- lme(RGR_MOFA_AGB_678 ~ 
	Treatment_MOFA_678 + AGB_MOFA_678,
	random=~1 | Tree_MOFA_678))
anova(lme_MOFA_AGB_RGR_201678)
emmeans(lme_MOFA_AGB_RGR_201678, list(pairwise ~ Treatment_MOFA_678), adjust = "tukey")
# ab ab b a

### 2016-2019

# Cumulative RGR
summary(lm_MOFA_AGB_RGR_2016_2019 <- lm(RGR_AGB_2016_2019 ~ 
	Treatment,
	data=dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_MOFA_AGB_RGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a ab b a

# Annual RGR
RGR_MOFA_AGB_6789 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$RGR_AGB_2016_2017,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$RGR_AGB_2017_2018,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$RGR_AGB_2018_2019
	)
Treatment_MOFA_6789 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$Treatment,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$Treatment,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$Treatment
	)
AGB_MOFA_6789 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGB_est_kg_01,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGB_est_kg_03,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGB_est_kg_04
	)
AGB_MOFA_6789[is.na(AGB_MOFA_6789)] <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,][is.na(AGB_MOFA_6789),]$AGB_est_kg_02
	)
Tree_MOFA_6789 <- 
	as.factor(c(
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$PID,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$PID,
	dat[dat$Species=="MOFA" & !is.na(dat$RGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$PID
	))
summary(lme_MOFA_AGB_RGR_2016789 <- lme(RGR_MOFA_AGB_6789 ~ 
	Treatment_MOFA_6789 + AGB_MOFA_6789,
	random=~1 | Tree_MOFA_6789))
anova(lme_MOFA_AGB_RGR_2016789)
emmeans(lme_MOFA_AGB_RGR_2016789, list(pairwise ~ Treatment_MOFA_6789), adjust = "tukey")
# a ab b a

##### AGR

### 2016-2017

summary(lm_MOFA_AGB_AGR_2016_2017 <- lm(log(AGR_AGB_2016_2017) ~ 
	Treatment,
	data=dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2016_2017) & 
	dat$Use_growth_03 >= 0.5,]))
emmeans(lm_MOFA_AGB_AGR_2016_2017, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

### 2016-2018

# Cumulative AGR
summary(lm_MOFA_AGB_AGR_2016_2018 <- lm(log(AGR_AGB_2016_2018) ~ 
	Treatment,
	data=dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2016_2018) & 
	dat$Use_growth_04 >= 0.5,]))
emmeans(lm_MOFA_AGB_AGR_2016_2018, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_MOFA_AGB_678 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018
	)
summary(lme_MOFA_AGB_AGR_201678 <- lme(log(AGR_MOFA_AGB_678) ~ 
	Treatment_MOFA_678 + AGB_MOFA_678,
	random=~1 | Tree_MOFA_678))
anova(lme_MOFA_AGB_AGR_201678)
emmeans(lme_MOFA_AGB_AGR_201678, list(pairwise ~ Treatment_MOFA_678), adjust = "tukey")
# a a a a

### 2016-2019

# Cumulative AGR
summary(lm_MOFA_AGB_AGR_2016_2019 <- lm(log(AGR_AGB_2016_2019) ~ 
	Treatment,
	data=dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2016_2019) & 
	dat$Use_growth_06 >= 0.5,]))
emmeans(lm_MOFA_AGB_AGR_2016_2019, list(pairwise ~ Treatment), adjust = "tukey")
# a a a a

# Annual AGR
AGR_MOFA_AGB_6789 <- 
	c(
	dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2016_2017) & 
		dat$Use_growth_03 >= 0.5,]$AGR_AGB_2016_2017,
	dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2017_2018) & 
		dat$Use_growth_04 >= 0.5,]$AGR_AGB_2017_2018,
	dat[dat$Species=="MOFA" & !is.na(dat$AGR_AGB_2018_2019) & 
		dat$Use_growth_06 >= 0.5,]$AGR_AGB_2018_2019
	)
summary(lme_MOFA_AGB_AGR_2016789 <- lme(log(AGR_MOFA_AGB_6789) ~ 
	Treatment_MOFA_6789 + AGB_MOFA_6789,
	random=~1 | Tree_MOFA_6789))
anova(lme_MOFA_AGB_AGR_2016789)
emmeans(lme_MOFA_AGB_AGR_2016789, list(pairwise ~ Treatment_MOFA_6789), adjust = "tukey")
# a a a a




#############################################################################
#############################################################################
#############################################################################