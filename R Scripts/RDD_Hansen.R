library(learnr)
library(haven)
library(tidyverse)
library(stargazer)
library(estimatr)
library(rdd)
library(rdrobust)
library(rddensity)
library(ggplot2)
library(dplyr)

dwi <- read.csv('https://raw.githubusercontent.com/hannahjonesut/causal-inference-class/master/Data/hansen_dwi.csv')

head(dwi)

#create a dummy variable for BAC >= 0.08
dwi <- dwi %>%
  mutate( 
    bac_high  = if_else(bac1 >= 0.08, 1, 0))

#sort on running variable to see if any manipulation
bac1_recid = dwi %>%
  group_by(bac1)%>%
  summarize(sum_recid = sum(recidivism))

ggplot(data = bac1_recid) +
  geom_col(aes(x = bac1, y = sum_recid)) +
  labs(
    title = "Frequency of Recidivism vs Blood Alcohol Content Histogram",
    caption = "Based on administrative records from the Washington State Impaired Driver Testing Program, 1999–2007.Vertical line at BAC = 0.08.",
    x="Blood Alcohol Content",
    y="Frequency") +
  geom_vline(xintercept = .08, colour = "red", linetype = 1) 

#density test
rdd<- rddensity(dwi$bac1, c = 0.08)
summary(rdd)

#Recreate Table 2 Panel A but only white, male, age and accident (acc) as dependent variables

lin_regr1 = lm_robust(bac_high ~ male+white+aged+acc, data = dwi)
lin_regr1

lin_regr2 = lm_robust(recidivism ~ male+white+aged+acc+bac_high+bac1+bac1*bac_high, data = dwi)
lin_regr2

#lm_malebac= lm_robust(male~bac_high, data = dwi)
#lm_whitebac= lm_robust(white~bac_high, data = dwi)
#lm_accbac= lm_robust(acc~bac_high, data = dwi)
#lm_agedbac= lm_robust(aged~bac_high, data = dwi)


#Recreate figure 2 **RIGHT ONE**
dwi_male = dwi %>%
  group_by(bac1)%>%
  summarize(n=n(), male_prop = sum(male)/n)
dwi_male

dwi_male_low <- dwi_male %>%
  filter(bac1 < 0.08)
lm_male_low = lm_robust(male_prop~bac1, data=dwi_male_low)
lm_male_low 

dwi_male_high <- dwi_male %>%
  filter(bac_high >= 0.08)
lm_male_high = lm_robust(male_prop~bac1, data=dwi_male_high)
lm_male_high

ggplot(data = dwi_male)+
  geom_point(aes(x=bac1, y= male_prop), alpha = 0.25)+
  geom_abline(aes(x = bac1, y = lm_male_low))+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)

dwi_white = dwi %>%
  group_by(bac1)%>%
  summarize(n=n(), white_prop = sum(white)/n)
dwi_white

dwi_acc = dwi %>%
  group_by(bac1)%>%
  summarize(n=n(), acc_prop = sum(acc)/n)
dwi_acc

dwi_aged = dwi %>%
  group_by(bac1)%>%
  summarize(n=n(), aged_prop = sum(aged)/n)
dwi_aged



ggplot(data = dwi_white)+
  geom_point(aes(x=bac1, y= white_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)

ggplot(data = dwi_acc)+
  geom_point(aes(x=bac1, y= acc_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)

ggplot(data = dwi_aged)+
  geom_point(aes(x=bac1, y= aged_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)


dat <- tibble(
  x = rnorm(1000, 50, 25)
) %>%
  mutate(
    x = if_else(x < 0, 0, x)
  ) %>%
  filter(x < 100)
# cutoff at x = 50
dat <- dat %>% 
  mutate(
    D  = if_else(x > 50, 1, 0),
    y1 = 25 + 0 * D + 1.5 * x + rnorm(n(), 0, 20)
  )
cli::cli_text("Counterfactual Potential Outcomes")
ggplot(aes(x, y1, colour = factor(D)), data = dat) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 50, colour = "grey", linetype = 2)+
  stat_smooth(method = "lm", se = F) +
  labs(x = "Test score (X)", y = "Potential Outcome (Y1)")



