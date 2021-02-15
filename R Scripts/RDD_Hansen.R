library(learnr)
library(haven)
library(tidyverse)
library(stargazer)
library(estimatr)
library(rdd)
library(rdrobust)
library(rddensity)
library(ggplot2)

dwi <- read.csv('https://raw.githubusercontent.com/hannahjonesut/causal-inference-class/master/Data/hansen_dwi.csv')

head(dwi)

#create a dummy variable for BAC >= 0.08
dwi <- dwi %>%
  mutate( 
    bac_high  = if_else(bac1 >= 0.08, 1, 0))

#sort on running variable to see if any manipulation
d1 = dwi %>%
  group_by(bac1)%>%
  summarize(sum_recid = sum(recidivism))

ggplot(data = d1) +
  geom_col(aes(x = bac1, y = sum_recid)) +
  labs(
    title = "Frequency of Recidivism vs Blood Alcohol Content Histogram",
    caption = "Based on administrative records from the Washington State Impaired Driver Testing Program, 1999–2007.Vertical line at BAC = 0.08.",
    x="Blood Alcohol Content",
    y="Frequency") +
  geom_vline(xintercept = .08, colour = "red", linetype = 1) 


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



