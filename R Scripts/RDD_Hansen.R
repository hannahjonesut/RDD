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
library(kableExtra)
library(ggpubr)

dwi <- read.csv('https://raw.githubusercontent.com/hannahjonesut/causal-inference-class/master/Data/hansen_dwi.csv')

#Number 3-- create a dummy variable for BAC >= 0.08
dwi <- dwi %>%
  mutate( 
    bac_high= if_else(bac1 >= 0.08, 1, 0))

#Number 4-- sort on running variable to see if any manipulation
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

#No 5 use RDD

cov_male <- RDestimate(male~bac1 | aged+acc+white+bac_high+bac1*bac_high, 
                       data = dwi,cutpoint = 0.08, bw = 0.05, kernel = "rectangular")

cov_age<- RDestimate(aged~bac1 | male+acc+white+bac_high+bac1*bac_high, 
                     data = dwi,cutpoint = 0.08, bw = 0.05, kernel = "rectangular")

cov_acc<- RDestimate(acc~bac1 | aged+male+white+bac_high+bac1*bac_high, 
                     data = dwi,cutpoint = 0.08, bw = 0.05, kernel = "rectangular")

cov_white<- RDestimate(white~bac1 | aged+acc+male+bac_high+bac1*bac_high, 
                       data = dwi,cutpoint = 0.08, bw = 0.05, kernel = "rectangular")

cov_bal<-data.frame(
  Estimates = c("Coefficient", "SE", "Z", "P-Value", "Confidence Interval"),
  malecov = c(cov_male$est[1], cov_male$se[1], cov_male$z[1], cov_male$p[1], cov_male$ci[1]),
  agecov = c(cov_age$est[1], cov_age$se[1], cov_age$z[1], cov_age$p[1], cov_age$ci[1]),
  acccov = c(cov_acc$est[1], cov_acc$se[1], cov_acc$z[1], cov_acc$p[1], cov_acc$ci[1]),
  whitecov = c(cov_white$est[1], cov_white$se[1], cov_white$z[1], cov_white$p[1], cov_white$ci[1])
)
  
cov_bal %>%
  kable(
    col.names = c("Characteristic", "Male", "Age", "Accident", "White"),
    digits = 3,
    caption = "Panel A Covariate Balance"
  )%>%
  kable_classic(full_width = F, html_font = "Cambria")
  


#Number 6

#Recreate figure 2
dwi_cov = dwi %>%
  group_by(bac1)%>%
  summarize(n=n(), male_prop = sum(male)/n, white_prop = sum(white)/n, acc_prop = sum(acc)/n, aged_prop = sum(aged)/n)

#MODEL setup
X <- dwi$bac1 - .08 # center the data
Y_male <- dwi$male

X2 <- (dwi$bac1 - .08)^2 

Xl <- (X < 0) * X # X if X < 0, o.w. 0
Xl2 <- Xl^2
Xr <- (X >= 0) * X # X if X >= 0, o.w. 0
Xr2 <- Xr^2
Tr <- as.integer(X >= 0) # dummy variable


weights <- rdd::kernelwts(X, center = 0, bw = .05, kernel = "rectangular")

#male
model_male <- lm(Y_male ~ Tr + Xl + Xr, weights = weights)

male_df <- data.frame(male_pred = predict(model_male, dwi), bac1 = dwi$bac1)
male_CI <- data.frame(male_CI = predict(model_male, dwi, interval = "confidence"),  bac1 = dwi$bac1)

male_graph <- ggplot(dwi_cov, aes(x = bac1))+
  geom_point(aes(y= male_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = male_df, aes(y = male_pred))+
  geom_line(data= male_CI, aes(y = male_CI.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= male_CI, aes(y = male_CI.lwr), color = "steelblue", linetype = 3 )+
  xlim(0, 0.2)+
  ylim(0.6,0.9)
male_graph 

quad_male <- lm(Y_male ~ Tr + Xl + Xl2 + Xr + Xr2, weights = weights)

male_df_q <- data.frame(male_pred_q = predict(quad_male, dwi), bac1 = dwi$bac1)
male_CI_q <- data.frame(male_CI_q = predict(quad_male, dwi, interval = "confidence"),  bac1 = dwi$bac1)

male_graph_q <- ggplot(dwi_cov, aes(x = bac1))+
  geom_point(aes(y= male_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = male_df_q, aes(y = male_pred_q))+
  geom_line(data= male_CI_q, aes(y = male_CI_q.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= male_CI_q, aes(y = male_CI_q.lwr), color = "steelblue", linetype = 3 )+
  xlim(0, 0.2)+
  ylim(0.6,0.9)
male_graph_q

#white
Y_white <- dwi$white

model_white <- lm(Y_white ~ Tr + Xl + Xr, weights = weights)
white_CI <- data.frame(white_CI = predict(model_white, dwi, interval = "confidence"),  bac1 = dwi$bac1)

white_df <- data.frame(white_pred = predict(model_white, dwi), bac1 = dwi$bac1)
white_graph <- ggplot(dwi_cov, aes(x = bac1))+
  geom_point(aes(y= white_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = white_df, aes(y = white_pred))+
  geom_line(data= white_CI, aes(y = white_CI.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= white_CI, aes(y = white_CI.lwr), color = "steelblue", linetype = 3 )+
  xlim(0, 0.2)+
  ylim(0.75,0.9)
white_graph 

#white quadratic

quad_white <- lm(Y_white ~ Tr + Xl + Xl2 + Xr + Xr2, weights = weights)

white_df_q <- data.frame(white_pred_q = predict(quad_white, dwi), bac1 = dwi$bac1)
white_CI_q <- data.frame(white_CI_q = predict(quad_white, dwi, interval = "confidence"),  bac1 = dwi$bac1)

white_graph_q <- ggplot(dwi_cov, aes(x = bac1))+
  geom_point(aes(y= white_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = white_df_q, aes(y =white_pred_q))+
  geom_line(data= white_CI_q, aes(y = white_CI_q.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= white_CI_q, aes(y = white_CI_q.lwr), color = "steelblue", linetype = 3 )+
  xlim(0, 0.2)+
  ylim(0.75,0.9)
white_graph_q


#acc
Y_acc <- dwi$acc

model_acc <- lm(Y_acc ~ Tr + Xl + Xr, weights = weights)
acc_CI <- data.frame(acc_CI = predict(model_acc, dwi, interval = "confidence"),  bac1 = dwi$bac1)

acc_df <- data.frame(acc_pred = predict(model_acc, dwi), bac1 = dwi$bac1)
acc_graph <- ggplot(dwi_cov, aes(x = bac1))+
  geom_point(aes(y= acc_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = acc_df, aes(y = acc_pred))+
  geom_line(data= acc_CI, aes(y = acc_CI.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= acc_CI, aes(y = acc_CI.lwr), color = "steelblue", linetype = 3 )+
  xlim(0, 0.2)+
  ylim(0, 0.3)
acc_graph 

#acc quadratic

quad_acc <- lm(Y_acc ~ Tr + Xl + Xl2 + Xr + Xr2, weights = weights)

acc_df_q <- data.frame(acc_pred_q = predict(quad_acc, dwi), bac1 = dwi$bac1)
acc_CI_q <- data.frame(acc_CI_q = predict(quad_acc, dwi, interval = "confidence"),  bac1 = dwi$bac1)

acc_graph_q <- ggplot(dwi_cov, aes(x = bac1))+
  geom_point(aes(y= acc_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = acc_df_q, aes(y =acc_pred_q))+
  geom_line(data= acc_CI_q, aes(y = acc_CI_q.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= acc_CI_q, aes(y = acc_CI_q.lwr), color = "steelblue", linetype = 3 )+
  xlim(0, 0.2)+
  ylim(0, 0.3)
acc_graph_q

#age
Y_aged <- dwi$aged

model_aged <- lm(Y_aged ~ Tr + Xl + Xr, weights = weights)
aged_CI <- data.frame(aged_CI = predict(model_aged, dwi, interval = "confidence"),  bac1 = dwi$bac1)

aged_df <- data.frame(aged_pred = predict(model_aged, dwi), bac1 = dwi$bac1)
aged_graph <- ggplot(dwi_cov, aes(x = bac1))+
  geom_point(aes(y= aged_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = aged_df, aes(y = aged_pred))+
  geom_line(data= aged_CI, aes(y = aged_CI.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= aged_CI, aes(y = aged_CI.lwr), color = "steelblue", linetype = 3 )+
  xlim(0, 0.2)+
  ylim(32,42)
aged_graph 

#aged quadratic

quad_aged <- lm(Y_aged ~ Tr + Xl + Xl2 + Xr + Xr2, weights = weights)

aged_df_q <- data.frame(aged_pred_q = predict(quad_aged, dwi), bac1 = dwi$bac1)
aged_CI_q <- data.frame(aged_CI_q = predict(quad_aged, dwi, interval = "confidence"),  bac1 = dwi$bac1)

aged_graph_q <- ggplot(dwi_cov, aes(x = bac1))+
  geom_point(aes(y= aged_prop), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = aged_df_q, aes(y =aged_pred_q))+
  geom_line(data= aged_CI_q, aes(y = aged_CI_q.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= aged_CI_q, aes(y = aged_CI_q.lwr), color = "steelblue", linetype = 3 )+
  xlim(0, 0.2)+
  ylim(32,42)
aged_graph_q

#combine graphs
ggarrange(male_graph, male_graph_q, 
          labels = c("Male Linear", "Male Quadratic"),
          ncol = 2, nrow = 1)

ggarrange(white_graph, white_graph_q,
          labels = c("White Linear", "White Quadratic"),
          ncol = 2, nrow = 1)

ggarrange(acc_graph, acc_graph_q, 
          labels = c("Acc Linear", "Acc Quadratic"),
          ncol = 2, nrow = 1)

garrange( aged_graph, aged_graph_q,
          labels = c("Age Linear", "Age Quadratic"),
          ncol = 2, nrow = 1)




#No 7
dwiA= dwi %>% 
  filter(bac1>0.03 & bac1 < 0.13)%>%
  mutate(bac2 = bac1*bac1, bac7 = bac1-0.08, dui = ifelse(bac7 >= 0,1,0))

dwiB= dwi %>% 
  filter(bac1>0.055 & bac1 < 0.105)%>%
  mutate(bac2 = bac1*bac1, bac7 = bac1-0.08, dui = ifelse(bac7 >= 0,1,0))

weightsA <- kernelwts(dwiA$bac7, center = 0, bw = .05, kernel = "rectangular")
weightsB <- kernelwts(dwiB$bac7, center = 0, bw = .05, kernel = "rectangular")

linearA <- lm_robust(recidivism ~ bac7 + dui + aged + white + male + acc,
              data = dwiA,
              weights = weightsA)
linearB <- lm_robust(recidivism ~ bac7 + dui + aged + white + male + acc,
              data = dwiB,
              weights = weightsB)


linearcutA <- lm_robust(recidivism ~ bac7*dui + aged + white + male + acc,
            data = dwiA,
            weights = weightsA)
linearcutB <- lm_robust(recidivism ~ bac7*dui + aged + white + male + acc,
              data = dwiB,
              weights = weightsB)


linquadcutA <- lm_robust(recidivism ~ dui*(bac7+bac2) + aged + white + male + acc,
                 data = dwiA,
                 weights = weightsA)
linquadcutB <- lm_robust(recidivism ~ dui*(bac7+bac2) + aged + white + male + acc,
                 data = dwiB,
                 weights = weightsB)


tab3_models<-data.frame(
  Panels = c("Bandwidth 0.03 to 0.13", "Std. Error", "Bandwidth 0.055 to 0.105", "Std. Error"),
  linear = c(linearA$coefficients[3], linearA$std.error[3], linearB$coefficients[3],linearA$std.error[3]),
  linearint = c(linearcutA$coefficients[3] , linearcutA$std.error[3], linearcutB$coefficients[3], linearcutB$std.error[3]),
  quadint = c(linquadcutA$coefficients[2], linquadcutA$std.error[2], linquadcutB$coefficients[2], linquadcutB$std.error[2])
  
)

tab3_models %>%
  kable(
    col.names = c("Panels", "Linear", "Linear with Interaction", "Quadratic with Interaction"),
    digits = 3,
    caption = "Treatment Effect Models"
  )%>%
  kable_classic(full_width = F, html_font = "Cambria")


#No 8

dwi8 = dwi %>% 
  filter(bac1 <= 0.15)%>%
  group_by(bac1)%>%
  summarize(n=n(), recid_prob = sum(recidivism)/n)


X <- dwi8$bac1 - .08 # center the data
Y <- dwi8$recid_prob
X2 <- X^2 

Xl <- (X < 0) * X # X if X < 0, o.w. 0
Xl2 <- Xl^2

Xr <- (X >= 0) * X # X if X >= 0, o.w. 0
Xr2 <- Xr^2

Tr <- as.integer(X >= 0) # dummy variable


weights <- rdd::kernelwts(X, center = 0, bw = .05, kernel = "rectangular")

model8 <- lm_robust(Y ~ Tr + Xl + Xr, data=dwi8, weights = weights)
quad8 <- lm_robust(Y ~ Tr +Xl2 + Xl +Xr2 + Xr, data=dwi8, weights = weights)

df <- data.frame(pred = predict(model8, dwi8), bac1 = dwi8$bac1)
CI <- data.frame(CI = predict(model8, dwi8, interval = "confidence"),  bac1 = dwi8$bac1)

df_q<- data.frame(pred_q = predict(quad8, dwi8), bac1 = dwi8$bac1)
CI_q <- data.frame(CI_q = predict(quad8, dwi8, interval = "confidence"),  bac1 = dwi8$bac1)

graph <- ggplot(dwi8, aes(x = bac1))+
  geom_point(aes(y= Y), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = df, aes(y = pred))+
  geom_line(data= CI, aes(y = CI.fit.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= CI, aes(y = CI.fit.lwr), color = "steelblue", linetype = 3 )


graph_q <- ggplot(dwi8, aes(x = bac1))+
  geom_point(aes(y= Y), alpha = 0.25)+
  geom_vline(xintercept = 0.08, colour = "blue", linetype = 1)+
  geom_line(data = df_q, aes(y = pred_q))+
  geom_line(data= CI_q, aes(y = CI_q.fit.upr), color = "steelblue", linetype = 3 )+
  geom_line(data= CI_q, aes(y = CI_q.fit.lwr), color = "steelblue", linetype = 3 )


#combine graphs
ggarrange(graph, graph_q,
          labels = c("Linear Model", "Quadratic Model"),
          ncol = 1, nrow = 2)


 

