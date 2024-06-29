library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(knitr)
library(GGally)
library(SciViews)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(binom)
library(DoseFinding)
library(rstanemax)
library(aplore3)





getwd()
setwd("~/Desktop")

#imput data
data <- read.csv("data_110.csv")
#Calculate the effect-variable
data$log_par_delta <- log(data$parasite_0h) - log(data$parasite_6h)
#Graphically examine the difference in the effect and side effect variable with respect to the study arm
data$study_arm <- as.factor(data$study_arm)

plot1 <- ggplot(data=data, aes(x=study_arm,y=log_par_delta)) +  
  geom_boxplot()
plot1

###From plot1, study arms 5 to 7 seem most promising as replacements for artemisinin (which is study arm 1) 

# compute side effects per group
if(!nzchar(system.file(package = "binom"))) install.packages("binom")
if(!nzchar(system.file(package = "dplyr"))) install.packages("dplyr")


bieffekt_data <- data %>% group_by(study_arm) %>%
  summarize(bieffekt=sum(bieffekt),
            binom.confint(x=bieffekt,n=n(), methods="wilson"),
            .groups="drop_last")

plot2 <- ggplot(data=bieffekt_data,aes(x=study_arm,y=mean)) + 
  geom_point(size=5) +
  geom_errorbar(aes(x=study_arm,y=mean,ymin=lower,ymax=upper),width = 0.2) +
  ylab("Proportion (and CI) with side effects") +
  xlab("Study arm")
plot2

data_AMP1050 <- subset(data, drug==2)

plot3 <- ggplot(data=data_AMP1050,aes(x=cmax,y=log_par_delta)) +  
  geom_point() 
plot3 
summary(plot3)
model_liner_reg <- lm( log_par_delta ~ cmax, data = data_AMP1050)
summary(model_liner_reg)
plot1 <- ggplot(data=data_AMP1050,aes(x=cmax,y=log_par_delta)) +  geom_point() 
plot1 + geom_smooth(method="lm")

model_liner_reg_2 <- lm( dose ~ cmax, data = data_AMP1050)
summary(model_liner_reg_2)

#y = 1.412x + 111.7

984  = mean + 256 * 1.96 
mean = 986 - 256 * 1.96 
mean

mean2 = 220 + 256 * 1.96 
mean2




#y = 0.0453x - 3.224
#p value
##as small as it can
#residuals
##as small as it can
#R^2

fit.emax <- stan_emax(log_par_delta ~ cmax, data = data_AMP1050)
fit.emax
extract_stanfit(fit.emax)
plot(fit.emax)
extract_stanfit(fit.emax)
summary(fit.emax)




nls_fit <- nls(log_par_delta ~ E0 + (EMAX*cmax^HILL)/(EC50^HILL + cmax^HILL), data = data_AMP1050, 
                 start = list( E0 = -6.68 , EMAX = 196.99, EC50 = 2918.79, HILL = 2), trace = TRUE)
nls_fit




arm1av <- data %>% filter(study_arm == 1)

#calculations
m = sd(arm1av$log_par_delta)/sqrt(120) * qnorm(0.975)
ci.lower.bound = mean(arm1av$log_par_delta) - m
ci.upper.bound = mean(arm1av$log_par_delta) + m


plot3 + geom_ribbon(aes(ymin = ci.lower.bound, ymax = ci.upper.bound), fill = "red", alpha = 0.5)

glm_cmax_fit <- glm(bieffekt ~ cmax, data = data_AMP1050, 
                   family = binomial(link = "logit"))
summary(glm_cmax_fit)



data_AMP1050 %>% 
  ggplot(aes(x=cmax, y=bieffekt)) + 
  geom_point(shape=1, 
             position=position_jitter(width=.05,height=.05)) + 
  stat_smooth(method="glm", 
              method.args=list(family="binomial"), 
              se=TRUE)

# probability of bieff (amp) = e ^ (-1.0285) = 0.3575 
# largest Cmax that is not inferior to the artemisinin study arm in terms of probability of side effects  
## log(0.291 / 1 - 0.291) = 0.0013 * x + 1.029
###  

##0.291 = 0.0013 * x - 1.029





