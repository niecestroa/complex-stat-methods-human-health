library(readr)
library(tidyverse)
library(survival)
library(asaur)

# Loading Data

NKI_cleaned <- read.csv("~/UTH Stochastic Processes in Biostatistics I/PHD1950 Project/PHD1950 Data/NKI_cleaned.csv")

attach(NKI_cleaned)

null.model.NKI <- coxph(Surv(survival, eventdeath) ~ 1 , data=NKI_cleaned)
summary(null.model.NKI)

### Step 2: Find the residuals for the Null Model

rr <- residuals(null.model.NKI)

# Numeric Variables

# need to look at bili copper albumin protime age

# Variable = age diam	posnodes

plot(rr ~ age, xlab="Age", ylab="Residual")
lines(lowess(rr ~ age), col="blue")
# maybe transformation needed - data in beginning is too small so residual not around 0
# try log transformation

plot(rr ~ diam, xlab="Diam", ylab="Residual")
lines(lowess(rr ~ diam), col="blue")
# maybe a transformation is needed

plot(rr ~ posnodes, xlab="Posnodes", ylab="Residual")
lines(lowess(rr ~ posnodes), col="blue")
# looks like transformation would work

plot(rr ~ timerecurrence, xlab="Time Recurrence", ylab="Residual")
lines(lowess(rr ~ timerecurrence), col="blue")
# try log transformation

# variable names - age timerecurrence	chemo	hormonal	amputation	histtype	diam	posnodes	grade	angioinv	lymphinfil
# Numeric - age timerecurrence diam	
# factors - chemo	hormonal	amputation histtype posnodes	grade	angioinv	lymphinfil

firstmodel <- coxph(Surv(survival, eventdeath) ~ age + log(timerecurrence) + diam + 
                      as.factor(chemo) + as.factor(hormonal) + 
                      as.factor(amputation) + as.factor(histtype) + 
                      as.factor(posnodes) + as.factor(grade) + 
                      as.factor(angioinv) + as.factor(lymphinfil))
summary(firstmodel)

# Variable Selection for first model with all variables
step(firstmodel, direction = "both", criterion = "AIC")
step(firstmodel, direction = "both", criterion = "BIC")

firstmodel_trans <- coxph(Surv(survival, eventdeath) ~  log(age) + log(timerecurrence) + diam + 
                      as.factor(chemo) + as.factor(hormonal) + 
                      as.factor(amputation) + as.factor(histtype) + 
                      as.factor(posnodes) + as.factor(grade) + 
                      as.factor(angioinv) + as.factor(lymphinfil))
summary(firstmodel_trans)

# Variable Selection for first model with all variables
step(firstmodel_trans, direction = "both", criterion = "AIC")
step(firstmodel_trans, direction = "both", criterion = "BIC")

# only timerecurrence and grade left in both AIC and BIC selection

basicmodel <- coxph(Surv(survival, eventdeath) ~ as.factor(grade) + timerecurrence )
summary(basicmodel)
# C-Statistic = 0.952
# PH assumption not held for model and timerecurrence

basicmodel_2 <- coxph(Surv(survival, eventdeath) ~ as.factor(grade) + log(timerecurrence) )
summary(basicmodel_2)
# C-Statistic = 0.949
# PH assumption held for model 

transmodel <- coxph(Surv(survival, eventdeath) ~ log(timerecurrence) + 
        diam + as.factor(chemo) + as.factor(grade) + as.factor(lymphinfil))
summary(transmodel)
# C-Statistic = 0.936
# PH assumption held for model 

secondmodel <- coxph(Surv(survival, eventdeath) ~ age + timerecurrence + diam + 
                       as.factor(chemo)*as.factor(hormonal)*as.factor(amputation) + 
                       as.factor(histtype) + as.factor(posnodes) + as.factor(grade) + 
                       as.factor(angioinv) + as.factor(lymphinfil) )
summary(secondmodel)
# C-Statistic = 0.933

# Variable Selection for second model with interaction terms
step(secondmodel, direction = "both", criterion = "AIC")
step(secondmodel, direction = "both", criterion = "BIC")

# variables left are timerecurrence chemo hormonal grade and chemo*grade

interactmodel <- coxph(Surv(survival, eventdeath) ~ log(timerecurrence) + 
                         as.factor(chemo)*as.factor(hormonal) +  
                         as.factor(grade) )
summary(interactmodel)
# C-Statistic = 0.945 and chemo and chemo*hormonal p-values > 0.05
# PH assumption not held for model and timerecurrence

secondmodel_trans <- coxph(Surv(survival, eventdeath) ~ log(age) + log(timerecurrence) + diam + 
                       as.factor(chemo)*as.factor(hormonal)*as.factor(amputation) + 
                       as.factor(histtype) + as.factor(posnodes) + as.factor(grade) + 
                       as.factor(angioinv) + as.factor(lymphinfil) )
summary(secondmodel_trans)
# C-Statistic = 0.918

# Variable Selection for second model with interaction terms
step(secondmodel_trans, direction = "both", criterion = "AIC")
step(secondmodel_trans, direction = "both", criterion = "BIC")
# too complicated to run alone -> too many variables

## Check proportional hazards assumption

(cox.prop.assump1 <- cox.zph(basicmodel) )
(cox.prop.assump1.1 <- cox.zph(basicmodel_2) )
(cox.prop.assump2 <- cox.zph(transmodel) )
(cox.prop.assump3 <- cox.zph(interactmodel) )

## Final Model Choices from Aaron Niecestro

finalmodel_1_AN <- coxph(Surv(survival, eventdeath) ~ log(timerecurrence) + as.factor(grade) )
summary(finalmodel_1_AN)
# C-Statistic = 0.949
# PH assumption held for model 
# All variables are significant and PH assumptions hold for each as well

finalmodel_2_AN <- coxph(Surv(survival, eventdeath) ~ log(timerecurrence) + as.factor(chemo) + as.factor(grade) + as.factor(lymphinfil))
summary(finalmodel_2_AN)
# C-Statistic = 0.935
# PH assumption held for model 
# All variables are significant and PH assumptions hold for each as well

# plots final model 

library(ggplot2)
library(survminer)

par(mfrow=c(1,1))
plot(survfit(finalmodel_1_AN), 
     ylab = "Probability of Survival", 
     xlab = "Time", col = c("black", "black", "black"),
     main = "Cox Proportional Hazards Model",
     conf.int = TRUE)

plot(cox.zph(finalmodel_1_AN)) # add to Google doc
ggcoxzph(cox.zph(finalmodel_1_AN))

# DFbetas for Cox Model 

# for Variable 1 - albunim

rr.dfbeta <- residuals(finalmodel_1_AN, type="dfbeta")
subject.index <- 1:nrow(rr.dfbeta)
plot(rr.dfbeta[,1] ~ subject.index, axes=F)
axis(1, at=seq(from=1, to=length(subject.index), by=5))
axis(2)
title("DFbeta for Final Model for Albunim")

# Identify outliers in variable log(timerecurrence)

# next 2 lines take too long to load

#identify(rr.dfbeta[,1] ~ subject.index)
#identify(resid.dfbeta[,1] ~ index.obs)

#which(resid.dfbeta[,1] > 0.05) # outliers identified here
#which(resid.dfbeta[,1] < -0.05)