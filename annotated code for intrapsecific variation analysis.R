

library(car)
library(MASS)
library(ggplot2)
library(plyr)
library(phytools)
library(metafor)
library(lmerTest)
library(emmeans)
library(raster)
library(sp)
library(MuMIn)

palo.clim <- read.csv('palo.clim.csv')

## test if Pachydiplax longipennis males are more likely to have pruinescence on both their abdomens and thoraxes in either warmer or drier parts of the species range
pru01 <- glmer(pru ~ MAT + arid + (1|State/County), data = palo.clim, na.action = na.exclude, family = 'binomial', control = glmerControl(optimizer = 'bobyqa')) 

summary(pru01)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
 # Family: binomial  ( logit )
# Formula: pru ~ MAT + arid + (1 | State/County)
   # Data: palo.clim
# Control: glmerControl(optimizer = "bobyqa")

     # AIC      BIC   logLik deviance df.resid 
   # 231.1    251.5   -110.5    221.1      431 

# Scaled residuals: 
    # Min      1Q  Median      3Q     Max 
# -4.5354 -0.1848 -0.0260  0.2660  3.2261 

# Random effects:
 # Groups       Name        Variance Std.Dev.
 # County:State (Intercept) 0.4672   0.6835  
 # State        (Intercept) 4.2620   2.0645  
# Number of obs: 436, groups:  County:State, 245; State, 44

# Fixed effects:
            # Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  3.85090    1.33432   2.886   0.0039 ** 
# MAT         -0.08799    0.07398  -1.189   0.2343    
# arid        -0.04415    0.00717  -6.157 7.41e-10 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Correlation of Fixed Effects:
     # (Intr) MAT   
# MAT  -0.803       
# arid -0.616  0.258


confint(pru01, parm = c('MAT','arid'))
           # 2.5 %      97.5 %
# MAT  -0.24262181  0.05954908
# arid -0.06021998 -0.03151886