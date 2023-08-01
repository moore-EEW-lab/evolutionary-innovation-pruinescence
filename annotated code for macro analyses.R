setwd('/Users/michaelmoore/desktop/Working Directory')

library(phytools)
library(car)
library(MASS)
library(nlme)
library(lme4)
library(lmerTest)
library(lsmeans)
library(ggplot2)
library(plyr)
library(phylolm)
library(phylosignal)
library(geiger)

# load data and tree
macro.dat <- read.csv('pruin.macro.dat.csv')
rownames(macro.dat) <- macro.dat$binom


na.drags <- read.tree('pruin.macro.dat.phylo.tre')



# compare the likelihood of pruinescence between reproductive strategies (f=fliers vs p=perchers)
repro.strag01 <- phyloglm(pruin.binom ~ repro.strag, data = macro.dat, boot = 1000, phy = na.drags, method = 'logistic_MPLE', btol = 50)
summary(repro.strag01)

# Call:
# phyloglm(formula = pruin.binom ~ repro.strag, data = macro.dat, 
    # phy = na.drags, method = "logistic_MPLE", btol = 50, boot = 1000)
       # AIC     logLik Pen.logLik 
    # 197.06     -95.53     -93.86 

# Method: logistic_MPLE
# Mean tip height: 0.6627907
# Parameter estimate(s):
# alpha: 53.70996 
      # bootstrap mean: 45.29445 (on log scale, then back transformed)
      # so possible downward bias.
      # bootstrap 95% CI: (1.581641,82.34826)

# Coefficients:
             # Estimate   StdErr  z.value lowerbootCI upperbootCI   p.value    
# (Intercept)  -3.72435  0.66273 -5.61968    -5.38982     -2.5784 1.913e-08 ***
# repro.stragp  2.23413  0.70221  3.18155     1.11591      3.7790  0.001465 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Note: Wald-type p-values for coefficients, conditional on alpha=53.70996
      # Parametric bootstrap results based on 1000 fitted replicates


## compare the likelihood of pruinescence among typical breeding habitats 
macro.dat$typical.habitat <- factor(macro.dat$typical.habitat, levels = c('forest', 'both', 'open'))

hab01 <- phyloglm(pruin.binom ~ typical.habitat, data = macro.dat, boot = 1000, phy = na.drags, method = 'logistic_MPLE', btol = 50)
summary(hab01)

# Call:
# phyloglm(formula = pruin.binom ~ typical.habitat, data = macro.dat, 
    # phy = na.drags, method = "logistic_MPLE", btol = 50, boot = 1000)
       # AIC     logLik Pen.logLik 
    # 208.32    -100.16     -96.89 

# Method: logistic_MPLE
# Mean tip height: 0.6627907
# Parameter estimate(s):
# alpha: 44.0903 
      # bootstrap mean: 38.0509 (on log scale, then back transformed)
      # so possible downward bias.
      # bootstrap 95% CI: (6.072752,82.22503)

# Coefficients:
                    # Estimate   StdErr  z.value lowerbootCI upperbootCI   p.value    
# (Intercept)         -2.72959  0.45603 -5.98552    -3.61215     -1.7059 2.157e-09 ***
# typical.habitatboth  0.83566  0.43883  1.90430     0.02920      1.5584   0.05687 .  
# typical.habitatopen  1.04025  0.41599  2.50067     0.23544      1.8938   0.01240 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Note: Wald-type p-values for coefficients, conditional on alpha=44.0903
      # Parametric bootstrap results based on 1000 fitted replicates
      
### test if species with warmer and drier ranges are more likely to have pruinescence
clim01 <- phyloglm(pruin.binom ~ z.MAT + z.arid, data = macro.dat, boot = 1000, phy = na.drags, method = 'logistic_MPLE', btol = 50)
summary(clim01)

# Call:
# phyloglm(formula = pruin.binom ~ z.MAT + z.arid, data = macro.dat, 
    # phy = na.drags, method = "logistic_MPLE", btol = 50, boot = 1000)
       # AIC     logLik Pen.logLik 
    # 203.34     -97.67     -93.48 

# Method: logistic_MPLE
# Mean tip height: 0.6627907
# Parameter estimate(s):
# alpha: 32.86044 
      # bootstrap mean: 16.57903 (on log scale, then back transformed)
      # so possible downward bias.
      # bootstrap 95% CI: (1.625135,69.41056)

# Coefficients:
            # Estimate   StdErr  z.value lowerbootCI upperbootCI   p.value    
# (Intercept) -2.48360  0.40341 -6.15651    -3.08229     -0.7544 7.437e-10 ***
# z.MAT        0.41641  0.23193  1.79540     0.01715      0.8475   0.07259 .  
# z.arid      -0.37327  0.19364 -1.92767    -0.71532     -0.0181   0.05390 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Note: Wald-type p-values for coefficients, conditional on alpha=32.86044
      # Parametric bootstrap results based on 1000 fitted replicates


## assess the relationship between species climatic conditions and pruinescence on different body segments

# create a couple of binomial response variables that correspond to pruinescence on each body segment
macro.dat$pruin.thor.yn <- ifelse(macro.dat$pruin.thorax == 0, 0, 1)
macro.dat$pruin.abdo.yn <- ifelse(macro.dat$pruin.abdo == 0, 0, 1)

## test thorax first
thorax01 <- phyloglm(pruin.thor.yn ~ z.MAT + z.arid, data = macro.dat, boot = 1000, phy = na.drags, method = 'logistic_MPLE', btol = 50)
summary(thorax01)

# Call:
# phyloglm(formula = pruin.thor.yn ~ z.MAT + z.arid, data = macro.dat, 
    # phy = na.drags, method = "logistic_MPLE", btol = 50, boot = 1000)
       # AIC     logLik Pen.logLik 
    # 175.88     -83.94     -79.67 

# Method: logistic_MPLE
# Mean tip height: 0.6627907
# Parameter estimate(s):
# alpha: 46.38327 
      # bootstrap mean: 0 (on log scale, then back transformed)
      # so possible downward bias.
      # bootstrap 95% CI: (2.643372,82.24536)

# Coefficients:
             # Estimate    StdErr   z.value lowerbootCI upperbootCI   p.value    
# (Intercept) -2.459168  0.361936 -6.794488   -2.954489     -0.8147 1.087e-11 ***
# z.MAT        0.319268  0.226167  1.411647   -0.044325      0.7222   0.15805    
# z.arid      -0.317871  0.193162 -1.645616   -0.664155      0.0309   0.09984 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## now abdomen
abdo01 <- phyloglm(pruin.abdo.yn ~ z.MAT + z.arid, data = macro.dat, boot = 1000, phy = na.drags, method = 'logistic_MPLE', btol = 50)
summary(abdo01)

# Call:
# phyloglm(formula = pruin.abdo.yn ~ z.MAT + z.arid, data = macro.dat, 
    # phy = na.drags, method = "logistic_MPLE", btol = 50, boot = 1000)
       # AIC     logLik Pen.logLik 
    # 205.45     -98.73     -94.58 

# Method: logistic_MPLE
# Mean tip height: 0.6627907
# Parameter estimate(s):
# alpha: 32.529 
      # bootstrap mean: 17.37839 (on log scale, then back transformed)
      # so possible downward bias.
      # bootstrap 95% CI: (1.535803,82.0948)

# Coefficients:
             # Estimate    StdErr   z.value lowerbootCI upperbootCI   p.value    
# (Intercept) -2.524936  0.409854 -6.160570   -3.125213     -0.8534 7.248e-10 ***
# z.MAT        0.464559  0.238581  1.947172    0.055269      0.9365   0.05151 .  
# z.arid      -0.380373  0.196572 -1.935035   -0.696560     -0.0277   0.05299 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1