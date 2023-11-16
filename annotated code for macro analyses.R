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

## prep data for pagel's test
pruin.yn <- ifelse(macro.dat$pruin.binom == 0, 'non', 'pruin')
names(pruin.yn) <- macro.dat$binom

repro.strag <- macro.dat$repro.strag
names(repro.strag) <- macro.dat$binom

# fit four possible models of correlated evolution
pagel.repro.xy <- fitPagel(tree = na.drags, x = pruin.yn, y = repro.strag)
pagel.repro.x <- fitPagel(tree = na.drags, x = pruin.yn, y = repro.strag, dep.var = 'x')
pagel.repro.y <- fitPagel(tree = na.drags, x = pruin.yn, y = repro.strag, dep.var = 'y')
repro.aic <- c(pagel.repro.xy$independent.AIC, pagel.repro.x$dependent.AIC, pagel.repro.y$dependent.AIC, pagel.repro.xy$dependent.AIC)

repro.aic
#397.0427 381.8378 398.4935 385.6343
# best model is that pruin depends on repro strag. 

trans.per.mill <- ((pagel.repro.x$dependent.Q)* 0.335548173)/200 # scale transition rates to known time dating for groups in analysis (Kohli et al. 2021, iScience)
pagel.repro.x$dependent.Q <- trans.per.mill

pagel.repro.x$dependent.Q # transition rates-
               # non|f        non|p       pruin|f      pruin|p
# non|f   -0.007956896  0.007814845  0.0001420513  0.000000000
# non|p    0.007509774 -0.041663280  0.0000000000  0.034153506
# pruin|f  0.111199104  0.000000000 -0.1190139484  0.007814845
# pruin|p  0.000000000  0.111199104  0.0075097738 -0.118708878


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


## prep data for Pagel's test
hab.open <- ifelse(macro.dat$typical.habitat == 'forest', 'closed', 'open')
names(hab.open) <- macro.dat$binom

# run 4 possible models of correlated evolution
pagel.hab.xy <- fitPagel(tree = na.drags, x = pruin.yn, y = hab.open)
pagel.hab.x <- fitPagel(tree = na.drags, x = pruin.yn, y = hab.open, dep.var = 'x')
pagel.hab.y <- fitPagel(tree = na.drags, x = pruin.yn, y = hab.open, dep.var = 'y')
hab.aic <- c(pagel.hab.xy$independent.AIC, pagel.hab.x$dependent.AIC, pagel.hab.y$dependent.AIC, pagel.hab.xy$dependent.AIC)

hab.aic
#656.2248 652.9862 652.6614 653.5145

hab.trans.per.mill <- ((pagel.hab.xy$dependent.Q)* 0.335548173)/200 # scale transition rates to known time dating for groups in analysis (Kohli et al. 2021, iScience)

pagel.hab.xy$dependent.Q <- hab.trans.per.mill
pagel.hab.xy$dependent.Q # transition rates

                   # non|closed   non|open pruin|closed  pruin|open
# non|closed   -0.11119910  0.1111991   0.00000000  0.00000000
# non|open      0.09511276 -0.1217264   0.00000000  0.02661363
# pruin|closed  0.11119910  0.0000000  -0.22239821  0.11119910
# pruin|open    0.00000000  0.1111991   0.03915986 -0.15035897
  
      
###test if reproductive strategy evolves in a correlated fashion with microhabitat       

## run 4 possible models      
pagel.hr.xy <- fitPagel(tree = na.drags, x = repro.strag, y = hab.open)
pagel.hr.x <- fitPagel(tree = na.drags, x = repro.strag, y = hab.open, dep.var = 'x')
pagel.hr.y <- fitPagel(tree = na.drags, x = repro.strag, y = hab.open, dep.var = 'y')
hr.aic <- c(pagel.hr.xy$independent.AIC, pagel.hr.x$dependent.AIC, pagel.hr.y$dependent.AIC, pagel.hr.xy$dependent.AIC)
hr.aic
#622.1985 622.7133 622.1504 622.9484
# models are all about equal, so go with the simplest. Which says that there's no correlated evolution between microhabitat and reproductive strategy
      
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


## test if mating strategy is correlated with macroclimate
macro.dat$percher.yn <- ifelse(macro.dat$repro.strag == 'p', 1, 0)

perch01 <- phyloglm(percher.yn ~ z.MAT + z.arid, data = macro.dat, boot = 1000, phy = na.drags, method = 'logistic_MPLE', btol = 50)
summary(perch01)

# Call:
# phyloglm(formula = percher.yn ~ z.MAT + z.arid, data = macro.dat, 
    # phy = na.drags, method = "logistic_MPLE", btol = 50, boot = 1000)
       # AIC     logLik Pen.logLik 
    # 179.66     -85.83     -79.35 

# Method: logistic_MPLE
# Mean tip height: 0.6627907
# Parameter estimate(s):
# alpha: 6.584329 
      # bootstrap mean: 0 (on log scale, then back transformed)
      # so possible downward bias.
      # bootstrap 95% CI: (1.591721,12.4183)

# Coefficients:
             # Estimate    StdErr   z.value lowerbootCI upperbootCI p.value
# (Intercept) -0.011651  0.506916 -0.022984   -0.689793      0.7515  0.9817
# z.MAT       -0.010895  0.058405 -0.186547   -0.079052      0.0572  0.8520
# z.arid       0.056697  0.053083  1.068083   -0.023959      0.1270  0.2855


## rate of gains and losses of pruinescence

pruin.mod <- fitDiscrete(phy = na.drags, dat = pruin.yn, model = 'ARD', control = list(niter = 100, hessian = TRUE, CI = 95))
pruin.mod$opt

                # non     pruin
    # non   -10.74312  10.74312
    # pruin  79.93264 -79.93264


(10.74312 * 0.335548173)/200 # gain of pruin = 0.01802417 transitions per million years 
(79.93264 * 0.335548173)/200 # loss of pruin = 0.1341063 transitions per million years 


1/(0.1341063/1000000) # 1 loss every 7,456,771 years
1/(0.01802417/1000000) # 1 gain every 55,481,057 years
