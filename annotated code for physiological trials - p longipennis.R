setwd('/Users/michaelmoore/desktop/Working Directory')

library(lme4)
library(lmerTest)
library(ggplot2)
library(phytools)
library(plyr)
library(metafor)
library(lsmeans)
library(visreg)
library(car)
library(nlme)
library(MuMIn)
library(ggplot2)
library(phylosignal)
library(phylobase)
library(viridis)

long.heat <- read.csv('long.heat.csv')
long.dry <- read.csv('long.dry.manip.csv')

#### non-linear model to test if removing male pruinescence ('m' for manipulated) causes different patterns of heating compared to males with intact pruinescence ('c' for controls) in P. longipennis
heat01 = nlme(temp.gain ~ SSasymp(time, Asym, R0, lrc),
	fixed = list(Asym ~ as.factor(Treatment) + z.mass, R0 ~ as.factor(Treatment) + z.mass, lrc ~ as.factor(Treatment) + z.mass), 
	random = Asym + R0 + lrc ~ 1|ID,
	data = long.heat, 
	start = list(fixed = c(Asym = c(8,4,1), R0 = c(0, 0.4, 0.05), lrc = c(-3, -3, -0.1))),
	control = list(returnObject = TRUE, 
		tolerance = 0.01, pnlsmaxIter = 1000000, minScale = 1/2000000000, msMaxIter = 10000000, msMaxEval = 10000000),
	na.action = na.exclude)

summary(heat01)

# Nonlinear mixed-effects model fit by maximum likelihood
  # Model: temp.gain ~ SSasymp(time, Asym, R0, lrc) 
  # Data: long.heat
       # AIC      BIC   logLik
  # 255.6681 317.4406 -111.834

# Random effects:
 # Formula: list(Asym ~ 1, R0 ~ 1, lrc ~ 1)
 # Level: ID
 # Structure: General positive-definite, Log-Cholesky parametrization
                 # StdDev     Corr         
# Asym.(Intercept) 2.66555834 As.(I) R0.(I)
# R0.(Intercept)   0.01477997 -0.165       
# lrc.(Intercept)  0.16265562 -0.165  1.000
# Residual         0.21539739              

# Fixed effects:  list(Asym ~ as.factor(Treatment) + z.mass, R0 ~ as.factor(Treatment) +      z.mass, lrc ~ as.factor(Treatment) + z.mass) 
                               # Value Std.Error  DF    t-value p-value
# Asym.(Intercept)           11.376847 0.6819232 311   16.68347  0.0000
# Asym.as.factor(Treatment)m  2.137314 0.9701619 311    2.20305  0.0283
# Asym.z.mass                 1.016169 0.4862783 311    2.08969  0.0375
# R0.(Intercept)             -0.027129 0.0496309 311   -0.54662  0.5850
# R0.as.factor(Treatment)m    0.202660 0.0708342 311    2.86105  0.0045
# R0.z.mass                   0.002260 0.0360069 311    0.06276  0.9500
# lrc.(Intercept)            -4.520413 0.0444213 311 -101.76231  0.0000
# lrc.as.factor(Treatment)m   0.076552 0.0623467 311    1.22785  0.2204
# lrc.z.mass                 -0.162798 0.0318048 311   -5.11864  0.0000
intervals(heat01, which = 'fixed')

# Approximate 95% confidence intervals

 # Fixed effects:
                                 # lower         est.       upper
# Asym.(Intercept)           10.05239409 11.376846781 12.70129948
# Asym.as.factor(Treatment)m  0.25303469  2.137313514  4.02159234
# Asym.z.mass                 0.07170377  1.016168661  1.96063355
# R0.(Intercept)             -0.12352378 -0.027129087  0.06926561
# R0.as.factor(Treatment)m    0.06508370  0.202660116  0.34023653
# R0.z.mass                  -0.06767397  0.002259707  0.07219339
# lrc.(Intercept)            -4.60668992 -4.520413489 -4.43413706
# lrc.as.factor(Treatment)m  -0.04453952  0.076552284  0.19764409
# lrc.z.mass                 -0.22456992 -0.162797574 -0.10102523

### get the estimated marginals
emmeans(heat01, param = 'Asym', specs = 'Treatment')
 # Treatment emmean    SE  df lower.CL upper.CL
 # c           11.4 0.673 311     10.1     12.7
 # m           13.5 0.672 311     12.2     14.8
 
emmeans(heat01, param = 'lrc', specs = 'Treatment')
 # Treatment emmean     SE  df lower.CL upper.CL
 # c          -4.52 0.0438 311    -4.61    -4.43
 # m          -4.44 0.0427 311    -4.53    -4.36
 
 
### linear mixed-effects model to test if removing male pruinesecence results in greater water loss in P longipennis
dry01 <- lmer(prop.loss ~ pre.mass + treatment + (1|date), data = long.dry)
summary(dry01)

# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: prop.loss ~ pre.mass + treatment + (1 | date)
   # Data: long.dry

# REML criterion at convergence: -280.7

# Scaled residuals: 
     # Min       1Q   Median       3Q      Max 
# -2.27299 -0.54046 -0.07522  0.63907  2.99140 

# Random effects:
 # Groups   Name        Variance  Std.Dev.
 # date     (Intercept) 0.0002895 0.01701 
 # Residual             0.0003883 0.01970 
# Number of obs: 62, groups:  date, 7

# Fixed effects:
            # Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)  0.06128    0.01732 55.04806   3.538 0.000826 ***
# pre.mass     0.07138    0.08011 54.98085   0.891 0.376770    
# treatmentm   0.02910    0.00508 54.04199   5.730  4.6e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Correlation of Fixed Effects:
           # (Intr) pr.mss
# pre.mass   -0.904       
# treatmentm -0.105 -0.024

emmeans(dry01, specs = 'treatment')
 # treatment emmean      SE   df lower.CL upper.CL
 # c         0.0751 0.00744 7.00   0.0575   0.0927
 # m         0.1042 0.00769 7.88   0.0865   0.1220
 
confint(dry01)
                   # 2.5 %     97.5 %
# .sig01       0.008746061 0.03121896
# .sigma       0.016232640 0.02355785
# (Intercept)  0.027686739 0.09474779
# pre.mass    -0.084340866 0.22904666
# treatmentm   0.019203757 0.03908245
