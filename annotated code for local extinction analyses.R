library(lme4)
library(MASS)
library(MuMIn)
library(visreg)
library(emmeans)
library(dplyr)
library(sf)
library(lattice)
library(sjPlot)
library(merTools)
library(ggeffects)
library(DHARMa)
library(ggplot2)
library(rlang)
library(MCMCglmm)
library(phytools)

# load data
persist.dat <- read.csv('pruin.persist.dat.csv')

### fit model that estimates how local extinctions (persist = 0) of pruinose and non-pruinose species were related to changing temperature and precipitation 
persist.glmer <- glmer(Persistence ~ Pruinosity + Zprecip + Ztemp +
                                 Pruinosity:Zprecip + Pruinosity:Ztemp +
                                 (1|Species) + (1|FID) + (1|Family) + (1|Genus), data = data, family = 'binomial', 
                               control = glmerControl('bobyqa'), na.action = na.fail)
                          
                          
# look at changes in the likelihood of local extinction in response to changes in PPT for pruinose vs non-pruinose species
emtrends(persist.glmer, specs = 'Pruinosity', var = 'Zprecip')
 # Pruinosity Zprecip.trend    SE  df asymp.LCL asymp.UCL
 # n                  0.318 0.149 Inf    0.0256      0.61
 # y                  0.869 0.290 Inf    0.3001      1.44

# look at changes in the likelihood of local extinction in response to changes in temperature for pruinose vs non-pruinose species
emtrends(persist.glmer, specs = 'Pruinosity', var = 'Ztemp')
 # Pruinosity Ztemp.trend    SE  df asymp.LCL asymp.UCL
 # n               -0.277 0.120 Inf    -0.511   -0.0421
 # y               -0.145 0.224 Inf    -0.584    0.2950
 

## examine if we should be fitting a full phylogenetic model or not. Given uncertainty in topology, etc, only fit if the phylo model fits much better than the naive model. No frequentist implementation for this, so need to use Bayesian MCMCglmm. Compare DIC of models with a full phylogeny vs with a nested random effects structure like that above

# load in tree, but sure the nodes have names
p.drags <- read.tree('pruin.persistence.tre')
p.drags$node.label <- seq(1:59)

# need to make separate column called "animal" so the model knows how to match up the rows with the phylogenetic variance-covariance matrix
persist.dat$animal <- persist.dat$Genus_species

prior.1 <- list(R = list(V=1, fix=1), G = list(G1 = list(V = 1, nu=0.0001), G2 = list(V = 1, nu=0.0001), G3 = list(V = 1, nu=0.0001)))
b.mod01 <- MCMCglmm(as.factor(Persistence) ~ Pruinosity + Ztemp + Zprecip + Pruinosity:Ztemp + Pruinosity:Zprecip, family = 'categorical', random =~Genus_species + FID + animal, prior = prior.1, pedigree = p.drags, node = 'tips', burnin = 50000, nitt = 200000, thin = 200, data = data, pr = TRUE, verbose = FALSE)

prior.2 <- list(R = list(V=1, fix =1), G = list(G1 = list(V = 1, nu = 0.0001), G2 = list(V = 1, nu=0.0001), G3 = list(V = 1, nu=0.0001), G4 = list(V = 1, nu=0.0001)))
b.mod02 <- MCMCglmm(as.factor(Persistence) ~ Pruinosity + Ztemp + Zprecip + Pruinosity:Ztemp + Pruinosity:Zprecip, family = 'categorical', random =~Genus_species + FID + Genus + Family, prior = prior.2, burnin = 50000, nitt = 200000, thin = 200, data = data, pr = TRUE, verbose = FALSE)
DIC(b.mod01, b.mod02)

      # df      DIC
# b.mod01 10 1125.971
# b.mod02 11 1124.515

# models are equal, so better to go with heirarchical nesting structure which makes fewer assumptions