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

### fit models that estimates how local extinctions (extinction: persist = 0) of pruinose and non-pruinose species were related to changing temperature
# first assess if we need to control for background temperature
temp.hist.glmer <- glmer(Persistence ~ Pruinosity + Ztemp + z.hist.temp + 
                                 Pruinosity:Ztemp +
                                 (1|Species) + (1|FID) + (1|Family) + (1|Genus), data = persist.dat, family = 'binomial', 
                               control = glmerControl('bobyqa'), na.action = na.fail)

temp.glmer <- glmer(Persistence ~ Pruinosity + Ztemp + 
                                 Pruinosity:Ztemp +
                                 (1|Species) + (1|FID) + (1|Family) + (1|Genus), data = persist.dat, family = 'binomial', 
                               control = glmerControl('bobyqa'), na.action = na.fail)

AICc(temp.hist.glmer, temp.glmer)
                # df     AICc
# temp.hist.glmer  9 1305.464
# temp.glmer       8 1304.883               

# background temperature doesn't improve the model
                          

# look at changes in the likelihood of local extinction in response to changes in temperature for pruinose vs non-pruinose species
emtrends(temp.glmer, specs = 'Pruinosity', var = 'Ztemp')
 # Pruinosity Ztemp.trend    SE  df asymp.LCL asymp.UCL
 # n               -0.338 0.122 Inf    -0.577   -0.0996
 # y               -0.365 0.224 Inf    -0.804    0.0746
 
### fit models that estimates how local extinctions (extinction: persist = 0) of pruinose and non-pruinose species were related to changing PPT
# first assess if we need to control for background ppt

ppt.hist.glmer <- glmer(Persistence ~ Pruinosity + Zprecip + z.hist.precip + 
                                 Pruinosity:Zprecip +
                                 (1|Species) + (1|FID) + (1|Family) + (1|Genus), data = persist.dat, family = 'binomial', 
                               		control = glmerControl('bobyqa'), na.action = na.fail)

ppt.glmer <- glmer(Persistence ~ Pruinosity + Zprecip +
                                 Pruinosity:Zprecip +
                                 (1|Species) + (1|FID) + (1|Family) + (1|Genus), data = persist.data, family = 'binomial', 
                               	 control = glmerControl('bobyqa'), na.action = na.fail)

AICc(ppt.hist.glmer, ppt.glmer)
               # # df     AICc
# # ppt.hist.glmer  9 1284.524
# # ppt.glmer       8 1298.738

# controlling for background PPT improves the model substantially

# estimated marginal trends
emtrends(ppt.hist.glmer, specs = 'Pruinosity', var = 'Zprecip')
 # Pruinosity Zprecip.trend    SE  df asymp.LCL asymp.UCL
 # n                  0.123 0.164 Inf    -0.197     0.444
 # y                  0.661 0.293 Inf     0.087     1.234


## examine if we should be fitting a full phylogenetic model or not. Given uncertainty in topology, etc, only fit if the phylo model fits much better than the naive model. No frequentist implementation for this, so need to use Bayesian MCMCglmm. Compare DIC of models with a full phylogeny vs with a nested random effects structure like that above

# load in tree, but sure the nodes have names
p.drags <- read.tree('pruin.persistence.tre')
p.drags$node.label <- seq(1:59)

persist.dat$Genus_species<- revalue(persist.dat$Genus_species, c("Phanogomphus_lividus" = "Gomphus_lividus", "Phanogomphus_exilis" = "Gomphus_exilis", "Cordulia_shurtleffii" = "Cordulia_shurtleffi")) # need to align a few species names to what they are called in the tree

# need to make separate column called "animal" so the model knows how to match up the rows with the phylogenetic variance-covariance matrix
persist.dat$animal <- persist.dat$Genus_species

prior.1 <- list(R = list(V=1, fix=1), G = list(G1 = list(V = 1, nu=0.0001), G2 = list(V = 1, nu=0.0001), G3 = list(V = 1, nu=0.0001)))
b.mod01 <- MCMCglmm(as.factor(Persistence) ~ Pruinosity + Ztemp + Pruinosity:Ztemp + z.hist.temp, family = 'categorical', random =~Genus_species + FID + animal, prior = prior.1, pedigree = p.drags, node = 'tips', burnin = 50000, nitt = 200000, thin = 200, data = persist.dat, pr = TRUE, verbose = FALSE)

prior.2 <- list(R = list(V=1, fix =1), G = list(G1 = list(V = 1, nu = 0.0001), G2 = list(V = 1, nu=0.0001), G3 = list(V = 1, nu=0.0001), G4 = list(V = 1, nu=0.0001)))
b.mod02 <- MCMCglmm(as.factor(Persistence) ~ Pruinosity + Ztemp + Pruinosity:Ztemp + z.hist.temp, family = 'categorical', random =~Genus_species + FID + Genus + Family, prior = prior.2, burnin = 50000, nitt = 200000, thin = 200, data = persist.dat, pr = TRUE, verbose = FALSE)

DIC(b.mod01, b.mod02)
        # df      DIC
# b.mod01  9 1123.463
# b.mod02 10 1124.100

# models are equal, so better to go with heirarchical nesting structure which makes fewer biological assumptions about correlation structure
