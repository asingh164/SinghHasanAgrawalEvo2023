#Packages required
require(car)
require(MCMCglmm)
require(QGglmm)
require(coda)

QGdata = read.csv(file = "Dataset_June_26.csv", header = T)
QGdata$block = as.factor(QGdata$block)
QGdata$agrawal.stock.number = as.factor(QGdata$agrawal.stock.number)
QGdata$DSPR.stock.number = as.factor(QGdata$DSPR.stock.number)
QGdata$prop.wt = QGdata$total.wt / QGdata$total.flies
QGdata$total.wtNOT = QGdata$total.flies - QGdata$total.wt


#This prior specification was the same for all models

prior2B = list(R = list(V = diag(2), nu = 1),
              G = list(G1 = list(V=diag(2), nu=1, alpha.mu = c(0, 0), alpha.V = 1000*diag(2)), G2 = list(V=diag(2), nu=1, alpha.mu = c(0, 0), alpha.V = 1000*diag(2))))


## For QG ANALYSIS OF VIAL AND CAGE DATA SEPERATELY  
vial.data = subset(QGdata, QGdata$mating.regime == "vial")
cage.data = subset(QGdata, QGdata$mating.regime == "cage")

## For QG ANALYSIS OF FEMALE AND MALE DATA SEPERATELY
female.bothEnv.data = subset(QGdata, QGdata$sex == "f")
male.bothEnv.data = subset(QGdata, QGdata$sex == "m")

#MCMCglmm code for analysis of vial data

MCMCmodel.vial.prop.B = MCMCglmm(cbind(total.wt, total.wtNOT) ~ sex - 1,
                               random = ~us(sex):DSPR.stock.number + idh(sex):block,
                               rcov = ~idh(sex):units,
                               data = vial.data,
                               prior = prior2B,
                               family = "multinomial2",
                               verbose = TRUE,
                               nitt = 5*10^5, burnin = 10000, thin = 100
                               
)



MCMCmodel.cage.prop.B = MCMCglmm(cbind(total.wt, total.wtNOT) ~ sex - 1,
                                 random = ~us(sex):DSPR.stock.number + idh(sex):block,
                                 rcov = ~idh(sex):units,
                                 data = cage.data,
                                 prior = prior2B,
                                 family = "multinomial2",
                                 verbose = TRUE,
                                 nitt = 5*10^5, burnin = 10000, thin = 100
                                 
)

MCMCmodel.female.bothEnv.prop.B = MCMCglmm(cbind(total.wt, total.wtNOT) ~ mating.regime - 1,
                                 random = ~us(mating.regime):DSPR.stock.number + idh(mating.regime):block,
                                 rcov = ~idh(mating.regime):units,
                                 data = female.bothEnv.data,
                                 prior = prior2B,
                                 family = "multinomial2",
                                 verbose = TRUE,
                                 nitt = 5*10^5, burnin = 10000, thin = 100
                                 
)

MCMCmodel.male.bothEnv.prop.B = MCMCglmm(cbind(total.wt, total.wtNOT) ~ mating.regime - 1,
                                         random = ~us(mating.regime):DSPR.stock.number + idh(mating.regime):block,
                                         rcov = ~idh(mating.regime):units,
                                         data = male.bothEnv.data,
                                         prior = prior2B,
                                         family = "multinomial2",
                                         verbose = TRUE,
                                         nitt = 5*10^5, burnin = 10000, thin = 100
                                         
)

#save(MCMCmodel.vial.prop.B, MCMCmodel.cage.prop.B, MCMCmodel.female.bothEnv.prop.B, MCMCmodel.male.bothEnv.prop.B, file = "/Users/aneil/Google\ Drive/AgrawalLab/Amardeep/GeneticArchitecture/QGanalysis/QGmodelsForRMF.Rdata")


# autocorr(MCMCmodel.vial.prop.B$VCV[,"sexf:sexf.DSPR.stock.number"])
# autocorr(MCMCmodel.vial.prop.B$VCV[,"sexm:sexm.DSPR.stock.number"])
# autocorr(MCMCmodel.cage.prop.B$VCV[,"sexf:sexf.DSPR.stock.number"])
# autocorr(MCMCmodel.cage.prop.B$VCV[,"sexm:sexm.DSPR.stock.number"])

#As a check, visually inspect the distribution of variance estimates to make sure that the distributions are not bunched up at 0, if not we can probably trust the HPDIntervals as as a test of significance
hist(mcmc(MCMCmodel.vial.prop.B$VCV)[,"sexf:sexf.DSPR.stock.number"])
hist(mcmc(MCMCmodel.vial.prop.B$VCV)[,"sexm:sexm.DSPR.stock.number"])
hist(mcmc(MCMCmodel.cage.prop.B$VCV)[,"sexf:sexf.DSPR.stock.number"])
hist(mcmc(MCMCmodel.cage.prop.B$VCV)[,"sexm:sexm.DSPR.stock.number"])


## values for use by QGglmm (see de Villemeruil, "How to use the QGglmm package")
b.mu.vial.prop = colMeans(MCMCmodel.vial.prop.B$Sol)
b.G.vial.prop = matrix(colMeans(MCMCmodel.vial.prop.B$VCV)[1:4], ncol = 2)
b.R.vial.prop = diag(colMeans(MCMCmodel.vial.prop.B$VCV)[5:6]) + diag(colMeans(MCMCmodel.vial.prop.B$VCV)[7:8])
b.P.vial.prop = b.G.vial.prop + b.R.vial.prop
## 
b.mu.cage.prop = colMeans(MCMCmodel.cage.prop.B$Sol)
b.G.cage.prop = matrix(colMeans(MCMCmodel.cage.prop.B$VCV)[1:4], ncol = 2)
b.R.cage.prop = diag(colMeans(MCMCmodel.cage.prop.B$VCV)[5:6]) + diag(colMeans(MCMCmodel.cage.prop.B$VCV)[7:8])
b.P.cage.prop = b.G.cage.prop + b.R.cage.prop
##
b.mu.female.bothEnv.prop = colMeans(MCMCmodel.female.bothEnv.prop.B$Sol)
b.G.female.bothEnv.prop = matrix(colMeans(MCMCmodel.female.bothEnv.prop.B$VCV)[1:4], ncol = 2)
b.R.female.bothEnv.prop = diag(colMeans(MCMCmodel.female.bothEnv.prop.B$VCV)[5:6]) + diag(colMeans(MCMCmodel.female.bothEnv.prop.B$VCV)[7:8])
b.P.female.bothEnv.prop = b.G.female.bothEnv.prop + b.R.female.bothEnv.prop
##
b.mu.male.bothEnv.prop = colMeans(MCMCmodel.male.bothEnv.prop.B$Sol)
b.G.male.bothEnv.prop = matrix(colMeans(MCMCmodel.male.bothEnv.prop.B$VCV)[1:4], ncol = 2)
b.R.male.bothEnv.prop = diag(colMeans(MCMCmodel.male.bothEnv.prop.B$VCV)[5:6]) + diag(colMeans(MCMCmodel.male.bothEnv.prop.B$VCV)[7:8])
b.P.male.bothEnv.prop = b.G.male.bothEnv.prop + b.R.male.bothEnv.prop


## average number of offspring per replicate
n.avg.counts.vial.f = mean(vial.data$total.flies[vial.data$sex =="f"], na.rm = T)
n.avg.counts.vial.m = mean(vial.data$total.flies[vial.data$sex =="m"], na.rm = T)

n.avg.counts.cage.f = mean(cage.data$total.flies[cage.data$sex =="f"], na.rm = T)
n.avg.counts.cage.m = mean(cage.data$total.flies[cage.data$sex =="m"], na.rm = T)

## female fitness - vial
b.vial.female.fit.ObsScale.GroupOf4 = QGparams(mu = b.mu.vial.prop[1], var.a = b.G.vial.prop[1,1], var.p = b.P.vial.prop[1,1], model = "binomN.logit", n.obs = n.avg.counts.vial.f)

## male fitness - vial
b.vial.male.fit.ObsScale.GroupOf4 = QGparams(mu = b.mu.vial.prop[2], var.a = b.G.vial.prop[2,2], var.p = b.P.vial.prop[2,2], model = "binomN.logit", n.obs = n.avg.counts.vial.m)

## female fitness - cage
b.cage.female.fit.ObsScale.GroupOf4 = QGparams(mu = b.mu.cage.prop[1], var.a = b.G.cage.prop[1,1], var.p = b.P.cage.prop[1,1], model = "binomN.logit", n.obs = n.avg.counts.cage.f)

## male fitness - cage
b.cage.male.fit.ObsScale.GroupOf4 = QGparams(mu = b.mu.cage.prop[2], var.a = b.G.cage.prop[2,2], var.p = b.P.cage.prop[2,2], model = "binomN.logit", n.obs = n.avg.counts.cage.m)

## Get G-matrix wrt to observed scale - vial
b.vial.bothSexes.fit.ObsScale.GroupOf4 = QGmvparams(mu = b.mu.vial.prop, vcv.G = b.G.vial.prop, vcv.P = b.P.vial.prop, model = c("binomN.logit", "binomN.logit"), n.obs = c(n.avg.counts.vial.f, n.avg.counts.vial.m))

## Get G-matrix wrt to observed scale - females (both environments)
b.female.bothEnv.fit.ObsScale.GroupOf4 = QGmvparams(mu = b.mu.vial.prop, vcv.G = b.G.vial.prop, vcv.P = b.P.vial.prop, model = c("binomN.logit", "binomN.logit"), n.obs = c(n.avg.counts.vial.f, n.avg.counts.vial.m))




##function below returns "per individual" estimates of 
## (1) mean number of offspring, wbar
## (2) phenotypic variance (vP)
## (3) additive genetic variance (vA)
## (4) h^2 = vA/vP
## (5) evolvability vA / (wbar^2)
makeAdjustmentsToSingleSexResultsForIndividualLevel<-function(QGglmm.output1sex){
  z<-as.numeric(c(QGglmm.output1sex[1]/4, QGglmm.output1sex[2]/4 - (3/16)*QGglmm.output1sex[3], QGglmm.output1sex[3]/8, NA, NA))
  z[4]<-z[3]/z[2]
  z[5]<- z[3]/(z[1]^2)
  return(z)
}

makeAdjustmentsToBothSexesResultsForIndividualLevel<-function(QGglmm.output2sexex){
  z<-as.numeric(c(QGglmm.output2sexex$vcv.G.obs[1,1]/8, QGglmm.output2sexex$vcv.G.obs[1,2]/8, QGglmm.output2sexex$vcv.G.obs[2,2]/8, NA))
  z[4]<-z[2]/sqrt(z[1]*z[3])
  return(z)
}


### The function below takes the output MCMCglmm and calculates
###  the posterior distributions for 
### individual level measures (rather than measures for "groups of 4 individuals")
### (1) mean offspring per focal female (assumes total number of offspring scored per vial is equal to average of observed values)
### (2) V(P) for females, i.e., phenotypic variance in offspring per focal female
### (3) V(A) for females
### (4) h^2 for females
### (5) evolvability (= V(A)/wbar^2) for females
### (6-10) as above for males
### (11-13) elements of G-matrix from 2-sex analysis: Gff, Gfm, Gmm
### (14) rmf
### (15) V(P_SA) "phenotypic variance in SA" (SA & SC are measured wrt to relative fitness)
### (16) V(P_SC)
### (17) V(A_SA)
### (18) V(A_SC)
### (19) pseudo h^2 for SA
### (20) pseudo h^2 for SC
getResultsFromQGglmm<-function(MCMCglmm.output, n.counts.f, n.counts.m){
  df<-data.frame(f.mean.obs = NA, f.var.obs = NA, f.var.a.obs = NA, f.h2.obs = NA, f.evolv.obs = NA,
                 m.mean.obs = NA, m.var.obs = NA, m.var.a.obs = NA, m.h2.obs = NA, m.evolv.obs = NA,
                 G11 = NA, G12 = NA, G22 = NA, rmf = NA,
                 SA.var.obs = NA, SC.var.obs = NA, SA.a.var.obs = NA, SC.a.var.obs = NA,
                 SA.h2.obs = NA, SC.h2.obs = NA)
  ## values for QGglmm (de Villemeruil "How to use the QGglmm package")
  for(i in 1:nrow(MCMCglmm.output$VCV)){ #nrow(MCMCglmm.output$VCV)
    z.mu.vial.prop = MCMCglmm.output$Sol[i,]
    z.G.vial.prop = matrix(MCMCglmm.output$VCV[i,1:4], ncol = 2)
    z.R.vial.prop = diag(MCMCglmm.output$VCV[i,5:6]) + diag(MCMCglmm.output$VCV[i,7:8])
    z.P.vial.prop = z.G.vial.prop + z.R.vial.prop

    z.vial.female.fit.ObsScale.GroupOf4 = QGparams(mu = z.mu.vial.prop[1], var.a = z.G.vial.prop[1,1], var.p = z.P.vial.prop[1,1], model = "binomN.logit", n.obs = n.counts.f, verbose = F)
    z.vial.male.fit.ObsScale.GroupOf4 = QGparams(mu = z.mu.vial.prop[2], var.a = z.G.vial.prop[2,2], var.p = z.P.vial.prop[2,2], model = "binomN.logit", n.obs = n.counts.m, verbose = F)
    z.vial.bothSexes.fit.ObsScale.GroupOf4= QGmvparams(mu = z.mu.vial.prop, vcv.G = z.G.vial.prop, vcv.P = z.P.vial.prop, model = c("binomN.logit", "binomN.logit"), n.obs = c(n.counts.f, n.counts.m), verbose = F)
    
    zz = c(makeAdjustmentsToSingleSexResultsForIndividualLevel(z.vial.female.fit.ObsScale.GroupOf4),
           makeAdjustmentsToSingleSexResultsForIndividualLevel(z.vial.male.fit.ObsScale.GroupOf4),
           makeAdjustmentsToBothSexesResultsForIndividualLevel(z.vial.bothSexes.fit.ObsScale.GroupOf4),
           NA, NA, NA, NA, NA, NA)
    
    zz[15] = (1/2)*(zz[2]/(zz[1]^2) + zz[7]/(zz[6]^2) ) - zz[12]/(zz[1]*zz[6])
    zz[16] = (1/2)*(zz[2]/(zz[1]^2) + zz[7]/(zz[6]^2) ) + zz[12]/(zz[1]*zz[6])
    zz[17] = (1/2)*(zz[11]/(zz[1]^2) + zz[13]/(zz[6]^2) ) - zz[12]/(zz[1]*zz[6])
    zz[18] = (1/2)*(zz[11]/(zz[1]^2) + zz[13]/(zz[6]^2) ) + zz[12]/(zz[1]*zz[6])
    zz[19] = zz[17]/zz[15]
    zz[20] = zz[18]/zz[16]
    
    df[i,] = zz
  }
  return(df)
}

b.posterior.obs.scale.vial = getResultsFromQGglmm(MCMCmodel.vial.prop.B, n.avg.counts.vial.f, n.avg.counts.vial.m)

b.posterior.obs.scale.cage = getResultsFromQGglmm(MCMCmodel.cage.prop.B, n.avg.counts.cage.f, n.avg.counts.cage.m)

b.posterior.obs.scale.female.bothEnv = getResultsFromQGglmm(MCMCmodel.female.bothEnv.prop.B, n.avg.counts.cage.f, n.avg.counts.vial.f)
b.posterior.obs.scale.male.bothEnv = getResultsFromQGglmm(MCMCmodel.male.bothEnv.prop.B, n.avg.counts.cage.m, n.avg.counts.vial.m)


summarizeResultsFromPosteriors<-function(posteriors.df){
  df = data.frame( avg = colMeans(posteriors.df), lower = NA, upper = NA)
  credibleIntervals = HPDinterval(as.mcmc(posteriors.df))
  df$lower = as.numeric(credibleIntervals[,1])
  df$upper = as.numeric(credibleIntervals[,2])
  return(df)
}

## Make summary of for analzyes of both sexes, (then display for both environments)
{vial.summary.b <-summarizeResultsFromPosteriors(b.posterior.obs.scale.vial)
colnames(vial.summary.b) = paste("vial", colnames(vial.summary.b), sep = ".")
cage.summary.b <-summarizeResultsFromPosteriors(b.posterior.obs.scale.cage)
colnames(cage.summary.b) = paste("cage", colnames(cage.summary.b), sep = ".")
full.summary.b <-cbind(vial.summary.b, cage.summary.b)}
full.summary.b

## Make summary of for analzyes of both environments, (then display for both sexes)
{female.bothEnv.summary.b <-summarizeResultsFromPosteriors(b.posterior.obs.scale.female.bothEnv)
  colnames(female.bothEnv.summary.b) = paste("female", colnames(female.bothEnv.summary.b), sep = ".")
  male.bothEnv.summary.b <-summarizeResultsFromPosteriors(b.posterior.obs.scale.male.bothEnv)
  colnames(male.bothEnv.summary.b) = paste("male", colnames(male.bothEnv.summary.b), sep = ".")
  full.bothEnv.summary.b <-cbind(female.bothEnv.summary.b, male.bothEnv.summary.b)
  rownames(full.bothEnv.summary.b) = c("cage.mean.obs", "cage.var.obs", "cage.var.a.obs", "cage.h2.obs", "meaningless1",
                                       "vial.mean.obs", "vial.var.obs", "vial.var.a.obs", "vial.h2.obs", "meaningless2",
                                       "G11", "G12", "G22", "r_env",
                                       "EnvAntag.var.obs", "EnvConcord.var.obs", "EnvAntag.a.var.obs", "EnvConcord.a.var.obs",
                                       "EnvAntag.h2.obs", "EnvConcord.h2.obs")}


full.bothEnv.summary.b[1:4,]


 #Save MCMCglmm output as RDS object
#saveRDS(MCMCmodel.vial, file = "MCMCglmm.output.vial.rds")



### Here I permute replicates within blocks (for a given sex & environment).
### Then re-run the QG analyses. RIL variances are MUCH lower (as expected)

get.One.Permuted.vial.data<-function(){
  perm.vial.data = vial.data[0,]
  m.data = vial.data[vial.data$sex == "m",]
  f.data = vial.data[vial.data$sex == "f",]
  for(this.block in levels(vial.data$block)){
    this.block.real.m = m.data[m.data$block == this.block,]
    n.m = nrow(this.block.real.m)
    if(n.m > 0){
      m.x = sample(1:n.m, n.m, replace = F)
      this.block.perm.m = this.block.real.m[m.x, ]
      this.block.perm.m$DSPR.stock.number = this.block.real.m$DSPR.stock.number
      this.block.perm.m$agrawal.stock.number = this.block.real.m$agrawal.stock.number
      perm.vial.data = rbind(perm.vial.data, this.block.perm.m)
    }
    this.block.real.f = f.data[f.data$block == this.block,]
    n.f = nrow(this.block.real.f)
    if(n.f > 0){
      f.x = sample(1:n.f, n.f, replace = F)
      this.block.perm.f = this.block.real.f[f.x, ]
      this.block.perm.f$DSPR.stock.number = this.block.real.f$DSPR.stock.number
      this.block.perm.f$agrawal.stock.number = this.block.real.f$agrawal.stock.number
      perm.vial.data = rbind(perm.vial.data, this.block.perm.f)
    }
  }
  return(perm.vial.data)
}



get.One.Permuted.cage.data<-function(){
  perm.cage.data = cage.data[0,]
  m.data = cage.data[cage.data$sex == "m",]
  f.data = cage.data[cage.data$sex == "f",]
  for(this.block in levels(cage.data$block)){
    this.block.real.m = m.data[m.data$block == this.block,]
    n.m = nrow(this.block.real.m)
    if(n.m > 0){
      m.x = sample(1:n.m, n.m, replace = F)
      this.block.perm.m = this.block.real.m[m.x, ]
      this.block.perm.m$DSPR.stock.number = this.block.real.m$DSPR.stock.number
      this.block.perm.m$agrawal.stock.number = this.block.real.m$agrawal.stock.number
      perm.cage.data = rbind(perm.cage.data, this.block.perm.m)
    }
    this.block.real.f = f.data[f.data$block == this.block,]
    n.f = nrow(this.block.real.f)
    if(n.f > 0){
      f.x = sample(1:n.f, n.f, replace = F)
      this.block.perm.f = this.block.real.f[f.x, ]
      this.block.perm.f$DSPR.stock.number = this.block.real.f$DSPR.stock.number
      this.block.perm.f$agrawal.stock.number = this.block.real.f$agrawal.stock.number
      perm.cage.data = rbind(perm.cage.data, this.block.perm.f)
    }
  }
  return(perm.cage.data)
}


permuted.vial.data.1 = get.One.Permuted.vial.data()
nrow(permuted.vial.data.1)

MCMCmodel.vial.prop.PERM.1 = MCMCglmm(cbind(total.wt, total.wtNOT) ~ sex - 1,
                                      random = ~us(sex):DSPR.stock.number + idh(sex):block,
                                      rcov = ~idh(sex):units,
                                      data = permuted.vial.data.1,
                                      prior = prior2B,
                                      family = "multinomial2",
                                      verbose = TRUE,
                                      nitt = 5*10^5, burnin = 10000, thin = 100
                                      #nitt = 1000000, burnin = 10000, thin = 1000
                                      
)

permuted.cage.data.1 = get.One.Permuted.cage.data()
nrow(permuted.cage.data.1)

MCMCmodel.cage.prop.PERM.1 = MCMCglmm(cbind(total.wt, total.wtNOT) ~ sex - 1,
                                      random = ~us(sex):DSPR.stock.number + idh(sex):block,
                                      rcov = ~idh(sex):units,
                                      data = permuted.cage.data.1,
                                      prior = prior2B,
                                      family = "multinomial2",
                                      verbose = TRUE,
                                      nitt = 5*10^5, burnin = 10000, thin = 100
                                      #nitt = 1000000, burnin = 10000, thin = 1000
                                      
)

b.posterior.obs.scale.vial.PERM1 = getResultsFromQGglmm(MCMCmodel.vial.prop.PERM.1, n.avg.counts.vial.f, n.avg.counts.vial.m)
b.posterior.obs.scale.cage.PERM1 = getResultsFromQGglmm(MCMCmodel.cage.prop.PERM.1, n.avg.counts.vial.f, n.avg.counts.vial.m)

{vial.summary.b.PERM1 <-summarizeResultsFromPosteriors(b.posterior.obs.scale.vial.PERM1)
  colnames(vial.summary.b.PERM1) = paste("vial", colnames(vial.summary.b.PERM1), sep = ".")
  cage.summary.b.PERM1 <-summarizeResultsFromPosteriors(b.posterior.obs.scale.cage.PERM1)
  colnames(cage.summary.b.PERM1) = paste("cage", colnames(cage.summary.b.PERM1), sep = ".")
  full.summary.b.PERM1 <-cbind(vial.summary.b.PERM1, cage.summary.b.PERM1)}


permuted.vial.data.2 = get.One.Permuted.vial.data()
nrow(permuted.vial.data.2)

MCMCmodel.vial.prop.PERM.2 = MCMCglmm(cbind(total.wt, total.wtNOT) ~ sex - 1,
                                      random = ~us(sex):DSPR.stock.number + idh(sex):block,
                                      rcov = ~idh(sex):units,
                                      data = permuted.vial.data.2,
                                      prior = prior2B,
                                      family = "multinomial2",
                                      verbose = TRUE,
                                      nitt = 5*10^5, burnin = 10000, thin = 100
                                      #nitt = 1000000, burnin = 10000, thin = 1000
                                      
)

permuted.cage.data.2 = get.One.Permuted.cage.data()
nrow(permuted.cage.data.2)

MCMCmodel.cage.prop.PERM.2 = MCMCglmm(cbind(total.wt, total.wtNOT) ~ sex - 1,
                                      random = ~us(sex):DSPR.stock.number + idh(sex):block,
                                      rcov = ~idh(sex):units,
                                      data = permuted.cage.data.2,
                                      prior = prior2B,
                                      family = "multinomial2",
                                      verbose = TRUE,
                                      nitt = 5*10^5, burnin = 10000, thin = 100
                                      #nitt = 1000000, burnin = 10000, thin = 1000
                                      
)

save(MCMCmodel.cage.prop.PERM.1, MCMCmodel.cage.prop.PERM.2, MCMCmodel.vial.prop.PERM.1, MCMCmodel.vial.prop.PERM.2, file = "/Users/aneil/Google\ Drive/AgrawalLab/Amardeep/GeneticArchitecture/QGanalysis/MCMCglmm.output.PERMUTE.Rdata")


b.posterior.obs.scale.vial.PERM2 = getResultsFromQGglmm(MCMCmodel.vial.prop.PERM.2, n.avg.counts.vial.f, n.avg.counts.vial.m)
b.posterior.obs.scale.cage.PERM2 = getResultsFromQGglmm(MCMCmodel.cage.prop.PERM.2, n.avg.counts.vial.f, n.avg.counts.vial.m)

{vial.summary.b.PERM2 <-summarizeResultsFromPosteriors(b.posterior.obs.scale.vial.PERM2)
  colnames(vial.summary.b.PERM2) = paste("vial", colnames(vial.summary.b.PERM2), sep = ".")
  cage.summary.b.PERM2 <-summarizeResultsFromPosteriors(b.posterior.obs.scale.cage.PERM2)
  colnames(cage.summary.b.PERM2) = paste("cage", colnames(cage.summary.b.PERM2), sep = ".")
  full.summary.b.PERM2 <-cbind(vial.summary.b.PERM2, cage.summary.b.PERM2)}

full.summary.b[11:13,]
full.summary.b.PERM1[11:13,]
full.summary.b.PERM2[11:13,]

full.summary.b[c(4,5,9,10),]
full.summary.b.PERM2[c(4,5,9,10),]

## genetic variances are much smaller and posterior is strongly skewed (as expected).
plot(as.mcmc(b.posterior.obs.scale.cage.PERM1$f.var.a.obs))
plot(as.mcmc(b.posterior.obs.scale.cage.PERM2$f.var.a.obs))



####### Mutational Burden analysis

invlogit<-function(x) 1/(1 + exp(-x))
logit<- function(p) log(p/(1-p))

QGdata.extended = read.csv("QGdata.extend.csv", header=TRUE)



prior.4 = list(R = list(V = diag(4), nu = 0.002),
               G = list(G1 = list(V=diag(4), nu=1, alpha.mu = rep(0, 4), alpha.V = 1000*diag(4)), G2 = list(V=diag(4), nu=1, alpha.mu = rep(0, 4), alpha.V = 1000*diag(4))))

## remove 2 genotypes with highly deviant numbers of radical amino acid changes
data.for.burden.analysis = na.omit(QGdata.extended[QGdata.extended$DSPR.stock.number !="11297" & QGdata.extended$DSPR.stock.number !="12258",
                                                   c("DSPR.stock.number", "block" , "sex", "mating.regime", 
                                                     "total.wt", "total.wtNOT",
                                                     "most.severe.variants.whole.genome", "zhang.radical.score.polarity.and.volume.whole.genome", "INDELs.whole.genome")])



dim(data.for.burden.analysis)
MCMCmodel.both.prop = MCMCglmm(cbind(total.wt, total.wtNOT) ~ 
                                 sex*mating.regime*(most.severe.variants.whole.genome + zhang.radical.score.polarity.and.volume.whole.genome + INDELs.whole.genome),
                               random = ~us(sex:mating.regime):DSPR.stock.number + idh(sex:mating.regime):block,
                               rcov = ~idh(sex:mating.regime):units,
                               data = data.for.burden.analysis,
                               prior = prior.4,
                               family = "multinomial2",
                               verbose = TRUE,
                               nitt = 5*10^5, burnin = 10000, thin = 100
                               #nitt = 200000, burnin = 10000, thin = 100
                               
)



########
########

n.msv.avg = mean(data.for.burden.analysis$most.severe.variants.whole.genome)
n.rad.avg = mean(data.for.burden.analysis$zhang.radical.score.polarity.and.volume.whole.genome)
n.indel.avg = mean(data.for.burden.analysis$INDELs.whole.genome)

## average fitnesses of each type on observed scale
logit.of.avg.fit.c.f = MCMCmodel.both.prop$Sol[,"(Intercept)"] + n.msv.avg*MCMCmodel.both.prop$Sol[, "most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"INDELs.whole.genome"]
avg.fit.c.f = invlogit(logit.of.avg.fit.c.f)
avg.fit.c.m = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("sexm")] + n.msv.avg*MCMCmodel.both.prop$Sol[, "sexm:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "sexm:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"sexm:INDELs.whole.genome"])
avg.fit.v.f = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("mating.regimevial")] + n.msv.avg*MCMCmodel.both.prop$Sol[, "mating.regimevial:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"mating.regimevial:INDELs.whole.genome"] )
avg.fit.v.m = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("sexm")] + n.msv.avg*MCMCmodel.both.prop$Sol[, "sexm:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "sexm:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"sexm:INDELs.whole.genome"]
                       + MCMCmodel.both.prop$Sol[,c("mating.regimevial")] + n.msv.avg*MCMCmodel.both.prop$Sol[, "mating.regimevial:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"mating.regimevial:INDELs.whole.genome"]
                       +MCMCmodel.both.prop$Sol[,c("sexm:mating.regimevial")] + n.msv.avg*MCMCmodel.both.prop$Sol[, "sexm:mating.regimevial:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "sexm:mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"sexm:mating.regimevial:INDELs.whole.genome"])
colMeans(cbind(avg.fit.c.f, avg.fit.c.m, avg.fit.v.f, avg.fit.v.m))
HPDinterval(as.mcmc(cbind(avg.fit.c.f, avg.fit.c.m, avg.fit.v.f, avg.fit.v.m)))

##### calculate selection on LoFs
## expected fitnesses of each type on observed scale with ONE additional LoF
logit.of.avg.fit.c.f.1msv = MCMCmodel.both.prop$Sol[,"(Intercept)"] + (n.msv.avg + 1)*MCMCmodel.both.prop$Sol[, "most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"INDELs.whole.genome"]
avg.fit.c.f.1msv = invlogit(logit.of.avg.fit.c.f.1msv)
avg.fit.c.m.1msv = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("sexm")] + (n.msv.avg + 1)*MCMCmodel.both.prop$Sol[, "sexm:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "sexm:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"sexm:INDELs.whole.genome"])
avg.fit.v.f.1msv = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("mating.regimevial")] + (n.msv.avg + 1)*MCMCmodel.both.prop$Sol[, "mating.regimevial:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"mating.regimevial:INDELs.whole.genome"] )
avg.fit.v.m.1msv = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("sexm")] + (n.msv.avg + 1)*MCMCmodel.both.prop$Sol[, "sexm:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "sexm:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"sexm:INDELs.whole.genome"]
                            + MCMCmodel.both.prop$Sol[,c("mating.regimevial")] + (n.msv.avg + 1)*MCMCmodel.both.prop$Sol[, "mating.regimevial:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"mating.regimevial:INDELs.whole.genome"]
                            +MCMCmodel.both.prop$Sol[,c("sexm:mating.regimevial")] + (n.msv.avg + 1)*MCMCmodel.both.prop$Sol[, "sexm:mating.regimevial:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "sexm:mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"sexm:mating.regimevial:INDELs.whole.genome"])

sel.c.f.1msv = (avg.fit.c.f.1msv - avg.fit.c.f)/avg.fit.c.f
sel.c.m.1msv = (avg.fit.c.m.1msv - avg.fit.c.m)/avg.fit.c.m
sel.v.f.1msv = (avg.fit.v.f.1msv - avg.fit.v.f)/avg.fit.v.f
sel.v.m.1msv = (avg.fit.v.m.1msv - avg.fit.v.m)/avg.fit.v.m

sel.avg.f.1msv = (sel.c.f.1msv + sel.v.f.1msv)/2
sel.avg.m.1msv = (sel.c.m.1msv + sel.v.m.1msv)/2
sex.diff.sel.avg.1msv = sel.avg.f.1msv - sel.avg.m.1msv
sel.c.avg.1msv = (sel.c.f.1msv + sel.c.m.1msv)/2
sel.v.avg.1msv = (sel.v.f.1msv + sel.v.m.1msv)/2
env.diff.sel.avg.1msv = sel.c.avg.1msv - sel.v.avg.1msv
sel.avg.1msv = (sel.avg.f.1msv + sel.avg.m.1msv)/2

{results.msv = as.data.frame(HPDinterval(as.mcmc(cbind(sel.c.f.1msv, sel.c.m.1msv, sel.v.f.1msv, sel.v.m.1msv, sel.avg.f.1msv, sel.avg.m.1msv, sex.diff.sel.avg.1msv, sel.c.avg.1msv, sel.v.avg.1msv, env.diff.sel.avg.1msv, sel.avg.1msv))))
  results.msv$post.mean = colMeans(as.mcmc(cbind(sel.c.f.1msv, sel.c.m.1msv, sel.v.f.1msv, sel.v.m.1msv, sel.avg.f.1msv, sel.avg.m.1msv, sex.diff.sel.avg.1msv, sel.c.avg.1msv, sel.v.avg.1msv, env.diff.sel.avg.1msv, sel.avg.1msv)))
  results.msv$label = rownames(results.msv)
  results.msv$placeholder = NA
  results.msv = results.msv[c(11, 5:10, 1:4),c(4, 5, 3,1,2)]}
results.msv

##### calculate selection on radical amino acid variants
logit.of.avg.fit.c.f.1rad = MCMCmodel.both.prop$Sol[,"(Intercept)"] + n.msv.avg*MCMCmodel.both.prop$Sol[, "most.severe.variants.whole.genome"] + (n.rad.avg + 1)*MCMCmodel.both.prop$Sol[, "zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"INDELs.whole.genome"]
avg.fit.c.f.1rad = invlogit(logit.of.avg.fit.c.f.1rad)
avg.fit.c.m.1rad = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("sexm")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "sexm:most.severe.variants.whole.genome"] + (n.rad.avg + 1)*MCMCmodel.both.prop$Sol[, "sexm:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"sexm:INDELs.whole.genome"])
avg.fit.v.f.1rad = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("mating.regimevial")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "mating.regimevial:most.severe.variants.whole.genome"] + (n.rad.avg + 1)*MCMCmodel.both.prop$Sol[, "mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"mating.regimevial:INDELs.whole.genome"] )

avg.fit.v.m.1rad = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("sexm")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "sexm:most.severe.variants.whole.genome"] + (n.rad.avg+1)*MCMCmodel.both.prop$Sol[, "sexm:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"sexm:INDELs.whole.genome"]
                            + MCMCmodel.both.prop$Sol[,c("mating.regimevial")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "mating.regimevial:most.severe.variants.whole.genome"] + (n.rad.avg+1)*MCMCmodel.both.prop$Sol[, "mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"mating.regimevial:INDELs.whole.genome"]
                            +MCMCmodel.both.prop$Sol[,c("sexm:mating.regimevial")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "sexm:mating.regimevial:most.severe.variants.whole.genome"] + (n.rad.avg+1)*MCMCmodel.both.prop$Sol[, "sexm:mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + n.indel.avg*MCMCmodel.both.prop$Sol[,"sexm:mating.regimevial:INDELs.whole.genome"])

sel.c.f.1rad = (avg.fit.c.f.1rad - avg.fit.c.f)/avg.fit.c.f
sel.c.m.1rad = (avg.fit.c.m.1rad - avg.fit.c.m)/avg.fit.c.m
sel.v.f.1rad = (avg.fit.v.f.1rad - avg.fit.v.f)/avg.fit.v.f
sel.v.m.1rad = (avg.fit.v.m.1rad - avg.fit.v.m)/avg.fit.v.m

sel.avg.f.1rad = (sel.c.f.1rad + sel.v.f.1rad)/2
sel.avg.m.1rad = (sel.c.m.1rad + sel.v.m.1rad)/2
sex.diff.sel.avg.1rad = sel.avg.f.1rad - sel.avg.m.1rad
sel.c.avg.1rad = (sel.c.f.1rad + sel.c.m.1rad)/2
sel.v.avg.1rad = (sel.v.f.1rad + sel.v.m.1rad)/2
env.diff.sel.avg.1rad = sel.c.avg.1rad - sel.v.avg.1rad
sel.avg.1rad = (sel.avg.f.1rad + sel.avg.m.1rad)/2

{results.rad = as.data.frame(HPDinterval(as.mcmc(cbind(sel.c.f.1rad, sel.c.m.1rad, sel.v.f.1rad, sel.v.m.1rad, sel.avg.f.1rad, sel.avg.m.1rad, sex.diff.sel.avg.1rad, sel.c.avg.1rad, sel.v.avg.1rad, env.diff.sel.avg.1rad, sel.avg.1rad))))
  results.rad$post.mean = colMeans(as.mcmc(cbind(sel.c.f.1rad, sel.c.m.1rad, sel.v.f.1rad, sel.v.m.1rad, sel.avg.f.1rad, sel.avg.m.1rad, sex.diff.sel.avg.1rad, sel.c.avg.1rad, sel.v.avg.1rad, env.diff.sel.avg.1rad, sel.avg.1rad)))
  results.rad$label = rownames(results.rad)
  results.rad$placeholder = NA
  results.rad = results.rad[c(11, 5:10, 1:4),c(4, 5, 3,1,2)]}
results.rad

##### calculate selection on indels
logit.of.avg.fit.c.f.1indel = MCMCmodel.both.prop$Sol[,"(Intercept)"] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "zhang.radical.score.polarity.and.volume.whole.genome"] + (n.indel.avg+1)*MCMCmodel.both.prop$Sol[,"INDELs.whole.genome"]
avg.fit.c.f.1indel = invlogit(logit.of.avg.fit.c.f.1indel)
avg.fit.c.m.1indel = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("sexm")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "sexm:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "sexm:zhang.radical.score.polarity.and.volume.whole.genome"] + (n.indel.avg+1)*MCMCmodel.both.prop$Sol[,"sexm:INDELs.whole.genome"])
avg.fit.v.f.1indel = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("mating.regimevial")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "mating.regimevial:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + (n.indel.avg+1)*MCMCmodel.both.prop$Sol[,"mating.regimevial:INDELs.whole.genome"] )

avg.fit.v.m.1indel = invlogit(logit.of.avg.fit.c.f + MCMCmodel.both.prop$Sol[,c("sexm")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "sexm:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "sexm:zhang.radical.score.polarity.and.volume.whole.genome"] + (n.indel.avg+1)*MCMCmodel.both.prop$Sol[,"sexm:INDELs.whole.genome"]
                              + MCMCmodel.both.prop$Sol[,c("mating.regimevial")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "mating.regimevial:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + (n.indel.avg+1)*MCMCmodel.both.prop$Sol[,"mating.regimevial:INDELs.whole.genome"]
                              +MCMCmodel.both.prop$Sol[,c("sexm:mating.regimevial")] + (n.msv.avg)*MCMCmodel.both.prop$Sol[, "sexm:mating.regimevial:most.severe.variants.whole.genome"] + n.rad.avg*MCMCmodel.both.prop$Sol[, "sexm:mating.regimevial:zhang.radical.score.polarity.and.volume.whole.genome"] + (n.indel.avg+1)*MCMCmodel.both.prop$Sol[,"sexm:mating.regimevial:INDELs.whole.genome"])
colMeans(cbind(avg.fit.c.f.1indel, avg.fit.c.m.1indel, avg.fit.v.f.1indel, avg.fit.v.m.1indel))
HPDinterval(as.mcmc(cbind(avg.fit.c.f.1indel, avg.fit.c.m.1indel, avg.fit.v.f.1indel, avg.fit.v.m.1indel)))

sel.c.f.1indel = (avg.fit.c.f.1indel - avg.fit.c.f)/avg.fit.c.f
sel.c.m.1indel = (avg.fit.c.m.1indel - avg.fit.c.m)/avg.fit.c.m
sel.v.f.1indel = (avg.fit.v.f.1indel - avg.fit.v.f)/avg.fit.v.f
sel.v.m.1indel = (avg.fit.v.m.1indel - avg.fit.v.m)/avg.fit.v.m
colMeans(cbind(sel.c.f.1indel, sel.c.m.1indel, sel.v.f.1indel, sel.v.m.1indel))

sel.avg.f.1indel = (sel.c.f.1indel + sel.v.f.1indel)/2
sel.avg.m.1indel = (sel.c.m.1indel + sel.v.m.1indel)/2
sex.diff.sel.avg.1indel = sel.avg.f.1indel - sel.avg.m.1indel
sel.c.avg.1indel = (sel.c.f.1indel + sel.c.m.1indel)/2
sel.v.avg.1indel = (sel.v.f.1indel + sel.v.m.1indel)/2
env.diff.sel.avg.1indel = sel.c.avg.1indel - sel.v.avg.1indel
sel.avg.1indel = (sel.avg.f.1indel + sel.avg.m.1indel)/2

{results.indel = as.data.frame(HPDinterval(as.mcmc(cbind(sel.c.f.1indel, sel.c.m.1indel, sel.v.f.1indel, sel.v.m.1indel, sel.avg.f.1indel, sel.avg.m.1indel, sex.diff.sel.avg.1indel, sel.c.avg.1indel, sel.v.avg.1indel, env.diff.sel.avg.1indel, sel.avg.1indel))))
  results.indel$post.mean = colMeans(as.mcmc(cbind(sel.c.f.1indel, sel.c.m.1indel, sel.v.f.1indel, sel.v.m.1indel, sel.avg.f.1indel, sel.avg.m.1indel, sex.diff.sel.avg.1indel, sel.c.avg.1indel, sel.v.avg.1indel, env.diff.sel.avg.1indel, sel.avg.1indel)))
  results.indel$label = rownames(results.indel)
  results.indel$placeholder = NA
  results.indel = results.indel[c(11, 5:10, 1:4),c(4, 5, 3,1,2)]}

MakeSimpleFigFromSummaryTable.For.s<-function(summaryTable, x.varLabel, y.varLabel){
  xvals = 1:nrow(summaryTable)
  ymin = min(summaryTable[,4])
  ymax = max(summaryTable[,5])
  yrange = ymax - ymin
  plot(xvals,
       summaryTable[,3], 
       ylim = c(ymin - 0.1*yrange, ymax + 0.1*yrange),
       xlab = x.varLabel, ylab = y.varLabel, pch = 19, xaxt = 'n', col = c("black", "grey", "grey", "black", "grey", "grey", "black"))
  axis(1, at = 1:nrow(summaryTable), labels = summaryTable[,1], las = 2, cex.axis = 0.8)
  # abline(v = 1.7, lwd = 3)
  # abline(v = 2.7, lwd = 0.5); abline(v = 3.7, lwd = 0.5); 
  arrows(x0=xvals, y0=summaryTable[,4],
         x1 = xvals, y1=summaryTable[,5], 
         code=3, angle = 90, length = 0.05,
         lwd=2,col = c("black", "grey", "grey", "black", "grey", "grey", "black"))
}

MakeSimpleFigFromSummaryTable.For.s(results.msv[1:7,], "msv", ""); abline(h = 0, lty = 2)
MakeSimpleFigFromSummaryTable.For.s(results.rad[1:7,], "rad", ""); abline(h = 0, lty = 2)
MakeSimpleFigFromSummaryTable.For.s(results.indel[1:7,], "indels", ""); abline(h = 0, lty = 2)

########
######## 
## repeating analysis using alternative parameterization gives same results (as expected)
data.for.burden.analysis$msv = scale(data.for.burden.analysis$most.severe.variants.whole.genome, center = T, scale = F)
data.for.burden.analysis$rad = scale(data.for.burden.analysis$zhang.radical.score.polarity.and.volume.whole.genome, center = T, scale = F)
data.for.burden.analysis$indels = scale(data.for.burden.analysis$INDELs.whole.genome, center = T, scale = F)

data.for.burden.analysis$sexADJ = as.numeric(sapply(data.for.burden.analysis$sex, function(x) ifelse(x=="f", -0.5, 0.5)))
data.for.burden.analysis$envADJ = as.numeric(sapply(data.for.burden.analysis$mating.regime, function(x) ifelse(x=="cage", -0.5, 0.5)))
data.for.burden.analysis$intADJ.sexXenv = data.for.burden.analysis$sexADJ*data.for.burden.analysis$envADJ


MCMCmodel.both.propV2 = MCMCglmm(cbind(total.wt, total.wtNOT) ~ 
                                   (sexADJ + envADJ + intADJ.sexXenv)*(msv + rad + indels),
                                 random = ~us(sex:mating.regime):DSPR.stock.number + idh(sex:mating.regime):block,
                                 rcov = ~idh(sex:mating.regime):units,
                                 data = data.for.burden.analysis,
                                 prior = prior.4,
                                 family = "multinomial2",
                                 verbose = TRUE,
                                 nitt = 5*10^5, burnin = 10000, thin = 100
                                 #nitt = 200000, burnin = 10000, thin = 100
                                 
)

save(MCMCmodel.both.prop, MCMCmodel.both.propV2, file = "/Users/aneil/Google\ Drive/AgrawalLab/Amardeep/GeneticArchitecture/QGanalysis/MCMC_BurdenAnalyses.Rdata")



AP.avg.fit.c.f = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] )
AP.avg.fit.c.m = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] )
AP.avg.fit.v.f = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] )
AP.avg.fit.v.m = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] )

AP.avg.fit.c.f.1msv = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] 
                               + MCMCmodel.both.propV2$Sol[,"msv"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:msv"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:msv"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:msv"])
AP.avg.fit.c.m.1msv = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"]
                               + MCMCmodel.both.propV2$Sol[,"msv"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:msv"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:msv"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:msv"])
AP.avg.fit.v.f.1msv = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] 
                               + MCMCmodel.both.propV2$Sol[,"msv"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:msv"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:msv"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:msv"])
AP.avg.fit.v.m.1msv = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] 
                               + MCMCmodel.both.propV2$Sol[,"msv"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:msv"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:msv"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:msv"])



AP.sel.c.f.1msv = (AP.avg.fit.c.f.1msv - AP.avg.fit.c.f)/AP.avg.fit.c.f
AP.sel.c.m.1msv = (AP.avg.fit.c.m.1msv - AP.avg.fit.c.m)/AP.avg.fit.c.m
AP.sel.v.f.1msv = (AP.avg.fit.v.f.1msv - AP.avg.fit.v.f)/AP.avg.fit.v.f
AP.sel.v.m.1msv = (AP.avg.fit.v.m.1msv - AP.avg.fit.v.m)/AP.avg.fit.v.m
colMeans(cbind(AP.sel.c.f.1msv, AP.sel.c.m.1msv, AP.sel.v.f.1msv, AP.sel.v.m.1msv))

AP.sel.avg.f.1msv = (AP.sel.c.f.1msv + AP.sel.v.f.1msv)/2
AP.sel.avg.m.1msv = (AP.sel.c.m.1msv + AP.sel.v.m.1msv)/2
AP.sex.diff.sel.avg.1msv = AP.sel.avg.f.1msv - AP.sel.avg.m.1msv
AP.sel.c.avg.1msv = (AP.sel.c.f.1msv + AP.sel.c.m.1msv)/2
AP.sel.v.avg.1msv = (AP.sel.v.f.1msv + AP.sel.v.m.1msv)/2
AP.env.diff.sel.avg.1msv = AP.sel.c.avg.1msv - AP.sel.v.avg.1msv
AP.sel.avg.1msv = (AP.sel.avg.f.1msv + AP.sel.avg.m.1msv)/2


##### radical AA selection
AP.avg.fit.c.f.1rad = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] 
                               + MCMCmodel.both.propV2$Sol[,"rad"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:rad"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:rad"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:rad"])
AP.avg.fit.c.m.1rad = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"]
                               + MCMCmodel.both.propV2$Sol[,"rad"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:rad"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:rad"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:rad"])
AP.avg.fit.v.f.1rad = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] 
                               + MCMCmodel.both.propV2$Sol[,"rad"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:rad"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:rad"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:rad"])
AP.avg.fit.v.m.1rad = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] 
                               + MCMCmodel.both.propV2$Sol[,"rad"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:rad"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:rad"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:rad"])

AP.sel.c.f.1rad = (AP.avg.fit.c.f.1rad - AP.avg.fit.c.f)/AP.avg.fit.c.f
AP.sel.c.m.1rad = (AP.avg.fit.c.m.1rad - AP.avg.fit.c.m)/AP.avg.fit.c.m
AP.sel.v.f.1rad = (AP.avg.fit.v.f.1rad - AP.avg.fit.v.f)/AP.avg.fit.v.f
AP.sel.v.m.1rad = (AP.avg.fit.v.m.1rad - AP.avg.fit.v.m)/AP.avg.fit.v.m

AP.sel.avg.f.1rad = (AP.sel.c.f.1rad + AP.sel.v.f.1rad)/2
AP.sel.avg.m.1rad = (AP.sel.c.m.1rad + AP.sel.v.m.1rad)/2
AP.sex.diff.sel.avg.1rad = AP.sel.avg.f.1rad - AP.sel.avg.m.1rad
AP.sel.c.avg.1rad = (AP.sel.c.f.1rad + AP.sel.c.m.1rad)/2
AP.sel.v.avg.1rad = (AP.sel.v.f.1rad + AP.sel.v.m.1rad)/2
AP.env.diff.sel.avg.1rad = AP.sel.c.avg.1rad - AP.sel.v.avg.1rad
AP.sel.avg.1rad = (AP.sel.avg.f.1rad + AP.sel.avg.m.1rad)/2



##### indel selection
AP.avg.fit.c.f.1indel = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] 
                                 + MCMCmodel.both.propV2$Sol[,"indels"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:indels"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:indels"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:indels"])
AP.avg.fit.c.m.1indel = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"]
                                 + MCMCmodel.both.propV2$Sol[,"indels"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:indels"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:indels"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:indels"])
AP.avg.fit.v.f.1indel = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] 
                                 + MCMCmodel.both.propV2$Sol[,"indels"] + (-1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:indels"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:indels"] + (-1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:indels"])
AP.avg.fit.v.m.1indel = invlogit(MCMCmodel.both.propV2$Sol[,"(Intercept)"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv"] 
                                 + MCMCmodel.both.propV2$Sol[,"indels"] + (1/2)*MCMCmodel.both.propV2$Sol[,"sexADJ:indels"] + (1/2)*MCMCmodel.both.propV2$Sol[,"envADJ:indels"] + (1/4)*MCMCmodel.both.propV2$Sol[,"intADJ.sexXenv:indels"])

AP.sel.c.f.1indel = (AP.avg.fit.c.f.1indel - AP.avg.fit.c.f)/AP.avg.fit.c.f
AP.sel.c.m.1indel = (AP.avg.fit.c.m.1indel - AP.avg.fit.c.m)/AP.avg.fit.c.m
AP.sel.v.f.1indel = (AP.avg.fit.v.f.1indel - AP.avg.fit.v.f)/AP.avg.fit.v.f
AP.sel.v.m.1indel = (AP.avg.fit.v.m.1indel - AP.avg.fit.v.m)/AP.avg.fit.v.m

AP.sel.avg.f.1indel = (AP.sel.c.f.1indel + AP.sel.v.f.1indel)/2
AP.sel.avg.m.1indel = (AP.sel.c.m.1indel + AP.sel.v.m.1indel)/2
AP.sex.diff.sel.avg.1indel = AP.sel.avg.f.1indel - AP.sel.avg.m.1indel
AP.sel.c.avg.1indel = (AP.sel.c.f.1indel + AP.sel.c.m.1indel)/2
AP.sel.v.avg.1indel = (AP.sel.v.f.1indel + AP.sel.v.m.1indel)/2
AP.env.diff.sel.avg.1indel = AP.sel.c.avg.1indel - AP.sel.v.avg.1indel
AP.sel.avg.1indel = (AP.sel.avg.f.1indel + AP.sel.avg.m.1indel)/2











