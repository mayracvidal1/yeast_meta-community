# This script includes all the statistical analyses for the MS titled ""
#
#

# required packages
require(ggplot2) # version 3.5.1
require(viridis) # version 0.6.5
require(gridExtra) # version 2.3
require(lme4) # version 1.1-35.5
require(MCMCglmm) # version 2.36
require(betapart) # version 1.6
require(MASS) # version 7.3-61
require(hillR) # version 0.5.2
require(glmmTMB) # version 1.1.10
library(tidyverse) # version 2.0.0

# 1. Diversity estimation

# 1.1 alpha diversity
#alpha diversity based on Hill numbers  species richness (q=0)
#Shannon diversity (q=1, the exponential of Shannon entropy) and
#Simpson diversity (q=2, the inverse of Simpson concentration).

# for exploited community
comm1<-read.csv("comm_exploited.csv",row.names = 1)
comm1_meta<-comm1[comm1$meta=="yes",]
comm1_meta<-comm1_meta[,3:11]
comm1_meta<-comm1_meta*10
q0<-hill_taxa(comm1_meta, q = 0)
q1<-hill_taxa(comm1_meta, q = 1)
q2<-hill_taxa(comm1_meta, q = 2)
alpha_meta1<-data.frame("q0"=q0,"q1"=q1,"q2"=q2)
write.csv(alpha_meta1,file="alpha_meta1.csv")
comm1_non_meta<-comm1[comm1$meta=="no",]
comm1_non_meta<-comm1_non_meta[,3:11]
comm1_non_meta<-comm1_non_meta*10
q0<-hill_taxa(comm1_non_meta, q = 0)
q1<-hill_taxa(comm1_non_meta, q = 1)
q2<-hill_taxa(comm1_non_meta, q = 2)
alpha_non_meta1<-data.frame("q0"=q0,"q1"=q1,"q2"=q2)

# for mutualists-only plate
comm2<-read.csv("comm_mutualists_only.csv",row.names = 1)
comm2_meta<-comm2[comm2$meta=="yes",]
comm2_meta<-comm2_meta[,3:10]
comm2_meta<-comm2_meta*10
q0<-hill_taxa(comm2_meta, q = 0)
q1<-hill_taxa(comm2_meta, q = 1)
q2<-hill_taxa(comm2_meta, q = 2)
alpha_meta2<-data.frame("q0"=q0,"q1"=q1,"q2"=q2)
write.csv(alpha_meta1,file="alpha_meta2.csv")
comm2_non_meta<-comm2[comm2$meta=="no",]
comm2_non_meta<-comm2_non_meta[,3:10]
comm2_non_meta<-comm2_non_meta*10
q0<-hill_taxa(comm2_non_meta, q = 0)
q1<-hill_taxa(comm2_non_meta, q = 1)
q2<-hill_taxa(comm2_non_meta, q = 2)
alpha_non_meta2<-data.frame("q0"=q0,"q1"=q1,"q2"=q2)

# 1.2 Beta diversity
# only calculate beta diversity for mutualists-only community
for (i in c("yes","no")){
  comm<-comm2[comm_mcv$meta==i,]
  for (j in 2:7){
    comm_x<-comm[comm$week==j,]
    comm_x<-comm_x[,4:11]*10
    tbd_x<-as.matrix(beta.pair.abund(comm_x,index.family = "bray")$beta.bray)
    #write.csv(tbd_x,file=paste0("tbd_",i,j,".csv"))
    tbd_x[!lower.tri(tbd_x)]<-NA
    dim(tbd_x)<-c(dim(tbd_x)[1]*dim(tbd_x)[2],1)
    tbd_x<-tbd_x[complete.cases(tbd_x),]
    write.csv(tbd_x,file=paste0("tbd_",i,j,".csv"))
  }
}


# 2. MCMCglmm models predicting the strain abundance by the type of mutualism,
#   type of mutualists, or yeast strains (Table 1)

# MCMCglmm with ordinal family
glm_abundance<-read.csv("abudance_week2.csv")
glm_abundance_adeop<-glm_abundance[glm_abundance$strain_type=="Adeop",]
glm_abundance_lysop<-glm_abundance[glm_abundance$strain_type=="Lysop",]
#
#set prior
a <- 1000
# Number of interations
nitt <- 2400000
# Length of burni
burnin <- 300000
# Amount of thinning
thin <- 1000

# Model a
prior_a <- list(R = list(V = diag(1), nu = 0.002,fix=1),
                G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                         G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
mcmc_abundance_adeop_lysop <- MCMCglmm(abundance ~ factor(plate)+factor(strain_type), random = ~comm+strains,
                                       data = glm_abundance, prior = prior_a, nitt = nitt,burnin=burnin,thin=thin)

plot(mcmc_abundance_adeop_lysop$Sol)
summary(mcmc_abundance_adeop_lysop)

# Model b
prior_b <- list(R = list(V = diag(1), nu = 0.002,fix=1),
                G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
mcmc_abundance_ade <- MCMCglmm(abundance ~ factor(plate)+factor(strains), random = ~comm,
                               data = glm_abundance_adeop, prior = prior_b, nitt = nitt,burnin=burnin,thin=thin)
plot(mcmc_abundance_ade$Sol)
summary(mcmc_abundance_ade)


#Model c
mcmc_abundance_lys <- MCMCglmm(abundance ~ factor(plate)+factor(strains), random = ~comm,
                               data = glm_abundance_lysop, prior = prior_b, nitt = nitt,burnin=burnin,thin=thin)
plot(mcmc_abundance_lys$Sol)
summary(mcmc_abundance_lys)

# 3. MCMCglmm for alpha diversity (hill number q=1)
# loading alpha diversity data
mcmc_alpha<-read.csv("alpha_diversity.csv") %>%
  filter(q==1)
#set prior
a <- 1000
prior3 <- list(R = list(V = diag(1), nu = 0.002,fix=1),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))
# Number of interations
nitt <- 2400000
# Length of burni
burnin <- 300000
# Amount of thinning
thin <- 1000

mcmc_alpha_model <- MCMCglmm(d_value ~ factor(meta_n)+factor(cheater), random = ~comm_id+plate_id,
                             data = mcmc_alpha, prior = prior3, nitt = nitt,burnin=burnin,thin=thin)

plot(mcmc_alpha_model$Sol)
summary(mcmc_alpha_model)

# 4.rescue effect of dispersal in mutualists-only vs. exploited communities
#glmm

# Fit GLMM with plate_id as random effect
intro<-read.csv("proportion_recovered_strains.csv")
model <- glmmTMB(p_strong ~ as.factor(community_type_id) + (1 | plate_id), data = intro)
summary(model)


# 4. glmm for beta diversity
beta_d<-read.csv("beta_diversity.csv")

# for dispersal_allowed communities
beta_d_meta<-beta_d%>%
  filter(dispersal_allowed=="yes") %>%
  filter(week >3) # analyse beta diversity pattern from week3.

ggplot(data=beta_d_meta,aes(x=as.factor(week),y=beta_diversity,fill=plate_id))+geom_boxplot()

model <- glmmTMB(beta_diversity ~ week + (1 | plate_id+pairwise_id), data = beta_d_meta)
summary(model)

# for isolated communities
beta_d_no_meta<-beta_d%>%
  filter(dispersal_allowed=="no")

ggplot(data=beta_d_no_meta,aes(x=as.factor(week),y=beta_diversity,fill=plate_id))+geom_boxplot()
model <- glmmTMB(beta_diversity ~ week + (1 | plate_id+pairwise_id), data = beta_d_no_meta)
summary(model)
