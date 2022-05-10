#MCMCglmm model from Palmer Foster 2022

#setwd("C:/Users/jdpal/Documents/R/ToxinStats")
library(phytools)
library(phangorn)
library(dispRity)
library(geiger)
library(lme4)
library(phyr)
library(dplyr)
library(INLA)
library(ggridges)
pacman::p_load(pacman, ggplot2, MCMCglmm, ape, coda, phytools, phangorn, dispRity,phyr,dplyr)

#Load phylogeny
filo <- read.nexus("Phylo268_Single.nex")
#Load spectrum and regulation data.  Single_DD0.csv for initial model.  Single_Strict.csv for strict model.
dader <- read.csv("Single_Strict.csv", header = TRUE)

#Convert phylogeny to chronogram
filo1 <- chronos(filo, lambda = 1, model = "correlated", quiet = TRUE,
                 calibration = makeChronosCalib(filo),
                 control = chronos.control())
#root midpoint, remove zeros, set inverse matrix
filo1$node.label = NULL
filo1 <- midpoint.root(filo1)
filo1 <- remove.zero.brlen(filo1)
inv.filo1 <- inverseA(filo1)

#model running specs
nitt <- 1000000
burnin <- 5000
thin <- 50

#Define priors
prior <- list(B=list(mu=c(0,0,0), V=gelman.prior(~as.factor(Spectr), 
              data=dader, scale=1+1+pi^2/3)), R=list(V=1, fix=1), 
              G=list(G1=list(V=diag(1)*0.1, nu=1)))

#Run Model
model_simple<-MCMCglmm(DDR~as.factor(Spectr),random=~Bacteria,
                       family="categorical",ginverse=list(Bacteria=inv.filo1$Ainv),verbose=FALSE,prior=prior,
                       data=dader, scale = F,slice = T, nitt = nitt, burnin = burnin, thin = thin)

#output results
summary(model_simple)
posterior.mode(model_simple$Sol)


