library(DredgeRainbows)

# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
# A: in stocking area), B:below stocking area
#Run section between hash marks to create parameters for each simulation
######
phiA <- c(.4, .9, .75, .9)
phiB <- phiA /4
psiAB <- 0.05
psiBA <- 0.01
pA <- 0.20
pB <- pA * (20/10.3)/(60/(6.7+9.4+6.4))
n.occasions <- 5 #Fall 18, Spring and Fall 19, Spring and Fall 20
n.states <- 3
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- c(1500, 1500, 0, 1500, 0)
marked[,2] <- marked[,3] <- rep(0, n.occasions)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time
# 1. State process matrix
PSI.STATE <- array(NA, dim=c(n.states, n.states, sum(marked), n.occasions-1))
for (i in 1:sum(marked)){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phiA[t]*(1-psiAB), phiA[t]*psiAB,     1-phiA[t],
      phiB[t]*psiBA,     phiB[t]*(1-psiBA), 1-phiB[t],
      0,              0,              1       ), nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, sum(marked), n.occasions-1))
for (i in 1:sum(marked)){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  1-pA,
      0,  pB, 1-pB,
      0,  0,  1       ), nrow = n.states, byrow = TRUE)
  } #t
} #i
######

# Initial values
inits <- function(x = rCH, f_0 = f){list(phiA = runif(n.occasions - 1, 0, 1),
                         phiB = runif(n.occasions - 1, 0, 1),
                         mean.psi = runif(2, 0, .2),
                         mean.p = runif(2, .1, .7),
                         z = ms_init_z(x, f_0))}


# Parameters monitored
parameters <- c("phiA", "phiB", "mean.psi", "mean.p")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 2000
nc <- 3

#####vary phiAB#####

#####phiB = phiA#####
#use sim_dat20pct
mean(sapply(sim_post20pct, function(x)  x$mean$mean.psi[1]))
quantile(sapply(sim_post20pct, function(x)  (x$q97.5$mean.psi[1] - x$q2.5$mean.psi[1])), probs = c(0, .1, .5, .9, 1))
mean(sapply(sim_post20pct, function(x)  x$q97.5$mean.psi[1] >= 0.05))

#####simulate dataset: phiA/2#####
#phiA <- c(.4, .9, .75, .9), phiB <- phiA / 2, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.20, pB <- pA * .73
#sim_dathalf <- replicate(30, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE)

#####run simulations: phiA/2#####
sim_posthalf <-
  lapply(sim_dathalf, function(x){
    jagsUI::jags(x,
                 inits = list(inits(x$y, x$f), inits(x$y, x$f), inits(x$y, x$f)),
                 parameters,
                 ".\\models\\model_ms2.txt",
                 n.chains = nc,
                 n.thin = nt,
                 n.iter = ni,
                 n.burnin = nb,
                 parallel = TRUE)
  })
mean(sapply(sim_posthalf, function(x)  x$mean$mean.psi[1]))
quantile(sapply(sim_posthalf, function(x)  (x$q97.5$mean.psi[1] - x$q2.5$mean.psi[1])), probs = c(0, .1, .5, .9, 1))
mean(sapply(sim_posthalf, function(x)  x$q97.5$mean.psi[1] >= 0.05))

#####simulate dataset: phiA/4#####
#phiA <- c(.4, .9, .75, .9), phiB <- phiA / 4, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.20, pB <- pA * .73
#sim_datquater <- replicate(30, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE)

#####run simulations: phiA/4#####
sim_postquater <-
  lapply(sim_datquater, function(x){
    jagsUI::jags(x,
                 inits = list(inits(x$y, x$f), inits(x$y, x$f), inits(x$y, x$f)),
                 parameters,
                 ".\\models\\model_ms2.txt",
                 n.chains = nc,
                 n.thin = nt,
                 n.iter = ni,
                 n.burnin = nb,
                 parallel = TRUE)
    })
mean(sapply(sim_postquater, function(x)  x$mean$mean.psi[1]))
quantile(sapply(sim_postquater, function(x)  (x$q97.5$mean.psi[1] - x$q2.5$mean.psi[1])), probs = c(0, .1, .5, .9, 1))
mean(sapply(sim_postquater, function(x)  x$q97.5$mean.psi[1] >= 0.05))


#####save phiAB sims#####
sims_phiAB <- list(sim_dathalf, sim_posthalf,
                   sim_datquater, sim_postquater)
#saveRDS(sims_phiAB, ".\\scripts\\sim_phiAB.rds")

######vary pA#####

#####simulate dataset: pA = .3#####
#phiA <- c(.4, .9, .75, .9), phiB <- phiA, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.30, pB <- pA *.73
#sim_dat30pct <- replicate(30, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE)

#####run simulations: pA = .3#####
sim_post30pct <-
  lapply(sim_dat30pct, function(x){
    jagsUI::jags(x,
                 inits = list(inits(x$y, x$f), inits(x$y, x$f), inits(x$y, x$f)),
                 parameters,
                 ".\\models\\model_ms2.txt",
                 n.chains = nc,
                 n.thin = nt,
                 n.iter = ni,
                 n.burnin = nb,
                 parallel = TRUE)
  })
apply(sapply(sim_post30pct, function(x)  x$mean$phiA), 1, function(x) mean(unlist(x)))
apply(sapply(sim_post30pct, function(x)  (x$q97.5$phiA - x$q2.5$phiA)), 1, function(x) quantile(unlist(x), probs = c(0, .1, .5, .9, 1)))
apply(sapply(sim_post30pct, function(x)  x$q97.5$phiA >= c(.4, .9, .75, .9) & x$q2.5$phiA <= c(.4, .9, .75, .9)), 1, mean)


#####simulate dataset: pA = .2#####
#phiA <- c(.4, .9, .75, .9), phiB <- phiA, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.20, pB <- pA * .73
#sim_dat20pct <- replicate(30, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE)
#example capture history
knitr::kable(table(paste0(sim_dat20pct[[1]]$y[,1], sim_dat20pct[[1]]$y[,2], sim_dat20pct[[1]]$y[,3], sim_dat20pct[[1]]$y[,4], sim_dat20pct[[1]]$y[,5])))

#####run simulations: pA = .2#####
sim_post20pct <-
  lapply(sim_dat20pct, function(x){
    jagsUI::jags(x,
                 inits = list(inits(x$y, x$f), inits(x$y, x$f), inits(x$y, x$f)),
                 parameters,
                 ".\\models\\model_ms2.txt",
                 n.chains = nc,
                 n.thin = nt,
                 n.iter = ni,
                 n.burnin = nb,
                 parallel = TRUE)
  })
apply(sapply(sim_post20pct, function(x)  x$mean$phiA), 1, function(x) mean(unlist(x)))
apply(sapply(sim_post20pct, function(x)  (x$q97.5$phiA - x$q2.5$phiA)), 1, function(x) quantile(unlist(x), probs = c(0, .1, .5, .9, 1)))
apply(sapply(sim_post20pct, function(x)  x$q97.5$phiA >= c(.4, .9, .75, .9) & x$q2.5$phiA <= c(.4, .9, .75, .9)), 1, mean)

#####simulate dataset: pA = .1#####
#phiA <- c(.4, .9, .75, .9), phiB <- phiA, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.10, pB <- pA * .73
#sim_dat10pct <- replicate(30, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE)

#####run simulations: pA = .1#####
sim_post10pct <-
  lapply(sim_dat10pct, function(x){
    jagsUI::jags(x,
                 inits = list(inits(x$y, x$f), inits(x$y, x$f), inits(x$y, x$f)),
                 parameters,
                 ".\\models\\model_ms2.txt",
                 n.chains = nc,
                 n.thin = nt,
                 n.iter = ni,
                 n.burnin = nb,
                 parallel = TRUE)
  })
apply(sapply(sim_post10pct, function(x)  x$mean$phiA), 1, function(x) mean(unlist(x)))
apply(sapply(sim_post10pct, function(x)  (x$q97.5$phiA - x$q2.5$phiA)), 1, function(x) quantile(unlist(x), probs = c(0, .1, .5, .9, 1)))
apply(sapply(sim_post10pct, function(x)  x$q97.5$phiA >= c(.4, .9, .75, .9) & x$q2.5$phiA <= c(.4, .9, .75, .9)), 1, mean)

#####save pA sims#####
sims_pA <- list(sim_dat30pct, sim_post30pct,
                sim_dat20pct, sim_post20pct,
                sim_dat10pct, sim_post10pct)
#saveRDS(sims_pA, ".\\scripts\\sim_pA.rds")

