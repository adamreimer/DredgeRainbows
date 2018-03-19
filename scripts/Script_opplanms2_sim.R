# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
# A: in stocking area), B:below stocking area
phiA <- c(.4, .9, .75, .9)
phiB <- phiA
psiAB <- 0.05
psiBA <- 0.001
pA <- 0.20
pB <- pA / 3
n.occasions <- 5 #Fall 18, Spring and Fall 19, Spring and Fall 20
n.states <- 3
n.obs <- 3
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
marked[,1] <- c(1400, 1400, 0, 1400, 0)
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

jags_data <- function(PSI.STATE, PSI.OBS, marked){
  # Execute simulation function
  sim <- simul_ms(PSI.STATE, PSI.OBS, marked)
  CH <- sim$CH
  f <- apply(CH, 1, get_first)
  
  # Recode CH matrix: note, a 0 is not allowed in WinBUGS!
  # 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
  rCH <- CH          # Recoded CH
  rCH[rCH==0] <- 3

  # Bundle data
  list(y = rCH, 
       f = f, 
       n.occasions = dim(rCH)[2], 
       nind = dim(rCH)[1], 
       z = known_state_ms(rCH, 3))
}

# Initial values 
inits <- function(x = rCH, f_0 = f){list(phiA = runif(n.occasions - 1, 0, 1), 
                         phiB = runif(n.occasions - 1, 0, 1),
                         mean.psi = runif(2, 0, .2), 
                         mean.p = runif(2, .1, .7), 
                         z = ms_init_z(x, f_0))}  


# Parameters monitored
parameters <- c("phiA", "phiB", "mean.psi", "mean.p")

# MCMC settings
ni <- 3000
nt <- 3
nb <- 1500
nc <- 3

#vary phiAB
#sim_dat <- replicate(25, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE) 
#phiA <- c(.4, .9, .75, .9), phiB <- phiA, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.60, pB <- pA /3

knitr::kable(table(paste0(sim_dat[[1]]$y[,1], sim_dat[[1]]$y[,2], sim_dat[[1]]$y[,3], sim_dat[[1]]$y[,4], sim_dat[[1]]$y[,5])))
sim_post <- 
  lapply(sim_dat, function(x){
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
mean(sapply(sim_post, function(x)  x$mean$mean.psi[1]))
quantile(sapply(sim_post, function(x)  (x$q97.5$mean.psi[1] - x$q2.5$mean.psi[1])), probs = c(0, .1, .5, .9, 1))
mean(sapply(sim_post, function(x)  x$q97.5$mean.psi[1] >= 0.05))

#sim_dathalf <- replicate(25, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE)
#phiA <- c(.4, .9, .75, .9), phiB <- phiA / 2, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.60, pB <- pA / 3
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

#sim_datquater <- replicate(25, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE)
#phiA <- c(.4, .9, .75, .9), phiB <- phiA / 4, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.60, pB <- pA / 3
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

# sims_phiAB <- list(sim_dat, sim_post,
#                    sim_dathalf, sim_posthalf,
#                    sim_datquater, sim_postquater)
# saveRDS(sims_phiAB, ".\\scripts\\sim_phiAB.rds")

#vary pA
#sim_dat <- replicate(25, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE) 
#phiA <- c(.4, .9, .75, .9), phiB <- phiA, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.60, pB <- pA / 3
#Can use sim_post for pA=0.6 run
apply(sapply(sim_post, function(x)  x$mean$phiA), 1, function(x) mean(unlist(x)))
apply(sapply(sim_post, function(x)  (x$q97.5$phiA - x$q2.5$phiA)), 1, function(x) quantile(unlist(x), probs = c(0, .1, .5, .9, 1)))
apply(sapply(sim_post, function(x)  x$q97.5$phiA >= c(.4, .9, .75, .9) & x$q2.5$phiA <= c(.4, .9, .75, .9)), 1, mean)

#sim_dat40pct <- replicate(25, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE)
#phiA <- c(.4, .9, .75, .9), phiB <- phiA, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.40, pB <- pA / 3
sim_post40pct <- 
  lapply(sim_dat40pct, function(x){
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
apply(sapply(sim_post40pct, function(x)  x$mean$phiA), 1, function(x) mean(unlist(x)))
apply(sapply(sim_post40pct, function(x)  (x$q97.5$phiA - x$q2.5$phiA)), 1, function(x) quantile(unlist(x), probs = c(0, .1, .5, .9, 1)))
apply(sapply(sim_post40pct, function(x)  x$q97.5$phiA >= c(.4, .9, .75, .9) & x$q2.5$phiA <= c(.4, .9, .75, .9)), 1, mean)

sim_dat20pct <- replicate(25, jags_data(PSI.STATE, PSI.OBS, marked), simplify = FALSE)
#phiA <- c(.4, .9, .75, .9), phiB <- phiA, psiAB <- 0.05, psiBA <- 0.001, pA <- 0.20, pB <- pA / 3
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

sims_pA <- list(sim_dat40pct, sim_post40pct,
                sim_dat20pct, sim_post20pct)
saveRDS(sims_pA, ".\\scripts\\sim_pA.rds")

