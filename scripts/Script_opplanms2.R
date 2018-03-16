# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
# A: in stocking area), B:below stocking area
phiA <- c(.4, .9, .75, .9)
phiB <- phiA
psiAB <- 0.05
psiBA <- 0.001
pA <- 0.60
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

# Execute simulation function
sim <- simul_ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH
f <- apply(CH, 1, get_first)
  
# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3 = not seen
rCH <- CH          # Recoded CH
rCH[rCH==0] <- 3
  
jags.data <- list(y = rCH, 
                  f = f, 
                  n.occasions = dim(rCH)[2], 
                  nind = dim(rCH)[1], 
                  z = known_state_ms(rCH, 3))

# Initial values 
inits <- function(x = rCH){list(phiA = runif(n.occasions - 1, 0, 1), 
                         phiB = runif(n.occasions - 1, 0, 1),
                         mean.psi = runif(2, 0, .2), 
                         mean.p = runif(2, .1, .7), 
                         z = ms_init_z(x, f))}  


# Parameters monitored
parameters <- c("phiA", "phiB", "mean.psi", "mean.p")

# MCMC settings
ni <- 4000
nt <- 3
nb <- 2000
nc <- 3

# Call JAGS from R
ms2 <- jagsUI::jags(jags.data, 
                    inits, 
                    parameters, 
                    ".\\models\\model_ms2.txt", 
                    n.chains = nc, 
                    n.thin = nt, 
                    n.iter = ni, 
                    n.burnin = nb, 
                    parallel = TRUE)

print(ms2, digits = 3)