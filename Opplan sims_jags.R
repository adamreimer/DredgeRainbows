# Define mean survival, transitions, recapture, as well as number of occasions, states, observations and released individuals
# A: Glacier/Moraine (best habitat above drainage control), B:Crystal (poor habitat above drainage control), C: Moose (below drainage control)
phiA <- 0.75
phiB <- 0.75
phiC <- 0.75
psiAB <- 0.2
psiAC <- 0.05
psiBA <- 0.5
psiBC <- 0.05
psiCA <- 0.05
psiCB <- 0.05
pA <- 0.20
pB <- 0.20
pC <- 0.10
n.occasions <- 5
#Fall 18, Spring and Fall 19, Spring and Fall 20
n.states <- 4
n.obs <- 4
marked <- matrix(NA, ncol = n.states, nrow = n.occasions)
# marked[,1] <- rep(500, n.occasions)
# marked[,2] <- rep(500, n.occasions)
# marked[,3] <- rep(500, n.occasions)
# marked[,1] <- marked[,2] <- marked[,3] <- c(500, rep(0, n.occasions - 1))
marked[,1] <- marked[,2] <- marked[,3] <- c(500, 500, 0, 500, 0)
marked[,4] <- rep(0, n.occasions)

# Define matrices with survival, transition and recapture probabilities
# These are 4-dimensional matrices, with 
# Dimension 1: state of departure
# Dimension 2: state of arrival
# Dimension 3: individual
# Dimension 4: time
# 1. State process matrix
totrel <- sum(marked)*(n.occasions-1)
PSI.STATE <- array(NA, dim=c(n.states, n.states, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.STATE[,,i,t] <- matrix(c(
      phiA*(1-psiAB-psiAC), phiA*psiAB,           phiA*psiAC,           1-phiA,
      phiB*psiBA,           phiB*(1-psiBA-psiBC), phiB*psiBC,           1-phiB,
      phiC*psiCA,           phiC*psiCB,           phiC*(1-psiCA-psiCB), 1-phiC,
      0,                    0,                    0,                    1), nrow = n.states, byrow = TRUE)
  } #t
} #i

# 2.Observation process matrix
PSI.OBS <- array(NA, dim=c(n.states, n.obs, totrel, n.occasions-1))
for (i in 1:totrel){
  for (t in 1:(n.occasions-1)){
    PSI.OBS[,,i,t] <- matrix(c(
      pA, 0,  0,  1-pA,
      0,  pB, 0,  1-pB,
      0,  0,  pC, 1-pC,
      0,  0,  0,  1), nrow = n.states, byrow = TRUE)
  } #t
} #i

# Define function to simulate multistate capture-recapture data
simul.ms <- function(PSI.STATE, PSI.OBS, marked, unobservable = NA){
  # Unobservable: number of state that is unobservable
  n.occasions <- dim(PSI.STATE)[4] + 1
  CH <- CH.TRUE <- matrix(NA, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- matrix(0, ncol = dim(PSI.STATE)[1], nrow = sum(marked))
  g <- colSums(marked)
  for (s in 1:dim(PSI.STATE)[1]){
    if (g[s]==0) next  # To avoid error message if nothing to replace
    mark.occ[(cumsum(g[1:s])-g[s]+1)[s]:cumsum(g[1:s])[s],s] <-
      rep(1:n.occasions, marked[1:n.occasions,s])
  } #s
  for (i in 1:sum(marked)){
    for (s in 1:dim(PSI.STATE)[1]){
      if (mark.occ[i,s]==0) next
      first <- mark.occ[i,s]
      CH[i,first] <- s
      CH.TRUE[i,first] <- s
    } #s
    for (t in (first+1):n.occasions){
      # Multinomial trials for state transitions
      if (first==n.occasions) next
      state <- which(rmultinom(1, 1, PSI.STATE[CH.TRUE[i,t-1],,i,t-1])==1)
      CH.TRUE[i,t] <- state
      # Multinomial trials for observation process
      event <- which(rmultinom(1, 1, PSI.OBS[CH.TRUE[i,t],,i,t-1])==1)
      CH[i,t] <- event
    } #t
  } #i
  # Replace the NA and the highest state number (dead) in the file by 0
  CH[is.na(CH)] <- 0
  CH[CH==dim(PSI.STATE)[1]] <- 0
  CH[CH==unobservable] <- 0
  id <- rep(NA, dim(CH)[1])
  z <- apply(CH, 1, function(x) min(which(x != 0)))
  id <- (z != dim(CH)[2])
  return(list(CH=CH[id,], CH.TRUE=CH.TRUE[id,]))
  # CH: capture histories to be used
  # CH.TRUE: capture histories with perfect observation
}

# Execute simulation function
sim <- simul.ms(PSI.STATE, PSI.OBS, marked)
CH <- sim$CH

###function to create capture history character strings (need for input to RMARK)
pasty<-function(x){
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {
    out[i]<-paste(x[i,],collapse="")
  }
  return(out)
}

# Compute vector with occasions of first capture
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

# Recode CH matrix: note, a 0 is not allowed in WinBUGS!
# 1 = seen alive in A, 2 = seen alive in B, 3, seen alive in C, 4 = not seen
rCH <- CH  # Recoded CH
rCH[rCH==0] <- 4

# Function to create known latent states z
known.state.ms <- function(ms, notseen){
  # notseen: label for ‘not seen’
  state <- ms
  state[state==notseen] <- NA
  for (i in 1:dim(ms)[1]){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
}
  return(state)
}

# Function to create initial values for unknown z
ms.init.z <- function(ch, f){
  for (i in 1:dim(ch)[1]){ch[i,1:f[i]] <- NA}
  states <- max(ch, na.rm = TRUE)
  known.states <- 1:(states-1)
  v <- which(ch==states)
  ch[-v] <- NA
  ch[v] <- sample(known.states, length(v), replace = TRUE)
  return(ch)
}

# Bundle data
jags.data <- list(y = rCH, 
                  f = f, n.occasions = dim(rCH)[2], 
                  nind = dim(rCH)[1], 
                  z = known.state.ms(rCH, 4))

# Initial values 
inits <- function(){list(phiA = runif(1, 0, 1), 
                         phiB = runif(1, 0, 1), 
                         phiC = runif(1, 0, 1), 
                         lpsiA = rnorm(2), 
                         lpsiB = rnorm(2), 
                         lpsiC = rnorm(2), 
                         pA = runif(1, 0, 1), 
                         pB = runif(1, 0, 1), 
                         pC = runif(1, 0, 1), 
                         z = ms.init.z(rCH, f))}  


# Parameters monitored
parameters <- c("phiA", "phiB", "phiC", "psiA", "psiB", "psiC", "pA", "pB", "pC")

# MCMC settings
ni <- 10000
nt <- 4
nb <- 5000
nc <- 3

# Call JAGS from R
ms3 <- jagsUI::jags(jags.data, inits, parameters, "ms3-multinomlogit.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE)

print(ms3, digits = 3)
