nThreads = 4
nblocks = 20*nThreads

t = c(0,0.5,1,1.5,2,2.5,3,4)
S = c(254,235,201,153,121,110,97,83)
I = c(7,14,22,29,20,8,8,0)
N = 254 + 7
R = N - S - I
n = 8

### Prior ####
logprior <- function(param){
  alpha = param[1]
  beta = param[2]
  aprior = dnorm(alpha, mean = 0, sd = 100, log = TRUE)
  bprior = dnorm(beta, mean = 0, sd = 100, log = TRUE)
  #aprior = dlnorm(log_alpha_power,
  #                meanlog = 0,   #alpha_prior_params_lognormal[1]
  #                sdlog = 0.21)+ #alpha_prior_params_lognormal[2]
  #  bprior = dlnorm(log_beta,
  #         meanlog = -5.41,#beta_prior_params_lognormal[1]
  #         sdlog = 0.76) #beta_prior_params_lognormal[2]

  return(aprior+bprior)
}

### Likelihood ####
loglik <- function(par) {

  alpha = exp(par[1])
  beta = exp(par[2])

  fun <- function(k){
    return(log(SIR_prob(t=t[k+1]-t[k],beta=beta,alpha = alpha,
                        S0=S[k],I0=I[k],
                        nSI=S[k]-S[k+1],nIR=R[k+1]-R[k],direction = "Forward",
                        nThreads = nThreads, nblocks = nblocks)[S[k]-S[k+1]+1,R[k+1]-R[k]+1]))
  }
  tmp = sapply(1:(n-1),fun)
  if (is.na(sum(tmp))) {
    print(c(alpha,beta))
  }
  return(sum(tmp))
}

### Derivatives ####
derivatives <- function(par) {

  alpha = exp(par[1])
  beta = exp(par[2])

  fun <- function(k){
    P = SIR_prob(t=t[k+1]-t[k],beta=beta,alpha=alpha,
                 S0=S[k],I0=I[k],
                 nSI=S[k]-S[k+1],nIR=R[k+1]-R[k],direction = "Forward",
                 nThreads = nThreads, nblocks = nblocks)[S[k]-S[k+1]+1,R[k+1]-R[k]+1]

    Pa = SIR_derivatives(t=t[k+1]-t[k],N,alpha,beta,S[k],I[k],
                         A=S[k]-S[k+1],B=R[k+1]-R[k],
                         derivative.order = "alpha", direction = "Forward",
                         nThreads = nThreads, nblocks = nblocks)[S[k]-S[k+1]+1,R[k+1]-R[k]+1]

    Pb = SIR_derivatives(t=t[k+1]-t[k],N,alpha,beta,S[k],I[k],
                         A=S[k]-S[k+1],B=R[k+1]-R[k],
                         derivative.order = "beta", direction = "Forward",
                         nThreads = nThreads, nblocks = nblocks)[S[k]-S[k+1]+1,R[k+1]-R[k]+1]
    return(c(Pa/P,Pb/P))
  }

  #   tmp_a = sapply(1:(n-1),fun_a)
  #   tmp_b = sapply(1:(n-1),fun_b)
  tmp = sapply(1:(n-1), fun)
  return(c(sum(tmp[1,])*alpha, sum(tmp[2,])*beta))
}


########## HMC sampler ####
require(MASS)

computeLoglikelihood <- function(par, truncation = FALSE, gradient = FALSE) {
  if (gradient) {
    return(derivatives(par))
  }
  else {
    return(loglik(par))
  }
}

Potential <- function(par, gradient=FALSE) {
  if (gradient) {
    logPriorGrad <- par/100^2
    logLikelihoodGrad <- derivatives(par)

    return(-(logPriorGrad + logLikelihoodGrad))
  }
  else {
    logPrior <- logprior(par)
    logLikelihood <- loglik(par)

    return(-(logPrior+logLikelihood))
  }
}

hmcsampler <- function(NumOfIterations, Trajectory, NumOfLeapfrog) {
  # Set up the parameters
  BurnIn = floor(0.2*NumOfIterations)
  StepSize = Trajectory/NumOfLeapfrog

  # Allocate output space
  LocationSaved = list()
  Target = vector()

  # Starting point
  locations <- c(log(3.2),log(0.0197))
  N = 1
  P = length(locations)

  Accepted = 0;
  Proposed = 0;
  # Initialize the location
  CurrentLocation = locations;
  CurrentU = Potential(locations)

  # Perform Hamiltonian Monte Carlo
  for (Iteration in 1:NumOfIterations) {
    ProposedLocation = CurrentLocation
    print(c("Iter:",Iteration))
    # Sample the marginal momentum
    CurrentMomentum = mvrnorm(N,rep(0,P),matrix(c(1,0,0,1),2,2))
    ProposedMomentum = CurrentMomentum
    #     print(ProposedMomentum)

    Proposed = Proposed + 1

    # Simulate the Hamiltonian Dynamics
    for (StepNum in 1:NumOfLeapfrog) {
      ProposedMomentum = ProposedMomentum - StepSize/2 * Potential(ProposedLocation,gradient=T)
      ProposedLocation = ProposedLocation + StepSize * ProposedMomentum
      ProposedMomentum = ProposedMomentum - StepSize/2 * Potential(ProposedLocation,gradient=T)
    }

    ProposedMomentum = - ProposedMomentum

    # Compute the Potential
    ProposedU = Potential(ProposedLocation)
    # Compute the Hamiltonian
    CurrentH = CurrentU + sum(CurrentMomentum^2)/2
    ProposedH = ProposedU + sum(ProposedMomentum^2)/2

    # Accept according to ratio
    Ratio = - ProposedH + CurrentH
    if (is.finite(Ratio) & (Ratio > min(0,log(runif(1))))) {
      CurrentLocation = ProposedLocation
      CurrentU = ProposedU
      Accepted = Accepted + 1
    }
    #     print(exp(Ratio))

    # Save if sample is required
    if (Iteration > BurnIn) {
      LocationSaved[[Iteration-BurnIn]] = CurrentLocation
      Target[Iteration-BurnIn] = CurrentU
    }

    # Show acceptance rate every 100 iterations
    if (Iteration%%100 == 0) {
      cat(Iteration, "iterations completed. Acceptance rate: ",Accepted/Proposed,"\n")

      Proposed = 0
      Accepted = 0
    }

    # Start timer after burn-in
    if (Iteration == BurnIn) {
      cat("Burn-in complete, now drawing samples ...\n")
      timer = proc.time()
    }


  }

  time = proc.time() - timer
  acprat = dim(LocationSaved[!duplicated(LocationSaved)])[1]/(NumOfIterations-BurnIn)
  return(list(samples = LocationSaved, target = Target, Time = time, acprat = acprat))

}


#### Run HMC ####
system.time(hmc <-hmcsampler(NumOfIterations = 10000, Trajectory = 0.2, NumOfLeapfrog = 2))
mcmc <- matrix(unlist(hmc$samples),byrow=T,ncol=2)
colnames(mcmc) = c("logalpha", "logbeta")

plot(mcmc[,1], type="l")
plot(mcmc[,2], type="l")
require(coda)
effectiveSize(mcmc)

#### Density plot ####
require(ggplot2)
x = as.vector(mcmc[,1])
y = as.vector(mcmc[,2])
df <- data.frame(x, y)

pdf("codes/HMC/HMC_density.pdf")
ggplot(df,aes(x=x,y=y)) +
  stat_density2d(aes(fill=..level..), geom="polygon", h = 0.2) +
  #   geom_density(adjust=5) +
  #   stat_density2d(aes(fill=..level..), geom="density2d")
  scale_fill_gradient(low="grey85", high="grey35", guide = FALSE) +
  xlab(expression(paste("log (",gamma,")",sep=""))) +
  ylab(expression(paste("log (",beta,")",sep=""))) +
  #   geom_point(aes(x = log(2.73), y = log(0.0178)), size = 5, shape = 3) +
  #   geom_point(aes(x = log(3.39), y = log(0.0212)), size = 5, shape = 4) +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22,face="bold"))
dev.off()

mean(exp(mcmc[,1]))
quantile(exp(mcmc[,1]), probs = c(0.025,0.975))
mean(exp(mcmc[,2]))
quantile(exp(mcmc[,2]), probs = c(0.025,0.975))

