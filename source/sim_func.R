#' Simulate an outbreak with a minimum number of cases
#'
#' Uses the simOutbreak function from the outbreaker package
#' to run a simulation given a set of parameters
#' true.outbreak will continue running simulations until a simulation
#' has the required minimum number of cases, defined by the parameter C
#' 
#' @param mu pathogen mutation rate in mutations/site/generation
#'           assumes that transition and transversion rates are equal (mu.transi=mu.transv=mu/2)
#' @param R pathogen reproductive number (R0)
#' @param D number of generations to run simulation for (duration)
#' @param S number of susceptible individuals in starting population (n.hosts)
#' @param gen distribution of individual infectiousness over time (infec.curve)
#' @param L length of sequences to simulation (seq.length)
#' @param rimport rate at which cases are imported at each time step (rate.import.case)
#' @param dimport the number of time steps to the MRCA of all imports (diverg.import)
#' @param C the minimum number of cases required to accept the simulated outbreak
#' @param X the maximum number of cases allowed to accept the simulated outbreak
#' 
#' @return a simOutbreak-like list describing a simulated outbreak with at least C cases
#' 
true.outbreak <- function(mu,R,D,S,gen,L,rimport,dimport,C,X) {
  
  repeat {
    
    # run a basic simOutbreak simulation
    x <- simOutbreak(R0 = R, infec.curve = gen, n.hosts = S, duration = D,
                     seq.length = L, mu.transi = mu/2, mu.transv = mu/2,
                     rate.import.case = rimport, diverg.import = dimport)
    
    # now test if we have achieved a 'true outbreak'
    if ((x$n >= C) & (x$n <= X)) {
      
      # store mu and R for future reference
      x$mu <- mu
      x$R <- R
      
      # return x and exit out of the loop
      # return as a list that roughly mirrors the simOutbreak ogject format
      x <- list("n"=x$n, "rho"=1, "mu"=x$mu, "R"=x$R, "dna"=x$dna,"id"=x$id,"ances"=x$ances)
      
      return(x)
    }
    
  }
  
}


#' Function to calculate a genetic distance matrix
#'
#' Calculates a distance matrix using the dist.dna function from the ape package
#' from simulated sequences stored in a simOutbreak or simOutbreak-like object
#' distances are in number of mutations between two sequences (of fixed length)
#' and no corrections are provided for missing data
#' 
#' @param simulation a simOutbreak object or simOutbreak-like list
#' 
#' @return a genetic distance matrix
#' 
get.dm <- function(simulation) {
  
  # first assign rownames to dnabin object to ensure dm is labeled correctly
  rownames(simulation$dna) <- simulation$id
  
  # calculate the distance matrix in number of mutations
  dm <- dist.dna(simulation$dna, model = "N", as.matrix = TRUE)
  
  return(dm)
  
}


#' Calculate a transmission link matrix
#'
#' Calculates a matrix representing true transmission links between cases
#' 
#' @param simulation a simOutbreak object or simOutbreak-like list
#' 
#' @return a symmetric matrix where '1' indicates an infector-infectee pair
#'         directionality of links (i.e. who infected whom) is not preserved
#' 
get.pm <- function(simulation) {
  
  # initialize a matrix of zeroes
  pm <- matrix(0,length(simulation$id),length(simulation$id))
  
  # set row and column names
  # if the simulation has been sampled, these may not be sequential
  rownames(pm) <- simulation$id
  colnames(pm) <- simulation$id
  
  # loop through cases and fill in matrix using ancestors 
  for (i in seq(length(simulation$id))){
    
    # ignore any NA values
    # should only occur if the first case is samples
    if (is.na(simulation$ances[i])){
      next
    }
    
    # fill in matrix based on matching infector/infectee pairs
    # only fill in if ancestor is in sample
    infector <- as.character(simulation$ances[i])
    infectee <- as.character(simulation$id[i])
    
    if (is.element(infector, simulation$id)){
      pm[infectee,infector] <- 1
      pm[infector,infectee] <- 1
    }
  }
  
  return(pm)
  
}


#' Subset a simulated outbreak given a sampling proportion
#'
#' Takes an outbreak simulated by true.outbreak (or directly by simOutbreak)
#' and randomly samples a number of cases defined by the sampling proportion
#' 
#' @param simulation a simOutbreak object or simOutbreak-like list
#' @param rho the desired proportion of cases to sample
#' 
#' @return a simOutbreak-like list containing only sampled cases
#'         some lists included in the original simOutbreak object are dropped:
#'         $onset, $dynam, $group, $nmut, $ngen, and $call are dropped for simplicity
#' 
#' @seealso subset.sim.M
#' 
subset.sim.rho <- function(simulation, rho){
  
  # use rho to determine the number of cases to sample
  n <- round(rho * simulation$n)
  
  # select the sampled cases from the list of ids
  ids <- sort(sample(simulation$id,n))
  
  # set up the list corresponding to the subsampled simulation
  # ensure that $dna and $ances include only sampled cases
  newsim <- list("n"=n, "rho"=rho, "mu"=simulation$mu, "R"=simulation$R, "dna"=simulation$dna[ids,],
                 "id"=ids, "ances"=simulation$ances[ids])
  
  return(newsim)
  
}


#' Subset a simulated outbreak given a sample size
#'
#' Takes an outbreak simulated by true.outbreak (or directly by simOutbreak)
#' and randomly samples a number of cases
#' 
#' @param simulation a simOutbreak object or simOutbreak-like list
#' @param M the desired sample size
#' 
#' @return a simOutbreak-like list containing only sampled cases
#'         some lists included in the original simOutbreak object are dropped:
#'         $onset, $dynam, $group, $nmut, $ngen, and $call are dropped for simplicity
#' 
#' @seealso subset.sim.rho
#' 
subset.sim.M <- function(simulation, M){
  
  # calculate the sampling proportion
  rho <- M / simulation$n
  
  # select the sampled cases from the list of ids
  ids <- sort(sample(simulation$id,M))
  
  # set up the list corresponding to the subsampled simulation
  # ensure that $dna and $ances include only sampled cases
  newsim <- list("n"=M, "rho"=rho, "mu"=simulation$mu, "R"=simulation$R, "dna"=simulation$dna[ids,],
                 "id"=ids, "ances"=simulation$ances[ids])
  
  return(newsim)
  
}


#' Calculate linkage metrics at a particular mutation threshold
#'
#' Calculates the sensitivity and specificity of a particular
#' genetic distance threshold (given in number of mutations)
#' 
#' @param dm a genetic distance matrix
#' @param pm a matrix indicating true transmission pairs
#' @param thresh numeric threshold at which to test linkage
#' 
#' @return a list of values related to sensititivy and specificity:
#'         the threshold used, sensitivity (eta), specificity (chi),
#'         the false discovery rate (fdr),
#'         the number of true positives and negatives,
#'         and the number of false positivies and negatives
#' 
get.metrics <- function(dm,pm,thresh) {
  
  # determine which cases pairs diverge by fewer mutations than the threshold
  # multiply the resulting true/false matrix by 2 to get 2/0 values
  link.dm <- (dm < thresh) * 2
  
  # add the resulting matrix to the links matrix
  # to get unique values corresponding to tp, tn, fp, fn
  link.dm <- link.dm + pm
  
  # turn unneeded values to -1 so they are not counted
  link.dm[upper.tri(link.dm)] <- -1
  diag(link.dm) <- -1
  
  # get desired values
  tp <- sum(link.dm==3)
  fp <- sum(link.dm==2)
  fn <- sum(link.dm==1)
  tn <- sum(link.dm==0)
  
  eta <- tp / (tp + fn)
  chi <- tn / (tn + fp)
  fdr <- fp / (fp + tp)
  
  values <- c(thresh,eta,chi,fdr,tp,tn,fp,fn)
  
  return(values)
  
}


#' Calculate linkage matrix for a simulation at all possible thresholds
#'
#' For a simulation, calculates the specificity, sensitivity, and related metrics
#' using all possible genetic distance thresholds to determine linkage
#' if subsampling is desired, calculations are performed on the subsampled simulation
#' 
#' @param simulation a simOutbreak object or simOutbreak-like list
#' @param M the desired sample size
#'        defaults to -1 if no subsampling based on sample size is desired
#' @param rho the desired sampling proportion
#'        defaults to 1 if no subsampling based on proportion is desired
#' 
#' @return a matrix containing all relevant inputs (mu, R, etc.)
#'         and all calculated metrics assessing the ability of
#'         genetic distance to predict transmission linkage
#' 
calc.params <- function(simulation,M=-1,rho=1){
  
  # get distance and links matrix for full simulation
  # if the full simulation has many cases, subset first to save time
  if (simulation$n > 2000){
    
    simprox <- subset.sim.M(simulation,2000)
    dm <- get.dm(simprox)
    pm <- get.pm(simprox)
    
  } else {
    
    dm <- get.dm(simulation)
    pm <- get.pm(simulation)
    
  }
  
  # calculate the maximum value of the distance matrix
  dm.max <- max(dm)
  
  # calculate the empirical reproductive number for this simulation
  # because the final generation does not infect anyone,
  # this will not match the value used to run the simulation
  Re <- mean(rowSums(pm)) - 1
  
  # get distances of linked and unlinked pairs
  dm.data <- lapply(seq(from = 0, to = dm.max+1), function(x) {
    
    # change diagonal and upper triangle to -1 to avoid repeats
    dm[upper.tri(dm)] <- -1
    diag(dm) <- -1
    
    if (x==0){
      
      dm.mod <- ifelse(dm==0,-2,dm)
      d.link <- sum((dm.mod*pm)==-2)
      d.nolink <- sum((dm.mod*(1-pm))==-2)
      
    } else {
      
      d.link <- sum((dm*pm)==x)
      d.nolink <- sum((dm*(1-pm))==x)
      
    }
    
    res <- c(d.link,d.nolink)
    
    return(res)
    
  })
  dm.data <- matrix(unlist(dm.data),ncol=2,byrow=TRUE)
  
  # get metrics for every possible threshold value
  full.thresh.data <- lapply(seq(from = 0, to = dm.max+1), function(iter) {get.metrics(dm,pm,thresh=iter)})
  full.thresh.data <- matrix(unlist(full.thresh.data),ncol=8,byrow=TRUE)
  
  # subset the simulation
  if (M<0) { subsim <- subset.sim.rho(simulation,rho)
  } else { subsim <- subset.sim.M(simulation,M)}
  
  # now calculate distance and links matrix for subsetted simulation
  dm.sub <- get.dm(subsim)
  pm.sub <- get.pm(subsim)
  
  # get metrics for all threshold values
  # we will use all thresholds from full data sets (may be more than in subset)
  sub.thresh.data <- lapply(seq(from = 0, to = dm.max+1), function(iter) {get.metrics(dm.sub,pm.sub,thresh=iter)})
  sub.thresh.data <- matrix(unlist(sub.thresh.data),ncol=8,byrow=TRUE)
  sub.thresh.data <- sub.thresh.data[,2:8] # do not repeat threshold value
  
  # return the full matrix of values we want
  thresh.data <- cbind(rep(simulation$n,dm.max+2),rep(subsim$n,dm.max+2),rep(subsim$rho,dm.max+2),
                       rep(subsim$mu,dm.max+2),rep(subsim$R,dm.max+2),
                       rep(Re,dm.max+2),full.thresh.data,sub.thresh.data,dm.data)
  colnames(thresh.data) <- 
    c("N","M","rho","mu","R","Re","t","eta.full","chi.full","fdr.full",
      "tp.full","tn.full","fp.full","fn.full","eta.sub","chi.sub","fdr.sub","tp.sub","tn.sub",
      "fp.sub","fn.sub","d.link","d.nolink")
  
  return(thresh.data)
  
}