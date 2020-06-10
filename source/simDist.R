#!/usr/bin/env Rscript
options(warn=1)

# Script to estimate the mean distance between cases
# Given R and D
## LOAD REQUIRED PACKAGES AND FUNCTIONS ##

library(ape)
library(outbreaker)
library(igraph)

# set up list of values to turn into data frame
R0 <- seq(1.3,18,0.1)
D <- rep(NA,length(R0))
prob.dist <- vector("list",length(R0))

# calculate the maximum distance for any simulation
for (i in 1:length(R0)){
  D[i] <- floor(log(1000)/log(R0[i]))
}
maxdist <- max(D)

# calculate values for all R0
for (i in 1:length(R0)){
  
  # initialize a counter
  c <- 0
  
  # get the current R value
  R <- R0[i]
  
  # get the number of generations for this R value
  gen <- D[i]

  # the maximum distance for this particular value of R
  Rdist <- 2*gen
  
  # initialize a list of length maxdist to store counts
  # dists[i] is the count for a distance of i
  dists <- rep(0,2*maxdist)
  
  # initialize value for the total number of comparisons
  N <- 0
  
  repeat {
    
    x <- simOutbreak(R0 = R, infec.curve = c(0,1), n.hosts = 100000, duration = gen, seq.length = 100,
                     rate.import.case = 0, diverg.import = 365)
    
    # make sure we have at least 3 cases in this outbreak
    if (length(x$id)>2){

      # prevent the outbreak from getting too big with large R
      if ((R<1.8) && (length(x$id))>2000){
         next
      }
      
      # turn the outbreak results into an undirected graph
      g <- graph_from_edgelist(cbind(x$ances,x$id)[2:length(x$id),])
      
      # get a matrix of all distances between nodes
      alldists <- distances(g)
      # only look at the bottom half of the matrix (top is redundant)
      # do not keep the diagonal values (0) because we do not care about 0 distance
      alldists[upper.tri(alldists, diag = TRUE)] <- NA
      
      # count how many of each value from 1 to the max dist are in the matrix
      # add this to the running tally for these values of R and D
      for (j in 1:Rdist){
        dists[j] <- dists[j] + sum(alldists==j,na.rm=TRUE)
      }
      
      # calculate the total number of pairs
      cases <- length(x$id)
      N <- N + (cases * (cases-1)) / 2
      
      # increase the counter
      c <- c+1
      
    }
    
    if (c==1000) {
      
      # now calculate the probability of getting any number from 1 to maxdist
      # divide the count for each number by the total number of pairs
      prob.dist[[i]] <- dists/N
      
      break
    }

  }
  
  print(paste(Sys.time()," Calculations finished for R = ",R,sep = ""))
  print(prob.dist[[i]])
  
}

# turn the results into a data frame
res <- data.frame("R"=R0,"gens"=D,do.call(rbind, prob.dist))

# save the data.frame as a csv file
 write.csv(res,"sim_distance_dist.csv",row.names=FALSE)
