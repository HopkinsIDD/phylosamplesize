#!/usr/bin/env Rscript
options(warn=1)

### Script to simulate outbreaks and calculate necessary parameters

## LOAD REQUIRED PACKAGES AND FUNCTIONS ##

library(ape)
library(outbreaker)
library(optparse)
library(igraph)
source("source/sim_func.R")


## SET UP COMMAND LINE ARGUMENTS ##

option_list <- list( 
  make_option(c("-m","--mu"), default=0.0000001, metavar="mu",
              help="Mutations per site per generation; can be a number or (min,max) vector"),
  make_option(c("-R","--R"), default=2, metavar="R",
              help = "Outbreak reproductive number; can be a number or (min,max) vector"),
  make_option(c("-N","--numsim"), default=100, metavar="number of simulations", 
              help="Number of simulations to run"),
  make_option(c("-p","--rho"), default=1, metavar="sampling proportion",
              help="Proportion of cases to sample from an outbreak"),
  make_option(c("-M","--samplesize"), default=-1, metavar="cases to sample",
              help="Number of cases to sample from an outbreak"),
  make_option(c("-S","--S"), default=100000, metavar="number of susceptibles",
              help="Initial number of susceptibles in the population"),
  make_option(c("-G","--gendist"), default=c(0,1), metavar="generation time distribution",
              help="Generation time distribution by timestep"),
  make_option(c("-L","--length"), default=1000, metavar="genome length",
              help="Length of pathogen genome"),
  make_option("--rimport", default=0, metavar="importation rate",
              help="Rate at which cases are imported at each time step"),
  make_option("--dimport", default=365, metavar="importation divergence",
              help="Number of time steps to the MRCA of all imported cases"),
  make_option(c("-C","--minsize"), default=100, metavar="minimum outbreak size",
              help="Number of cases required to constitute an outbreak"),
  make_option(c("-X","--maxsize"), default=2000, metavar="maximum outbreak size",
              help="Maximum number of cases allowed in an accepted outbreak"),
  make_option(c("-F","--finalsize"), default=1000, metavar="final outbreak size",
              help="Desired number of infected individuals at end of outbreak"),
  make_option(c("-o","--outdir"), help="File path to csv file for saving results")
)

## SIMULATE OUTBREAKS ##

# get command line options
# if help option encountered print help and exit
# if options not found on command line then use default
opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

# deal with R and mu vectors if applicable
if(is.character(opt$mu)) { opt$mu <- as.numeric(unlist(strsplit(opt$mu,split = ","))) }
if(is.character(opt$R)) { opt$R <- as.numeric(unlist(strsplit(opt$R,split = ","))) }

## SIMULATE OUTBREAKS FOR VARIABLE MU AND R ##
  
stopifnot(length(opt$R)==2)
  
steps <- ceiling(sqrt(opt$numsim)) # ensures that we have at least N simulations
mus <- seq(opt$mu[1],opt$mu[2],length.out = steps)
Rs <- seq(opt$R[1],opt$R[2],length.out = steps)
  
# get all combinations of mu and R to use for simulations
params <- expand.grid(mu = mus, R = Rs)
params <- cbind("i"=seq(steps^2),params)
print(params)
  
# get the maximum number of generations possible given these R values
maxgens <- floor( log(1000)/log(min(Rs)) )
#maxgens <- floor( log(opt$finalsize)/log(min(Rs)) )
maxgens <- 2*maxgens

# store all data in a list before combining
simdata_list = list()
gensdata_list = list()

# simulate outbreaks for all combinations in params
for(i in 1:nrow(params)) {
  p <- params[i,]

  # RUN SIMULATION
  
  i <- as.numeric(p["i"])
  mu <- as.numeric(p["mu"])
  R <- as.numeric(p["R"])
  D <- floor( log(opt$finalsize) / log(R) )
  x <- true.outbreak(mu=mu,R=R,D=D,S=opt$S,gen=opt$gendist,L=opt$length,
	       rimport=opt$rimport,dimport=opt$dimport,C=opt$minsize,X=opt$maxsize)
  print(paste("SIMULATION",i,sep=" "))
  print(paste("true outbreak found, simulation size:",x$n,sep=" "))

  # CALCULATE PARAMETERS
  
  data <- cbind("sim"=i,calc.params(simulation=x,M=opt$samplesize,rho=opt$rho),row.names=NULL)
  simdata_list[[i]] <- data
  
  # DETERMINE THE GENERATION DISTRIBUTION
  
  # turn the outbreak results into an undirected graph
  g <- graph_from_edgelist(cbind(x$ances,x$id)[2:length(x$id),])
  
  # get a matrix of all distances between nodes
  alldists <- distances(g)
  
  # only look at the bottom half of the matrix (top is redundant)
  # do not keep the diagonal values (0) because we do not care about 0 distance
  alldists[upper.tri(alldists, diag = TRUE)] <- NA
  
  # count how many of each value from 1 to the max dist (52) are in the matrix
  dists <- rep(0,maxgens)
  for (j in 1:maxgens){ dists[j] <- dists[j] + sum(alldists==j,na.rm=TRUE) }
  
  dists <- dists/sum(dists)
  dists <- cbind("sim"=i,"R"=R,data.frame(matrix(dists,nrow=1)),row.names=NULL)
  gensdata_list[[i]] <- dists

  }

simdata = do.call(rbind, simdata_list)
gensdata = do.call(rbind, gensdata_list)
write.csv(simdata,paste(opt$outdir,"simdata_var_gen_N",opt$numsim,"_rho",100*opt$rho,".csv",sep=""),row.names=FALSE)
write.csv(gensdata,paste(opt$outdir,"gendata_var_sim_N",opt$numsim,"_rho",100*opt$rho,".csv",sep=""),row.names=FALSE)
