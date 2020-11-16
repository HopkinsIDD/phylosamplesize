##' Function to plot the relationship between sample size and sampling proportion
##' versus the false discovery rate or expected number of identified links
##'
##' @param chi number giving the specificity of the linkage criteria
##' @param eta number giving the sensititivy of the linkage criteria
##' @param R number giving the effective reproductive number
##' @param rho.disc vector giving discrete values of the sampling proportion to evaluate
##' @param rho.cont vector giving continuous values of the sampling proportion to evaluate
##' @param M.disc vector giving discrete values of the sample size to evaluate
##' @param M.cont vector giving continuous values of the sample size to evaluate
##' @param assumption string indicating which assumptions about transmission and linkage criteria to use

plt.fdr.pairs <- function(chi,eta,R=NULL,rho.disc,rho.cont,M.disc,M.cont,assumption){
  
  if (!assumption %in% c("stsl","mtsl","mtml")) { stop("Incorrect assumption argument") }
  obs_pairs_func = paste("obs_pairs_",assumption,sep="")
  prob_trans_func = paste("prob_trans_",assumption,sep="")
  
  # set up the dataframe to be used in plotting vs rho
  g1 <- g3 <- expand.grid(rho.cont)
  names(g1) <- names(g3) <- c('x')
  
  for (i in seq(1,length(M.disc))){
    cname <- M.disc[i] # set name for column to be added
    
    # expected number of links
    g1 <- cbind(g1, doCall(get(obs_pairs_func),chi = chi, eta = eta, rho = g1$x, M = M.disc[i], R = R))
    colnames(g1)[length(colnames(g1))] <- cname
    g1_melted <- melt(g1,id='x')
    
    # false discovery rate
    g3 <- cbind(g3, 1-doCall(get(prob_trans_func),chi = chi, eta = eta, rho = g3$x, M = M.disc[i], R = R))
    colnames(g3)[length(colnames(g3))] <- cname
    g3_melted <- melt(g3,id='x')
    
  }
  
  # set up the dataframe to be used in plotting vs M
  g2 <- g4 <- expand.grid(M.cont)
  names(g2) <- names(g4) <- c('x')
  
  for (i in seq(1,length(rho.disc))){
    cname <- rho.disc[i] # set name for column to be added
    
    # expected number of links
    g2 <- cbind(g2, doCall(get(obs_pairs_func),chi = chi, eta = eta, rho = rho.disc[i], M = g2$x, R = R))
    colnames(g2)[length(colnames(g2))] <- cname
    g2_melted <- melt(g2,id='x')
    
    # false discovery rate
    g4 <- cbind(g4, 1-doCall(get(prob_trans_func),chi = chi, eta = eta, rho = rho.disc[i], M = g4$x, R = R))
    colnames(g4)[length(colnames(g4))] <- cname
    g4_melted <- melt(g4,id='x')
  }
  
  # set up colors for plotting
  orange = brewer.pal(n = 9, "Oranges")[5:9]
  blue = brewer.pal(n = 9, "Blues")[5:9]
  
  p1 <- ggplot(g1_melted, aes(x = x, y = value, colour = variable)) +
    geom_line(show.legend = FALSE, size = 1) +
    xlab("proportion sampled") +
    ylab("expected # of links") +
    scale_y_continuous(breaks = c(0,100,200,300,400,500), limits = c(0,500)) +
    theme_classic() +
    theme(axis.title=element_text(size=10,face="bold")) +
    scale_colour_manual(values=orange)
  
  p2 <- ggplot(g2_melted, aes(x = x, y = value, colour = variable)) +
    geom_line(show.legend = FALSE, size = 1) +
    xlab("sample size") +
    ylab("expected # of links") +
    scale_y_continuous(breaks = c(0,100,200,300,400,500), limits = c(0,500)) +
    theme_classic() +
    theme(axis.title=element_text(size=10,face="bold")) +
    scale_colour_manual(values=blue)
  
  p3 <- ggplot(g3_melted, aes(x = x, y = value, colour = variable)) +
    geom_line(size = 1) +
    xlab("proportion sampled") +
    ylab("false discovery rate") +
    theme_classic() +
    theme(legend.direction = "horizontal",
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.8,"line"),
          legend.text = element_text(size = 9),
          legend.justification = 'right',
          axis.title=element_text(size=10,face="bold")) +
    scale_colour_manual(values = orange, name = "sample size")

  p4 <- ggplot(g4_melted, aes(x = x, y = value, colour = variable)) +
    geom_line(size = 1) +
    xlab("sample size") +
    ylab("false discovery rate") +
    ylim(0,1) +
    theme_classic() +
    theme(legend.direction = "horizontal",
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.8,"line"),
          legend.text = element_text(size = 9),
          legend.justification = 'right',
          axis.title=element_text(size=10,face="bold")) +
    scale_colour_manual(values = blue, name = "% sampled")
  
  # arrange the plots
  
  # extract legends from plots
  leg3 <- cowplot::get_legend(p3)
  leg4 <- cowplot::get_legend(p4)
  
  # save a title for the plot
  if (assumption=="stsl") { title <- paste("specificity = ",chi,", sensitivity = ",eta) }
  else { title <- paste("specificity = ",chi,", sensitivity = ",eta,", R = ",R) }
  
  # arrange the plot
  blankPlot <- ggplot()+geom_blank(aes(1,1)) + cowplot::theme_nothing()
  grid.arrange(p1,blankPlot,p2,
               p3+theme(legend.position = "none"),blankPlot,p4+theme(legend.position = "none"),
               leg3,blankPlot,leg4,
               nrow=3,heights=c(4,4,1),widths=c(9,1,9),
               vp=viewport(width = 0.95, height = 0.95),
               top=textGrob(title,vjust=-1,gp=gpar(fontface="bold",fontsize=11)))
  
}


##' Function to plot a heatmap of false discovery rate values
##' at varying specificity and sensitivity values and a fixed sample number/proportion
##'
##' @param chi vector giving the range of the specificity of the linkage criteria
##' @param eta vector giving the range of the sensititivy of the linkage criteria
##' @param rho number giving the sampling proportion
##' @param M number giving the sample size
##' @param R number giving the effective reproductive number
##' @param assumption string indicating which assumptions about transmission and linkage criteria to use

plt.heatmap.fdr <- function(chi,eta,rho,M,R=NULL,assumption){
  
  if (!assumption %in% c("stsl","mtsl","mtml")) { stop("Incorrect assumption argument") }
  prob_trans_func = paste("prob_trans_",assumption,sep="")
  
  # create a modified chi to store 1-specificity
  chi.mod <- 1-chi
  
  # create matrix with modified chi and eta values
  g <- expand.grid(chi.mod,eta)
  names(g) <- c('chi.mod','eta')
  
  # use probability equation to add false discovery probability
  g <- cbind(g, 1-doCall(get(prob_trans_func),chi = 1-g$chi.mod, eta = g$eta, rho = rho, M = M, R = R))
  colnames(g)[length(colnames(g))] <- "FDR"
  
  # save the plot title
  if (assumption=="stsl") { plt.title <- paste("proportion sampled = ",rho,", sample size = ", M) }
  else { plt.title <- paste("proportion sampled = ",rho,", sample size = ", M,", R = ", R) }
  
  # make the plot
  
  p <- ggplot(g, aes(x=chi.mod, y=eta, z=FDR))
  
  p + stat_contour()
  
  p + geom_raster(aes(fill = FDR), hjust=0, vjust=0.2, interpolate=FALSE) + 
    scale_x_log10(limits = c(0.00001,0.1), expand = c(0.005,0)) +
    scale_y_continuous(limits = c(0,1.01), expand = c(0.005, 0)) +
    scale_fill_viridis(begin = 0.1, end = 0.8, direction = -1, limits = c(0,1),
                       breaks = c(0.0,0.25,0.5,0.75,1.0), labels = c("0.0",0.25,0.5,0.75,"1.0")) +
    geom_contour(aes(z = FDR), color='black', size=0.2, alpha = 0.5) +
    geom_text_contour(aes(z = FDR), nudge_y=0.0175) +
    theme_bw() +
    ylab('sensitivity of the linkage criteria') + 
    xlab('1 - specificity of the linkage criteria') + 
    ggtitle(plt.title) +
    theme(axis.text.x=element_text(size=10,vjust=-0.6),
          axis.text.y=element_text(size=10),
          axis.title.x=element_text(size=10,face="bold",vjust=-2),
          axis.title.y=element_text(size=10,face="bold",vjust=4),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0.2,0.2,0.2,0.2),"in"),
          plot.title = element_text(size=11,face="bold",hjust=0.5))
}


##' Function to get plot distance distribution from a simulated outbreak
##' and corresponding ROC curve
##'
##' @param outbreak path to simulated outbreak SimOutbreak object

plt.gendist.sim <- function(outbreak){
  
  # load a saved outbreak
  load(outbreak)
  
  # subsample this outbreak
  # choose cases 10-60
  # these are early in the outbreak
  ids <- seq(11,60)
  subsim <- list("n"=50, "dna"=sim_data$dna[ids,],"id"=ids, "ances"=sim_data$ances[ids])
  
  # calculate a distance matrix from outbreak data
  sim_dm <- get.dm(subsim)
  
  # get the transmission link matrix from the outbreak data
  sim_pm <- get.pm(subsim)
  
  # get the genetic distance and list status for all pairs
  gendist <- sim_dm[upper.tri(sim_dm)]
  linkstat <- sim_pm[upper.tri(sim_pm)]
  link_dat <- data.frame("gendist"=gendist,"link"=linkstat)
  
  # calculate ROC curve
  max_thresh <- max(sim_dm)+1
  sens <- rep(NA,max_thresh)
  spec <- rep(NA,max_thresh)
  for (i in 0:max_thresh){
    val <- get.metrics(sim_dm,sim_pm,i)
    sens[i+1] <- val[5]/(val[5]+val[8])
    spec[i+1] <- val[6]/(val[6]+val[7])
  }
  
  roc_data <- data.frame("thresh"=0:max_thresh,"sensitivity"=sens,"specificity"=1-spec)
  
  # get the optimal value for the ROC plot
  # use point closest to corner
  
  roc_optim <- roc_data
  roc_optim$dist <- sqrt((1-roc_optim$sensitivity)^2 + (roc_optim$specificity)^2)
  roc_optim <- roc_optim[2:(dim(roc_optim)[1]),] # remove first row with zero threshold
  optim_value <- min(roc_optim$dist)
  threshold <- roc_optim[roc_optim$dist==optim_value,]$thresh
  optim_point <- c(roc_optim[roc_optim$thresh==threshold,]$specificity,
                   roc_optim[roc_optim$thresh==threshold,]$sensitivity)
  
  # plot the distributions of genetic distance
  
  pal <- brewer.pal(5, "PuOr")
  
  p1 <- ggplot(link_dat,aes(x=gendist,color="white",fill=as.character(link))) +
    geom_histogram(aes(y=..density..),alpha=0.5,binwidth=1,boundary=-1.5,position="identity") +
    scale_fill_manual(name="", values=c(pal[2], pal[5])) +
    scale_color_manual(name="", values=c(alpha("white",0.5),alpha("white",0.5))) +
    scale_x_continuous(expand = c(0,0),limits = c(-1,35)) +
    scale_y_continuous(expand = c(0,0),limits = c(0,0.4)) +
    xlab('Genetic distance') + ylab('Density') +
    theme_classic() + theme(legend.position='none') +
    geom_vline(xintercept=threshold, linetype=2, size=0.5, color = "black")
  
  # plot the ROC curve
  
  r1 <- ggplot(roc_data, aes(x=specificity, y=sensitivity)) +
    geom_line(size=1.5) + xlim(0,1) + ylim(0,1) +
    geom_point(x=optim_point[1], y=optim_point[2],
               size = 3, stroke=1, shape=21, fill = "chartreuse") +
    geom_abline(slope=1, intercept = 0, linetype=2, alpha=0.5) +
    xlab("1-specificity") +
    theme_bw()
  
  ggdraw(p1) + cowplot::draw_plot(r1, hjust=-0.16, vjust=-0.16, scale=0.6)
  
}


##' Function to plot theoretical genetic distance distribution given a mutation rate
##' and corresponding ROC curve
##'
##' @param mut_rate mean number of mutations per generation (assumed to be poisson distributed)
##' @param mean_gens_pdf the density distribution of the mean number of generations between cases
##'       (assumes the index of the vector is the discrete distance between cases)
##' @param max_link_gens the maximium generations of separation for linked pairs
##' @param max_gens the maximum number of generations to consider, defaults to the highest
##'        number of generations in mean_gens_pdf with a non-zero probability
##' @param max_dist the maximum distance to calculate, defaults to max_gens * 99.9th percentile
##'       of mut_rate poisson distribution
##'

plt.gendist.mut <- function(mut_rate, mean_gens_pdf, max_link_gens=1,
                            max_gens=which(mean_gens_pdf!=0)[length(which(mean_gens_pdf!=0))],
                            max_dist = max_gens*qpois(.999,mut_rate)){
  
  # get theoretical distributions based on mutation rate and generation parameters
  dists <- as.data.frame(
    phylosamp::gen_dists(mut_rate = mut_rate, mean_gens_pdf = mean_gens_pdf,
                         max_link_gens = max_link_gens, max_gens = max_gens, max_dist = max_dist))
  # reshape dataframe for plotting
  dists <- melt(dists, id.vars = "dist", variable.name = "status", value.name = "prob")
  
  # set up the ROC curve based on these same parameters
  roc_calc <- phylosamp::sens_spec_roc(cutoff = 1:(max(dists$dist)-1), mut_rate = mut_rate,
                                       mean_gens_pdf = mean_gens_pdf,max_link_gens = max_link_gens,
                                       max_gens = max_gens, max_dist = max_dist)
  
  # get the optimal value for the ROC plot
  # use point closest to corner
  
  roc_optim <- roc_calc
  roc_optim$dist <- sqrt((1-roc_optim$sensitivity)^2 + (roc_optim$specificity)^2)
  roc_optim <- roc_optim[2:(dim(roc_optim)[1]),] # remove first row with zero threshold
  optim_value <- min(roc_optim$dist)
  threshold <- roc_optim[roc_optim$dist==optim_value,]$cutoff
  optim_point <- c(roc_optim[roc_optim$cutoff==threshold,]$specificity,
                   roc_optim[roc_optim$cutoff==threshold,]$sensitivity)
  
  # plot the distributions of genetic distance
  
  pal <- brewer.pal(5, "PuOr")
  
  p1 <- ggplot(dists, aes(x=dist, y=prob, fill=status)) +
    geom_bar(alpha=0.5, stat="identity", position="identity") +
    scale_fill_manual(name="", values=c(pal[5], pal[2])) +
    scale_x_continuous(limits = c(-1,35), expand = c(0,0)) +
    scale_y_continuous(limits = c(0,0.4), expand = c(0,0)) +
    xlab('Genetic distance') + ylab('Density') +
    theme_classic() + theme(legend.position='none') +
    geom_vline(xintercept=threshold, linetype=2, size=0.5, color = "black")
  
  # add ROC curve on top
  
  r1 <- ggplot(data=roc_calc, aes(x=specificity, y=sensitivity)) +
    geom_line(size=1.5) + xlim(0,1) + ylim(0,1) +
    geom_point(x=optim_point[1], y=optim_point[2],
               size = 3, stroke=1, shape=21, fill = "chartreuse") +
    geom_abline(slope=1, intercept = 0, linetype=2, alpha=0.5) +
    xlab("1-specificity") +
    theme_bw()
  
  ggdraw(p1) + cowplot::draw_plot(r1, hjust=-0.16, vjust=-0.16, scale=0.6)
  
}


##' Function to calculate the theoretical false discovery rate given simulation data
##'
##' @param simdata generalized path to simulation data (except rho value)
##' @param rho_values vector of rho values  (length=4)
##' @param max_sim_size maximum simulation size to allow (to prevent over-sparse outbreaks)
##' @param sens_spec_method method to calculate sensitivity and specificity: "sim" or "mutrate"
##' @param mgd simulated distribution of generations for different R values
##' @param mgd_type indicate if R value or simulation number should be used to select the distribution
##' @param outdir because calculations can take a while with a large number of simulations,
##            save the final full data frame to a prespecified path

calc.tfdr <- function(simdata,rho_values,max_sim_size,sens_spec_method,mgd,mgd_type="R",outdir){
  
  if (length(rho_values)!=4){ stop("length of rho must equal 4") }
  if (mgd_type=="sim"){ stopifnot(length(mgd)==4) } # ensure mgd_type matches the data provided
  if (mgd_type=="R"){ stopifnot(length(mgd)>4) } # ensure mgd_type matches the data provided

  # create empty list in which to store data
  data_list <- vector(mode = "list", length = 4)
    
  # load simulation data
  for (i in 1:4){
    simpath <- paste(simdata,"_rho",as.character(100*rho_values[i]),".csv",sep="")
    data_list[[i]] <- read.csv(simpath)
    }
    
  # loop through all available data (different rho values) and calculate tfdr
  for (i in 1:length(data_list)){
      
    tmp <- data_list[[i]]
      
    # filter out simulations with more than 2000 cases
    # this makes true positives too sparse for calculations
    # shouldn't do anything if max_sim_size used in true.outbreak
    tmp <- filter(tmp,N<=max_sim_size)
    
    # get the maximum threshold value for each simulation
    # do this first in case we subset the data later
    tmp <- tmp %>% group_by(sim) %>% mutate(max.t = max(t)) %>% ungroup()
    tmp <- data.frame(tmp)
    
    # ignore rows with FDR = NA
    tmp <- filter(tmp,!is.na(fdr.sub))
    
    # FOR TESTING ONLY: randomly sample rows of the dataframe
    #tmp <- tmp[sample(nrow(tmp), 1000), ]
    
    # determine method for calculating specificity and sensitivity
    if (sens_spec_method == "sim"){
      
      # calculate theoretical fdr using sensitivity and specificity from simulations
      tmp <- mutate(tmp,tfdr = 1-phylosamp::prob_trans_mtml(
        chi = chi.full, eta = eta.full, rho = rho, M = M, R = Re))
      
    } else if (sens_spec_method == "mutrate"){
        
      # calculate theoretical fdr using sensitivity and specificity based on mutation rate
      if (mgd_type=="R"){ newcols <- t(pbapply(tmp,1,function(x)
        sens_spec_calc(cutoff = x["t"],mut_rate = 1000*x["mu"],
                       mean_gens_pdf = mgd[mgd["R"]==round(x["R"],1)][3:dim(mgd)[2]],
                       max_link_gens = 1,max_dist = x["max.t"]+1)))[,2:3]
      } else if (mgd_type=="sim") { 
        mgd_rho <- mgd[[i]]
        newcols <- t(pbapply(tmp,1,function(x)
        sens_spec_calc(cutoff = x["t"],mut_rate = 1000*x["mu"],
                       mean_gens_pdf = mgd_rho[mgd_rho["sim"]==x["sim"]][3:dim(mgd_rho)[2]],
                       max_link_gens = 1,max_dist = x["max.t"]+1)))[,2:3]
      } else { stop("mgd_type must be 'R' or 'sim'") }
        
      colnames(newcols) <- c("t.eta","t.chi")
      tmp <- cbind(tmp,newcols) # add theoretical sensitivity and specificity
        
      # now calculate the mutation-rate based fdr
      tmp <- mutate(tmp,tfdr = 1-phylosamp::prob_trans_mtml(
        chi = t.chi, eta = t.eta, rho = rho, M = M, R = 1))
        
    }
      
    else { stop("sens_spec_method must be specified and must be 'sim' or 'mutrate'") }
      
    # calculate the bias and error of the calculated fdr
    tmp <- mutate(tmp, bias = tfdr-fdr.sub)
    tmp <- mutate(tmp, err = abs(tfdr-fdr.sub))
      
    data_list[[i]] <- tmp # save new data
  }
    
  # save this data to file for faster plot manipulation
  save(data_list,file=outdir)
}


##' Function to plot theoretical versus simulated false discovery rate
##'
##' @param simdata saved simulation data with tfdr
##' @param rho_values vector of rho values
##' @param max_sim_size maximum simulation size allowed in output
##' @param filter_by data column to use for filtering values
##' @param filter_value value to use to filter data
##' @param xparam column data to put on x-axis of plot
##' @param yparam column data to put on y-axis of plot
##' @param axis_labels vector of length 2 where first value is x-axis label, second value is y-axis label
##' @param include_titles true/false on whether or not to include panel titles

plt.sim.data <- function(simdata,rho_values,max_sim_size,filter_by=NA,filter_value=NA,
                         xparam,yparam,axis_labels,include_titles=TRUE){
  
  load(simdata)
  
  # set up colors for plotting
  pal1 = brewer.pal(n = 9, "Blues")[2:9]
  pal2 = brewer.pal(n = 9, "Greens")[2:9]
  pal3 = brewer.pal(n = 9, "Oranges")[2:9]
  pal4 = brewer.pal(n = 9, "Greys")[2:9]
  pals <- list(pal1,pal2,pal3,pal4)
  
  # set up plots
  plots <- vector(mode = "list", length = length(rho_values))
  
  # set up mapping of rho_values to data_list indexes
  data_list_setup <- c(0.1,0.25,0.5,0.75)
  
  for (j in 1:length(rho_values)){
    
    # get the correct data_list for this value of rho
    i <- which(data_list_setup==rho_values[j])
    
    # filter data frame by specified values
    if (!is.na(filter_by)) { data_list[[i]] <- data_list[[i]] %>%
      filter(get(filter_by)>=filter_value) }
    
    # set up legend breaks for plot
    top_lim <- rho_values[j] * max_sim_size
    breaks <- c(0,top_lim*0.25,top_lim*0.5,top_lim*0.75,top_lim)
    
    # create a plot
    plots[[j]] <- ggplot() + 
      geom_point(data = data_list[[i]] %>% arrange(M),
                 aes(x = get(xparam), y = get(yparam), colour = M), shape = 3) +
      scale_color_gradientn(colours = pals[[i]],
                            breaks = breaks, labels = breaks, limits = c(0,top_lim)) +
      xlab(axis_labels[1]) + xlim(0,1) +
      ylab(axis_labels[2]) + ylim(0,1) +
      {if(include_titles)ggtitle(paste("Sampling proportion:",as.character(rho_values[j])))} +
      theme_classic() +
      theme(plot.title = element_text(size = 11, hjust = 0.65, face = "bold"),
            axis.title=element_text(size=10),
            axis.text=element_text(size=9),
            legend.title = element_text(size = 9),
            legend.title.align = 0.1,
            legend.key.size = unit(0.5,"cm"),
            legend.text = element_text(size = 8)) +
      geom_smooth(data = data_list[[i]] %>% arrange(M),
                  aes(x = get(xparam), y = get(yparam)),color="white",span=0.1) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray70")
  }
  
  # arrange plots based on variable number
  nrows <- floor(sqrt(length(rho_values)))
  do.call("grid.arrange", c(plots, nrow=nrows))

}


##' Function to get information for bias table
##'
##' @param saved_data previously-generated full data frame
##' @param var data column to use for table; 'bias' or 'error'
##' @param break_by data column to use for determining rows
##' @param row_breaks breakpoints to use for rows
##' @param row_names names to use for rows
##' @param filter_by data column to use for filtering values
##' @param filter_value value to use to filter data

fdr.make.tables <- function(saved_data,var,break_by,row_breaks,row_names,filter_by=NA,filter_value=NA){
  
  # load saved data
  load(saved_data)
  
  # set up a matrix
  # rows: categories of split variable
  # columns: usually rho categories - 0.1,0.25,0.5,0.75,all
  # extra column to store number of points in that row
  # extra row to store number of points in that column
  df <- matrix(0,length(row_breaks)+1,6)
  all_data <- data_frame()
  
  # fill in columns for each category
  for (i in 1:4){
    dat <- data_list[[i]]
    
    # subset the entire dataframe using the filtering value
    if (!is.na(filter_by)) { dat <- dat %>% filter(get(filter_by)>=filter_value) }
    
    # bind filtered dataset for this sampling proportion to all data
    all_data <- rbind(all_data,dat)
    
    for (j in 1:(length(row_breaks)-1)){
      # start (inclusive) to end (exclusive)
      if (j<(length(row_breaks)-1)) { 
        tmp <- dat %>% 
          filter(get(break_by)>=row_breaks[j]) %>% filter(get(break_by)<row_breaks[j+1]) }
      # start (inclusive) to end (also inclusive) for last value
      else { tmp <- dat %>% 
        filter(get(break_by)>=row_breaks[j]) %>% filter(get(break_by)<=row_breaks[j+1]) }
      
      df[j,i] <- mean(tmp[[var]])
    }
    
    df[length(row_breaks),i] <- mean(dat[[var]]) # value over all rows
    df[length(row_breaks)+1,i] <- length(dat[[var]]) # number of points in this row
  }
  
  # determine the averages for all values of rho
  for (j in 1:length(row_breaks)){
    # start (inclusive) to end (exclusive)
    if (j<length(row_breaks)) { 
      tmp <- all_data %>% 
        filter(get(break_by)>=row_breaks[j]) %>% filter(get(break_by)<row_breaks[j+1]) }
    # start (inclusive) to end (also inclusive) for last value
    else { tmp <- all_data %>% 
      filter(get(break_by)>=row_breaks[j]) %>% filter(get(break_by)<=row_breaks[j+1]) }
    
    df[j,5] <- mean(tmp[[var]]) # value over all rho
    df[j,6] <- length(tmp[[var]]) # number of points in this row
  }
  
  df[length(row_breaks),5] <- mean(all_data[[var]])
  df[length(row_breaks),6] <- length(all_data[[var]])
  df[length(row_breaks)+1,5] <- sum(df[length(row_breaks)+1,1:4])
  df[length(row_breaks)+1,6] <- NA
  
  # turn this into nice-looking data frames
  df <- data.frame(df)
  names(df) <- c("0.10","0.25","0.50","0.75","Overall","N")
  
  df1 <- df[1:length(row_breaks),] %>% select(-N) %>% mutate_if(is.numeric, ~sprintf("%.4f",.))
  df1 <- rbind(df1,format(select(df[length(row_breaks)+1,],-N),big.mark=","))
  row.names(df1) = c(row_names,"All","N")
  dfN <- format(df %>% select(N),big.mark=",")
  
  # return dataframe for any desired manipulation
  return(cbind(df1,dfN))
  
}


##' Function to get information for bias table
##'
##' @param saved_data previously-generated full data frame
##' @param var data column to use for table; 'bias' or 'err'
##' @param rho rho value to use for this table
##' @param row_break_by data column to use for determining rows
##' @param row_breaks breakpoints to use for rows
##' @param row_names names to use for rows
##' @param col_break_by data column to use for determining rows
##' @param col_breaks breakpoints to use for rows
##' @param col_names names to use for rows
##' @param filter_by data column to use for filtering values
##' @param filter_value value to use to filter data


fdr.make.tables.var <- function(saved_data,var,rho,row_break_by,row_breaks,row_names,
                                col_break_by,col_breaks,col_names,filter_by=NA,filter_value=NA){
  
  # load saved data
  load(saved_data)
  
  # set up a matrix
  # extra column to store number of points in that row
  # extra row to store number of points in that column
  df <- matrix(0,length(row_breaks)+1,length(col_breaks)+1)
  
  # get the data for a specific value of rho
  if (rho==0.1) { dat <- data_list[[1]]
  } else if (rho==0.25) { dat <- data_list[[2]]
  } else if (rho==0.5) { dat <- data_list[[3]]
  } else if (rho==0.75) { dat <- data_list[[4]]
  } else { stop("not a valid value of rho") }
  
  # calculate error and bias of sensitivity and specificy
  # in case we want to plot these values instead of fdr error
  #dat <- mutate(dat, eta.err = abs(t.eta-eta.full))
  #dat <- mutate(dat, chi.err = abs(t.chi-chi.full))
  #dat <- mutate(dat, eta.bias = (t.eta-eta.full))
  #dat <- mutate(dat, chi.bias = (t.chi-chi.full))
  
  # subset the entire dataframe using the filtering value
  if (!is.na(filter_by)) { dat <- dat %>% filter(get(filter_by)>=filter_value) }
  
  for (i in 1:(length(col_breaks)-1)){
    # start (inclusive) to end (exclusive)
    if (i<(length(col_breaks)-1)) { 
      cval <- dat %>% 
        filter(get(col_break_by)>=col_breaks[i]) %>% filter(get(col_break_by)<col_breaks[i+1]) }
    # start (inclusive) to end (also inclusive) for last value
    else { cval <- dat %>% 
      filter(get(col_break_by)>=col_breaks[i]) %>% filter(get(col_break_by)<=col_breaks[i+1]) }
    
    for (j in 1:(length(row_breaks)-1)){
      # start (inclusive) to end (exclusive)
      if (j<(length(row_breaks)-1)) { 
        rval <- cval %>% 
          filter(get(row_break_by)>=row_breaks[j]) %>% filter(get(row_break_by)<row_breaks[j+1]) }
      # start (inclusive) to end (also inclusive) for last value
      else { rval <- cval %>% 
        filter(get(row_break_by)>=row_breaks[j]) %>% filter(get(row_break_by)<=row_breaks[j+1]) }
      
      df[j,i] <- mean(rval[[var]])
    }
    
    df[length(row_breaks),i] <- mean(cval[[var]]) # value over all rows
    df[length(row_breaks)+1,i] <- length(cval[[var]]) # number of points in this row
  }
  
  # get values over all columns
  for (j in 1:(length(row_breaks)-1)){
    # start (inclusive) to end (exclusive)
    if (j<length(row_breaks)) { 
      tmp <- dat %>% 
        filter(get(row_break_by)>=row_breaks[j]) %>% filter(get(row_break_by)<row_breaks[j+1]) }
    # start (inclusive) to end (also inclusive) for last value
    else { tmp <- dat %>% 
      filter(get(row_break_by)>=row_breaks[j]) %>% filter(get(row_break_by)<=row_breaks[j+1]) }
    
    
    df[j,length(col_breaks)] <- mean(tmp[[var]]) # value over all columns
    df[j,length(col_breaks)+1] <- length(tmp[[var]]) # number of points in this column
  }
  
  # get summary values over all rows and columns
  df[length(row_breaks),length(col_breaks)] <- mean(dat[[var]])
  df[length(row_breaks),length(col_breaks)+1] <- length(dat[[var]])
  df[length(row_breaks)+1,length(col_breaks)] <- sum(df[length(row_breaks)+1,1:(length(col_breaks)-1)])
  df[length(row_breaks)+1,length(col_breaks)+1] <- NA
  
  # turn this into nice-looking data frames
  df <- data.frame(df)
  names(df) <- c(col_names,"All","N")
  
  df1 <- df[1:length(row_breaks),] %>% select(-N) %>% mutate_if(is.numeric, ~sprintf("%.4f",.))
  df1 <- rbind(df1,format(select(df[length(row_breaks)+1,],-N),big.mark=","))
  row.names(df1) = c(row_names,"All","N")
  dfN <- format(df %>% select(N),big.mark=",")

  # return dataframe for any desired manipulation
  return(cbind(df1,dfN))
  
}


##' Function to make error and count heatmaps given sensitivity and specificity
##'
##' @param saved_data previously-generated full data frame
##' @param rho rho value to use for the plots
##' @param sens sensitivity values to use as bins
##' @param spec specificity values to use as bins
##' @param sens_breaks y-axis breaks
##' @param spec_breaks x-axis breaks

make.error.heatmaps <- function(saved_data,rho,sens,spec,sens_breaks,spec_breaks){

  load(saved_data)
  # get the data for a specific value of rho
  if (rho==0.1) { dat <- data_list[[1]]
  } else if (rho==0.25) { dat <- data_list[[2]]
  } else if (rho==0.5) { dat <- data_list[[3]]
  } else if (rho==0.75) { dat <- data_list[[4]]
  } else { stop("not a valid value of rho") }
  
  # definte sensitivity and specificity bins
  vals <- expand.grid(sens_start=sens[1:(length(sens)-1)],spec_start=spec[1:(length(spec)-1)])
  sens_stop <- rep(NA,length(vals$sens_start))
  for (i in 1:length(vals$sens_start)) { sens_stop[i] <- sens[which(sens==vals$sens_start[i])+1] }
  spec_stop <- rep(NA,length(vals$spec_start))
  for (i in 1:length(vals$spec_start)) { spec_stop[i] <- spec[which(spec==vals$spec_start[i])+1] }
  vals$sens_stop <- sens_stop
  vals$spec_stop <- spec_stop
  
  # define a function that filters the dataset and calculates the error and counts
  # allows us to avoid a for loop
  get.error.counts <- function(df,sens_start,sens_stop,spec_start,spec_stop){
    
    # if sens or spec is equal to one
    # capture this is the upper bound of that bin
    if (!(sens_stop==1 || spec_stop==1)){
      df <- df %>% filter(chi.full>=spec_start) %>% filter(chi.full<spec_stop)
      df <- df %>% filter(eta.full>=sens_start) %>% filter(eta.full<sens_stop)
    } else if (sens_stop==1 && spec_stop==1) {
      df <- df %>% filter(chi.full>=spec_start) %>% filter(chi.full<=spec_stop)
      df <- df %>% filter(eta.full>=sens_start) %>% filter(eta.full<=sens_stop)
    } else if (sens_stop==1) {
      df <- df %>% filter(chi.full>=spec_start) %>% filter(chi.full<spec_stop)
      df <- df %>% filter(eta.full>=sens_start) %>% filter(eta.full<=sens_stop)
    } else if (spec_stop==1) {
      df <- df %>% filter(chi.full>=spec_start) %>% filter(chi.full<=spec_stop)
      df <- df %>% filter(eta.full>=sens_start) %>% filter(eta.full<sens_stop)
    } else { stop("something went wrong!") }
    
    c <- length(df$err)
    if (c > 0) { err <- mean(df$err) }
    else { err <- NA }
    return(c(err,c))
    
  }
  
  res <- t(apply(vals,1,function(x) 
    get.error.counts(df=dat,sens_start=x["sens_start"],sens_stop=x["sens_stop"],
                     spec_start=x["spec_start"],spec_stop=x["spec_stop"])))
  colnames(res) <- c("mean_err","count")
  res <- cbind(vals,res)
  
  # change specificity to 1-specificity
  res$spec_start = 1-res$spec_start
  res$spec_stop = 1-res$spec_stop
  
  # change error to % for easy interpretation
  res$mean_err <- res$mean_err * 100
  
  # make the plots
  p_err <- ggplot(data = res, aes(x=spec_start, y=sens_start, fill=mean_err)) + 
    geom_tile() +
    xlab("1-specificity") +
    ylab("sensitivity") + 
    scale_x_continuous(expand = c(0, 0), breaks = spec_breaks) +
    scale_y_continuous(expand = c(0, 0), breaks = sens_breaks) +
    scale_fill_gradient(low="white",high="darkblue",na.value = "grey90",name="% error",limits=c(0,12),breaks=c(3,6,9)) +
    theme(legend.direction = "vertical",
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.8,"line"),
          legend.text = element_text(size = 9),
          legend.justification = 'right',
          axis.title=element_text(size=10,face="bold"))
  
  p_count <- ggplot(data = res, aes(x=spec_start, y=sens_start, fill=count)) + 
    geom_tile() +
    xlab("1-specificity") +
    ylab("sensitivity") + 
    scale_x_continuous(expand = c(0, 0), breaks = spec_breaks) +
    scale_y_continuous(expand = c(0, 0), breaks = sens_breaks) +
    scale_fill_gradient(low="grey90",high="coral3",na.value = "coral4",name="count",limits=c(0,200),breaks=c(50,100,150)) +
    theme(legend.direction = "vertical",
          legend.title = element_text(size = 9),
          legend.key.size = unit(0.8,"line"),
          legend.text = element_text(size = 9),
          legend.justification = 'right',
          axis.title=element_text(size=10,face="bold"))
  
  grid.arrange(p_err,p_count,nrow=2)
  
}

##' Function to make histograms showing deviation between calculated and simulated values
##'
##' @param saved_data previously-generated full data frame
##' @param single_value should only the optimal sens/spec be plotted for each simulation

plot.error.hist <- function(saved_data,single_value){
  
  # load saved simulation data
  load(saved_data)
  
  fdr_plots <- vector(mode = "list", length = 4)
  eta_plots <- vector(mode = "list", length = 4)
  chi_plots <- vector(mode = "list", length = 4)
  
  pal1 = brewer.pal(n = 9, "Blues")[7]
  pal2 = brewer.pal(n = 9, "Greens")[7]
  pal3 = brewer.pal(n = 9, "Oranges")[6]
  pal4 = brewer.pal(n = 9, "Greys")[8]
  colors <- c(pal1,pal2,pal3,pal4)
  
  for (i in 1:4){
    dat <- data_list[[i]]
    
    if (single_value==TRUE){
      tmp <- dat %>% mutate(corner.dist = sqrt(((1-t.chi)^2)+((1-t.eta)^2)))
      tmp <- tmp %>% group_by(sim) %>% filter(corner.dist == min(corner.dist)) %>% ungroup()
      tmp <- data.frame(tmp)
      dat <- tmp
    }
    
    fdr_plots[[i]] <- ggplot(dat,aes(x=tfdr-fdr.sub)) +
      geom_histogram(aes(y=..count../sum(..count..)),binwidth=0.01,position="identity",
                     fill=colors[i],color='white',size=0.2) +
      scale_x_continuous(expand = c(0,0),limits = c(-0.15,0.15),breaks = c(-0.10,0.00,0.10)) +
      scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
      theme_classic() + 
      theme(legend.position='none',
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=8),
            axis.text.y = element_text(size=8))
    
    eta_plots[[i]] <- ggplot(dat,aes(x=t.eta-eta.full)) +
      geom_histogram(aes(y=..count../sum(..count..)),binwidth=0.01,position="identity",
                     fill=colors[i],color='white',size=0.2) +
      scale_x_continuous(expand = c(0,0),limits = c(-0.15,0.15),breaks = c(-0.10,0.00,0.10)) +
      scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
      theme_classic() + 
      theme(legend.position='none',
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=8),
            axis.text.y = element_text(size=8))
    
    chi_plots[[i]] <- ggplot(dat,aes(x=t.chi-chi.full)) +
      geom_histogram(aes(y=..count../sum(..count..)),binwidth=0.01,position="identity",
                     fill=colors[i],color='white',size=0.2) +
      scale_x_continuous(expand = c(0,0),limits = c(-0.15,0.15),breaks = c(-0.10,0.00,0.10)) +
      scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
      theme_classic() + 
      theme(legend.position='none',
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=8),
            axis.text.y = element_text(size=8))
  }
  
  # return all plots in order
  plots <- c(fdr_plots,eta_plots,chi_plots)
  return(plots)
  
}

##' Function to make histograms showing deviation between calculated and simulated values
##' this time split up by sample size
##'
##' @param saved_data previously-generated full data frame
##' @param rho proportion sampled
##' @param col_break_by values to use to split up data
##' @param col_breaks breakpoints to use for rows

plot.error.hist.breaks <- function(saved_data,rho,col_break_by,col_breaks){
  
  # load saved simulation data
  load(saved_data)
  
  # get the data for a specific value of rho
  if (rho==0.1) { 
    dat <- data_list[[1]]
    hist_color <- brewer.pal(n = 9, "Blues")[7]
  } else if (rho==0.25) { 
    dat <- data_list[[2]]
    hist_color <- brewer.pal(n = 9, "Greens")[7]
  } else if (rho==0.5) { 
    dat <- data_list[[3]]
    hist_color <- brewer.pal(n = 9, "Oranges")[6]
  } else if (rho==0.75) { 
    dat <- data_list[[4]]
    hist_color <- brewer.pal(n = 9, "Greys")[8]
  } else { stop("not a valid value of rho") }
    
  plots <- vector(mode = "list", length = (length(col_breaks)-1))
    
  for (j in 1:(length(col_breaks)-1)){
    # start (inclusive) to end (exclusive)
    if (j<(length(col_breaks)-1)) { 
      cval <- dat %>% 
        filter(get(col_break_by)>=col_breaks[j]) %>% filter(get(col_break_by)<col_breaks[j+1]) }
    # start (inclusive) to end (also inclusive) for last value
    else { cval <- dat %>% 
      filter(get(col_break_by)>=col_breaks[j]) %>% filter(get(col_break_by)<=col_breaks[j+1]) }
    
    plots[[j]] <- ggplot(cval,aes(x=t.chi-chi.full)) +
      geom_histogram(aes(y=..count../sum(..count..)),binwidth=0.01,position="identity",
                     fill=hist_color,color='white',size=0.2) +
      scale_x_continuous(expand = c(0,0),limits = c(-0.25,0.25),
                         breaks = c(-0.20,-0.10,0.00,0.10,0.20)) +
      scale_y_continuous(expand = c(0,0),limits = c(0,0.4)) +
      theme_classic() + 
      theme(legend.position='none',
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=8),
            axis.text.y = element_text(size=8))
  }
  
  return(plots)
  
}