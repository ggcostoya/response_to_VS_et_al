
## Re-runnning simulations for VS et al. response ##

## Packages ----

library(mvtnorm)
library(tidyverse)
library(parallel)
library(foreach)
library(doParallel)
library(data.table)
library(DescTools)

registerDoParallel(detectCores()) # Register parallel cores

## Functions to generate populations ----

# modified randomized tpts function
randomize_tpts <- function(tpts, corr, s){
  
  # build basic genetic correlations matrix
  gmtx <- diag(length(tpts))
  
  # Simulating data from a multivariate analysis
  if(corr %in% c("gsto", "both")){gmtx[5,1] <- gmtx[1,5] <- gmtx[5,2] <- gmtx[2,5] <- s}
  if(corr %in% c("gsto", "both")){gmtx[5,4] <- gmtx[4,5] <- -s}
  
  # tde
  if(corr %in% c("tde", "both")){gmtx[5,3] <- gmtx[3,5] <- s}
  
  # define change in traits due to genetic correlations
  # note that the change is now calculated via 
  # note that suppressWarnings is added to suppress messages related to the
  # definitiveness of the matrix when S = 0.75, but thisissue is not ignored. 
  change <- suppressWarnings(rmvnorm(1, rep(0, 5), gmtx))
  
  # new tpt values
  new_tpts <- as.vector(tpts + change)
  
  # correct for impossible values
  
  # mid
  new_tpts[2] <- ifelse(new_tpts[2] > new_tpts[3] - 1, new_tpts[3] - 1, new_tpts[2])
  
  # ctmin
  new_tpts[1] <- ifelse(new_tpts[1] > new_tpts[2] - 1, new_tpts[1] - 1, new_tpts[1])
  
  # ctmax
  new_tpts[4] <- ifelse(new_tpts[4] < new_tpts[3] + 1, new_tpts[3] + 1, new_tpts[4])
  
  # round final values to two decimal places
  new_tpts <- round(new_tpts, digits = 2)
  
  return(new_tpts)
  
  
}

# build thermal-performance dataset from tpts
build_tpd <- function(tpts){
  
  # temperature
  t <- round(tpts[1:4], 1)
  
  # performance
  pmax <- tpts[5]
  p <- c(0, pmax * 0.6, pmax, 0)
  
  return(data.frame(t,p))
  
}

# generate a tpc from thermal-performance dataset
generate_tpc <- function(tpd){
  
  # temperature
  t <- round(seq(0, 50, by = 0.1), digits = 1)
  
  # holder performance
  p <- rep(NA, length(t))
  
  # set zero p for below and above ctmin and ctmax
  p[1:which(t == tpd$t[1])] <- p[which(t == tpd$t[4]):length(t)] <- 0
  
  # get indexes for critical points
  ctmin <- which(t == tpd$t[1])
  mid <- which(t == tpd$t[2])
  topt <- which(t == tpd$t[3])
  ctmax <- which(t == tpd$t[4])
  
  # set performance values for each vector
  p[ctmin:mid] <- seq(0, tpd$p[2], length.out = length(p[ctmin:mid]))
  p[mid:topt] <- seq(tpd$p[2], tpd$p[3], length.out = length(p[mid:topt]))
  p[topt:ctmax] <- seq(tpd$p[3], 0, length.out = length(p[topt:ctmax]))
  
  return(tibble(t = t, p = p))
  
}

# generate a population combining the previous 3 functions
generate_population <- function(reps, n, tpts, corr, s){
  
  # Loop for all replicates
  pops <- foreach(i = 1:reps) %do% {
    
    # loop for each replicate
    pop <- foreach(j = 1:n,
                   combine = rbind,
                   .packages = c("tidyverse","mvtnorm"),
                   .export = c("randomize_tpts", "build_tpd", "generate_tpc")) %dopar% {
                     
                     tpc <- nest(generate_tpc(build_tpd(randomize_tpts(tpts = tpts, corr = corr, s = s))))
                     
                   }
    
    # reshape population object
    do.call(rbind, pop) %>% mutate(tpc = data) %>% select(tpc)
    
  }
  
  return(pops)
}

## Set base to generate populations parameters ----

reps <- 10 # Number of replicates

topt <- 25 # Thermal optimum

pmax <- 10 # Maximum performance

ctmin <- 10 # Critical thermal minimum

ctmax <- 32.5 # Critical thermal maximum

midpoint <- mean(c(topt, ctmin)) # Fourth polygon point

tpts <- c(ctmin, midpoint, topt, ctmax, pmax)

## Generate new populations ----

# With S = 0.75
none_0.75 <- generate_population(reps, tpts = tpts, n = 50, corr = "none", s = 0.75)
gsto_0.75 <- generate_population(reps, tpts = tpts, n = 50, corr = "gsto", s = 0.75)
tde_0.75 <- generate_population(reps, tpts = tpts, n = 50, corr = "tde", s = 0.75)
both_0.75 <- generate_population(reps, tpts = tpts, n = 50, corr = "both", s = 0.75)

# With S = 0.5
none_0.5 <- generate_population(reps, tpts = tpts, n = 50, corr = "none", s = 0.5)
gsto_0.5 <- generate_population(reps, tpts = tpts, n = 50, corr = "gsto", s = 0.5)
tde_0.5 <- generate_population(reps, tpts = tpts, n = 50, corr = "tde", s = 0.5)
both_0.5 <- generate_population(reps, tpts = tpts, n = 50, corr = "both", s = 0.5)

## Function to generate temperature sequences ----

generate_tseq <- function(reps, gen, days, mean_0, sd_0, mean_c, sd_c, burnin){
  
  # sequence of mean temperature values
  means <- seq(mean_0, mean_0 + mean_c, length.out = gen)
  burnin_means <- rep(mean_0, burnin)
  means <- c(burnin_means, means)
  
  # determine standard deviation from breath
  sds <- seq(sd_0, sd_0 + sd_c, length.out = gen)
  burnin_sds <- rep(sd_0, burnin)
  sds <- c(burnin_sds, sds)
  
  # loop for all replicates
  tseqs <- foreach(i = 1:reps) %do% {
    
    # loop for each replicate
    foreach(j  = 1:(gen + burnin)) %do% {
      
      round(rnorm(n = days, mean = means[j], sd = sds[j]), 1)}
    
  }
  
  return(tseqs)
}

## Generate temperature sequences ----

reps <- 10 # Number of replicates

gen <- 80 # Number of generations

days <- 150 # Number of days per

mean_0 <- 25 # Initial env. mean temperature

sd_0 <- 1 # Initial SD

burnin <- 5 # Burn-in period

sd_rel_factor <- 0.15 # SD change relative to Mean change factor

rcp8.5 <- generate_tseq(reps, gen, days, mean_0, sd_0, burnin, mean_c = 2.96, sd_c = sd_rel_factor * 2.96)

## Define simulation function ----

# parallel simulation function in order to run multiple simulations at the same time
parallel_simulation <- function(pops, tseqs, corr_name, tseq_name){
  
  # define n_k
  n_k <- nrow(pops[1])
  
  # simulation function
  simulation <- function(pop, tseq){
    
    # get the initial population size and set that as cap
    popsize_cap <- nrow(pop)
    
    # set initial population as current population
    current_pop <- pop %>% mutate(gen = rep(0, nrow(.)), r = rep(0, nrow(.)))
    
    # object to store list of population data each generation
    sim_data <- list(current_pop)
    
    # start of the generation loop
    for(i in 1:length(tseq)){
      
      # select the temperature sequence belonging to the current generation & index it
      idx_t <- round((tseq[[i]]) / 0.1)
      
      # loop to determine next population
      next_pop <- foreach(k = 1:nrow(current_pop)) %do% {
        
        # extract individual tpc
        tpc <- sim_data[[i]]$tpc[[k]]
        
        # determine daily performance at the env temperature through indexing
        p_et <- tpc$p[idx_t]
        
        # determine daily probability of survival
        p_s <- 1/(1 + exp(-(-5 + 1 * p_et)))
        
        # determine daily survival
        days_s <- as.numeric(runif(length(p_s)) < p_s)
        
        # determine day of death (minimum day with survival = 0)
        min_idx <- suppressWarnings(min(which(days_s == 0)))
        
        # reshape days survived after individual death
        if(min_idx != Inf){days_s[(min_idx + 1):length(days_s)] <- 0}
        
        # determine the sum of days alive
        days_alive <- sum(days_s, na.rm = T)
        
        # determine offspring produced based on days alive
        n_off <- floor(days_alive/10)
        
        # record reproductive output
        sim_data[[i]]$r[[k]] <- n_off
        
        # nest individual tpc data
        tpc <- tpc %>% nest(tpc = c(t,p))
        
        # generate individual offspring
        off <- do.call("rbind", replicate(n_off, tpc, simplify = FALSE))
        
      }
      
      # reshape the next population object
      next_pop <- as_tibble(rbindlist(next_pop))
      
      # break loop if length of the next population is zero (extinction)
      if(nrow(next_pop) == 0){break}
      
      # if the population exceeds the carrying capacity select up to it
      if(nrow(next_pop) > popsize_cap){next_pop <- sample_n(next_pop, popsize_cap)}
      
      # update current population
      current_pop <- next_pop
      
      # add generation and reproduction columns to new offspring
      next_pop <- next_pop %>% mutate(gen = rep(i,), r = rep(0,))
      
      # bind to simulation data
      sim_data <- c(sim_data, list(next_pop))
      
    }
    
    # reshape simulation data output
    sim_data <- suppressWarnings(as_tibble(rbindlist(sim_data)))
    
    # set starting values for TPTs
    topt <- ctmin <- ctmax <- pmax <- mid <- NA
    
    # generate holder tpt columns
    sim_data <- mutate(sim_data, topt, ctmin, ctmax, pmax, mid)
    
    # loop to assign tpt values from TPC
    for(i in 1:nrow(sim_data)){
      
      # extract TPC
      tpc <- sim_data$tpc[[i]]
      
      # extract tpts
      sim_data$pmax[i] <- max(tpc$p)
      sim_data$topt[i] <- tpc$t[tpc$p == max(tpc$p)]
      sim_data$ctmin[i] <- max(tpc$t[tpc$p == 0 & tpc$t < tpc$t[tpc$p == max(tpc$p)]])
      sim_data$ctmax[i] <- min(tpc$t[tpc$p == 0 & tpc$t > tpc$t[tpc$p == max(tpc$p)]])
      mid <- suppressWarnings(tpc$t[tpc$p == Closest(tpc$p, max(tpc$p) / 2)])
      sim_data$mid[i] <- mid[1]
      
    }
    
    # indicate columns to keep
    sim_data <- sim_data %>% select(topt, pmax, ctmin, ctmax, mid, gen, r)
    
    return(sim_data)
    
  }
  
  # define packages needed
  packages <- c("tidyverse", 'parallel', "doParallel",
                "foreach", "data.table", "DescTools")
  
  # get raw simulation output
  raw_sim <- foreach(i = 1:length(pops), .packages = packages) %dopar% {
    
    foreach(j = 1:length(tseqs)) %do% {
      
      sim <- simulation(pop = pops[[i]], tseq = tseqs[[j]])
      
      sim %>% mutate(pop_n = rep(i, nrow(.)), tseq_n = rep(j, nrow(.)))
      
    }
  }
  
  # post process raw sim object
  sim <- do.call(rbind, do.call(c, raw_sim)) %>%
    mutate(corr = rep(corr_name, nrow(.)),
           tseq = rep(tseq_name, nrow(.)),
           n_k = rep(n_k, nrow(.)))
  
  return(sim)
  
}

## Run simulations ----

# run and bind simulations with S = 0.75
both_0.75_sim <- parallel_simulation(pops = both_0.75, tseqs = rcp8.5, corr_name = "GSTO + TDE", tseq_name = "RCP 8.5")
none_0.75_sim <- parallel_simulation(pops = none_0.75, tseqs = rcp8.5, corr_name = "None", tseq_name = "RCP 8.5")
tde_0.75_sim <- parallel_simulation(pops = tde_0.75, tseqs = rcp8.5, corr_name = "TDE", tseq_name = "RCP 8.5")
gsto_0.75_sim <- parallel_simulation(pops = gsto_0.75, tseqs = rcp8.5, corr_name = "GSTO", tseq_name = "RCP 8.5")
sims_0.75 <- rbind(none_0.75_sim, gsto_0.75_sim, tde_0.75_sim, both_0.75_sim)

# run and bind simulations with S = 0.5
both_0.5_sim <- parallel_simulation(pops = both_0.5, tseqs = rcp8.5, corr_name = "GSTO + TDE", tseq_name = "RCP 8.5")
none_0.5_sim <- parallel_simulation(pops = none_0.5, tseqs = rcp8.5, corr_name = "None", tseq_name = "RCP 8.5")
tde_0.5_sim <- parallel_simulation(pops = tde_0.5, tseqs = rcp8.5, corr_name = "TDE", tseq_name = "RCP 8.5")
gsto_0.5_sim <- parallel_simulation(pops = gsto_0.5, tseqs = rcp8.5, corr_name = "GSTO", tseq_name = "RCP 8.5")
sims_0.5 <- rbind(none_0.5_sim, gsto_0.5_sim, tde_0.5_sim, both_0.5_sim)

## Plotting ----

# loadding and processing original simulation results
load("data/sum_sim_results.RData")
sims_o_0.75 <- sum_sim_results %>% filter(tseq == "high", n0_k == 50) %>%
  rename(avg_pop_size = n) %>% 
  mutate(for_what = "OP", s = 0.75) %>% 
  dplyr::select(for_what, s, corr, gen, avg_pop_size)
sims_o_0.75$corr <- as.factor(sims_o_0.75$corr)
levels(sims_o_0.75$corr) <- list(`GSTO + TDE` = "both", GSTO = "gsto", None = "none", TDE = "tde")

# adding columns for joint plotting 
sims_0.5 <- sims_0.5 %>% mutate(for_what = "TN", s = 0.5)
sims_0.75 <- sims_0.75 %>% mutate(for_what = "TN", s = 0.75)

# combining simulation results
sims <- rbind(sims_0.5, sims_0.75)

# generate data with all generations
gen <- seq(0,85,1)
pop_n <- seq(1,10,1)
tseq_n <- seq(1,10,1)
corrs <- c("None", "GSTO", "TDE", "GSTO + TDE")
for_what <- "TN"
s <- c(0.5, 0.75)
all_gens <- expand.grid(for_what, s, corrs, gen, pop_n, tseq_n)
colnames(all_gens) <- c("for_what", "s", "corr", "gen", "pop_n", "tseq_n")

# get population sizes
n_plot <- sims %>% group_by(for_what, s, corr, gen, pop_n, tseq_n) %>% 
  summarise(pop_size = n())

# merge with all generations data & set N = 0 in generations where no data was recorded
n_plot <- merge(n_plot, all_gens, by = c("for_what", "s", "corr", "gen", "pop_n", "tseq_n"), all = TRUE)
n_plot$pop_size <- ifelse(is.na(n_plot$pop_size), 0, n_plot$pop_size)

# find average population size by correlation and generation
n_plot <- n_plot %>% group_by(for_what, s, corr, gen) %>% summarise(avg_pop_size = mean(pop_size)) 

# bind with original simulation results
n_plot <- rbind(n_plot, sims_o_0.75)

# set correlations as a factor and define colors
n_plot$corr <- factor(n_plot$corr, levels = c("None", "GSTO","TDE","GSTO + TDE"))
colors <- c("None" = "#009E73", "GSTO" = "#56B4E9","TDE" = "#CC79A7", "GSTO + TDE" = "#D55E00")

# plot
n_plot %>% 
  mutate(where = ifelse(for_what == "TN", "Re-run simulation", "Original simulation")) %>%
  mutate(s_name = paste("Correlation strength (S) =", s)) %>%
  #filter(s == 0.75) %>% # Filter for only original vs re-run comparison
  filter(for_what == "TN") %>% # Filter for only re-run comparison between S
  ggplot(aes(x = gen, y = avg_pop_size)) +
  geom_rect(xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf, alpha = 0.5, fill = "lightgray") +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 1.25) +
  geom_line(aes(col = corr), linewidth = 1.25, alpha = 0.75) +
  scale_x_continuous(expand = c(0,0)) +
  scale_color_manual(values = colors) +
  xlab("Generation") +
  ylab("Mean population size (N)") +
  #facet_grid(cols = vars(where)) + # add or remove for original vs re-run comparisons
  facet_grid(cols = vars(s_name)) + # add or remove for S comparisons
  #facet_wrap(~where*s_name) + # add or remove to show all results
  theme_minimal() +
  theme(axis.line = element_line(),
        axis.ticks = element_line(),
        legend.position = c(0.15, 0.25),
        legend.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.5), 
        strip.text = element_text(size = 12)) +
  #ggtitle("Correlation strenght (S) = 0.75") # title for original vs re-run comparisons
  ggtitle("Using reviewer's methodology") # title for S comparisons














