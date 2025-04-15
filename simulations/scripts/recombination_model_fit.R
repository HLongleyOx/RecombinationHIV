########################################################################################
## This code analyzes a dataset describing linkage over time for diverse sites in simulated sequence data
## Uses method from Romero and Feder 2024 to estimate recombination rates
########################################################################################


### Load libraries ####
library(tidyverse)    
library(ggpubr)      
library(RColorBrewer) 
library(minpack.lm)   

#Helper functions
#function to fit data to curve as described in Romero and feder paper
f <- function(x, c0, c1, c2) {
  c0 + c1 * (1 - exp(-c2 * x))
}

#Set colour palette
cbPalette <- brewer.pal(6, "Set1")

#Function for labelling on plots
scientific_10 <- function(x) {
  parse(text = gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Load in simulation data
setwd("linkage") # Change working directory to where linkage files are stored. Replace file path if necessary
recom = c("7.5e-06", "1e-05", "2.5e-05", "5e-05")  # Different recombination rates to analyze
df <- data.frame()  # Initialise empty dataframe to store results

# Process each recombination rate file
for (file in recom) {
  # Read in file and format variables
  sel_linkage <- 
    read_csv(paste0("linkage_", file, ".csv"), show_col_types = FALSE) %>%
    mutate(
      days = as.numeric(sub("_([0-9]+)_", "\\1", date)),  # Extract days from date string
      siteA = as.numeric(siteA),  # Convert site positions to numeric
      siteB = as.numeric(siteB),
      pA = as.numeric(pA),        # Convert allele frequencies to numeric
      pB = as.numeric(pB),
      ddash = as.numeric(ddash),  # Convert linkage metrics to numeric
      ld = as.numeric(ld)
    ) %>%
    drop_na() %>%               # Remove rows with NA values
    group_by(ID, siteA, siteB, baseA, baseB) %>%  # Group by unique site pairs
    filter(ld >= 0.0001)        # Filter for meaningful linkage disequilibrium values
  
  
  # Only analyze data points when LD first exceeds 0.2 and all time points after 
  sel_first_timepoint <-  
    sel_linkage %>% 
    filter(ddash > 0.2) %>%     # Find site pairs with strong LD (D' > 0.2)
    group_by(ID, siteA, siteB, baseA, baseB) %>% 
    filter(days == min(days)) %>%  # Get first occurrence of strong LD
    mutate(first_day = days) %>%   # Mark the day when LD was first detected
    ungroup() %>% 
    select(siteA, siteB, first_day, ID) %>%  # Select key columns
    right_join(sel_linkage) %>%     # Join back with all LD data
    filter(days >= first_day)  %>%  # Keep only timepoints after LD was first detected
    group_by(ID, siteA, siteB, baseA, baseB) %>%
    arrange(days) %>%             
    mutate(
      last_days = days - lag(days),  # Time since previous measurement
      dist = abs(siteA - siteB),     # Physical distance between sites
      ddash_lag = lag(ddash),        # LD at previous timepoint
      lagPa = lag(pA),               # Previous allele frequencies
      lagPb = lag(pB),
      distt = (last_days) * dist,    # Time-scaled distance
      ddratio = -1 * log(ddash / ddash_lag),  # Calculate LD decay ratio (log scale)
      first_ld = first(ddash)        # Record initial LD value
    ) %>%
    drop_na()  
  
  
  # Fits the recombination model to data filtered for different maximum time-scaled distances
  for (d in seq(20000, 150000, by = 10000)) {
    # Filter data by time-scaled distance threshold
    sel <- sel_first_timepoint %>% 
      drop_na() %>% 
      filter(distt < d) %>%
      mutate(grp = round(distt / 1) * 1)  # Group by rounded time-scaled distance (set at size 1 so no meaningful groups)
    
    # Calculate mean LD decay ratio for each distance group
    x = sel %>% 
      group_by(grp) %>% 
      summarise(mn = mean(ddratio), N = n())  
    
    a = x$mn      # LD decay values
    b = x$grp     # Distance groups
    n = sqrt(x$N) # Weight by square root of sample size
    
    # Non-linear least squares fit
    fit <-
      nlsLM(
        a ~ f(b, c0, c1, c2),  
        start = list(c0 = 0.2, c1 = 0.3, c2 = 0.0001),  
        weights = n  # Weight by sample size
      )
    
    # Calculate recombination rate estimate (product of fitted parameters)
    p = ((coef(fit)[[2]] * coef(fit)[[3]]))
    print(p)  # Print the estimated recombination rate
    
    # Store results
    df = rbind(df, data.frame(
      r = p,                  # Estimated recombination rate
      dist = d,               # Maximum time-scaled distance used
      true_recom = file       # True recombination rate from simulation
    ))
  }
}

# Create plot of results

sim_plot <-
  ggplot(
    df %>% filter(true_recom %in% c("1e-05", "2.5e-05", "5e-05", "7.5e-06")),
    aes(
      x = dist,
      y = r,
      group = as.factor(true_recom),
      colour = as.factor(true_recom)
    )
  ) + 
  geom_line() +  
  theme_classic() +  
  scale_y_log10(  
    labels = scientific_10, 
    breaks = c(1e-5, 3e-5, 5e-5, 1e-4),  
    limits = c(5 * 10^-6, 7.1 * 10^-5)   
  ) +
  ylab("Measured recombination rate\n(per site per generation)") +  # Y-axis label
  theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.title.x = element_text(margin = margin(t = 20)),
    legend.text =  element_text(size = 20),
    legend.title =  element_text(size = 20),
    legend.position = "right"
  ) +
  xlab("Maximum time-scaled distance") +  # X-axis label
  scale_colour_manual(  
    values = cbPalette,  
    name = expression(paste("True recomb. rate (1e-5)")),
    labels = c("0.5", "0.75", "1", "2.5", "5")
  ) +
  scale_x_continuous(  
    breaks = c(50000, 75000, 100000, 150000),
    labels = c("50000", "75000", "100000", "150000")
  ) 
