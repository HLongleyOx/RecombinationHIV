########################################################################################
##This code take a dataset describing linkage over time for diverse sites in simualated sequence data
##Uses method from Romero and Feder 2024
########################################################################################


### Load libraries ####
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(minpack.lm)


#Load in simulation data
setwd("linkage") #Replace file path if necessary
recom = c("7.5e-06", "1e-05", "2.5e-05",  "5e-05")
df <- data.frame()

for (file in recom) {
  sel_linkage  <-   #Read in file and format variables if necessary
    read_csv(paste0("linkage_", file, ".csv"), show_col_types = FALSE) %>%
    mutate(
      days = as.numeric(sub("_([0-9]+)_", "\\1", date)),
      siteA = as.numeric(siteA),
      siteB = as.numeric(siteB),
      pA = as.numeric(pA),
      pB = as.numeric(pB),
      ddash = as.numeric(ddash),
      ld = as.numeric(ld)
    ) %>%
    drop_na() %>% group_by(ID, siteA, siteB, baseA, baseB) %>% filter(ld >=
                                                                        0.0001)
  
  
  sel_first_timepoint <-  #Only looks at data points when LD first exceeds 0.2 and all time points after 
    sel_linkage %>% filter(ddash > 0.2) %>% group_by(ID, siteA, siteB, baseA, baseB) %>% filter(days ==
                                                                                                  min(days)) %>% mutate(first_day = days) %>%
    ungroup() %>% select(siteA, siteB, first_day, ID) %>% right_join(sel_linkage) %>% filter(days >=
                                                                                               first_day)  %>%
    group_by(ID, siteA, siteB, baseA, baseB) %>%
    arrange(days) %>%
    mutate(
      last_days = days - lag(days),
      dist = abs(siteA - siteB),
      ddash_lag = lag(ddash),
      lagPa = lag(pA),
      lagPb = lag(pB),
      distt = (last_days) * dist,
      ddratio = -1 * log(ddash / ddash_lag),  #Calculate D ratio
      first_ld = first(ddash)
    ) %>%
    drop_na()
  
  

  for (d in seq(20000, 150000, by = 10000)) {    #Fits the recombination model to data filtered for different 
    sel <- sel_first_timepoint %>% drop_na() %>% filter(distt < d) %>%
      mutate(grp = round(distt /
                           1) * 1)
    
    x = sel %>% group_by(grp) %>% summarise(mn = mean(ddratio), N = n())
    a = x$mn
    b = x$grp
    n = sqrt(x$N)
    
    fit <-
      nlsLM(
        a ~ f(b, c0, c1, c2),
        start = list(c0 = 0.2, c1 = 0.3, c2 = 0.0001),
        weights = n
      )
    p = ((coef(fit)[[2]] * coef(fit)[[3]]))
    print(p)
    df = rbind(df, data.frame(
      r = p,
      dist = d,
      true_recom = file
    ))
  }
  
}

#Create plot of results
sim_plot <-
  ggplot(
    df %>% filter(true_recom %in% c("1e-05", "2.5e-05", "5e-05", "7.5e-06")) ,
    aes(
      x = dist,
      y = r,
      group = as.factor(true_recom),
      colour = as.factor(true_recom)
    )
  ) + geom_line() +
  theme_classic() + scale_y_log10(
    labels = scientific_10,
    breaks = c(1e-5, 3e-5, 5e-5, 1e-4),
    limits = c(5 * 10 ^ -6, 7.1 * 10 ^ -5)
  ) +
  ylab("Measured recombination rate\n(per site per generation)") + theme(
    axis.text = element_text(size = 20),
    axis.title = element_text(size =
                                20),
    axis.title.x = element_text(margin = margin(t = 20)),
    legend.text =  element_text(size =
                                  20),
    legend.title =  element_text(size =
                                   20),
    legend.position = "right"
  ) +
  xlab("Maximum time-scaled distance") +
  scale_colour_manual(
    values = cbPalette,
    name = expression(paste("True recomb. rate (1e-5)")),
    labels = c("0.5", "0.75", "1", "2.5", "5")
  ) +
  scale_x_continuous(breaks = c(50000, 75000, 100000, 150000),
                     labels = c("50000", "75000", "100000", "150000")) +
  geom_hline(
    data = df_50,
    aes(yintercept = as.numeric(true_recom), colour = true_recom),
    linetype = 2,
    alpha = 0.7,
    size = 2
  ) + guides(colour = "none")



