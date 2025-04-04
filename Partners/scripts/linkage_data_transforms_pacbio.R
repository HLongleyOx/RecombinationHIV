##Script filters and processes dataset of linkage disequilibrium values outputted from linkageDisequilibrium_X.py
##Due to the size of the dataset, each individual is processed seperatley. 
##Script outputs a file with the D'ratios for every pairs of sites that are obsevred to have a D' > 0.2 at some time point 
##viral load data is also added  


library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
couple_w <- as.numeric(args[1])
type <- (args[2])
if (type=="donor"){
id = paste(couple_w, type)}else {id = paste(couple_w, "recip") }


#Read in necessary data for viral loads 
partners_features <- read_csv("partners_features.csv") %>%
  mutate(id_temp = sub("UW-", "", pt_id), couple=sub("\\d{2}$", "", id_temp)) %>%
  mutate(date=vl_dt, couple=as.numeric(couple), type=ifelse(status=="source","donor","recip"), ID=paste(couple, type))

if (type=="donor"){
  pacbio <- read_csv("linkage_donor_all.csv",
                       col_names = c("date", "siteA", "siteB", "ddash","baseA","baseB", "pA", "pB", "N", "ld","couple", "window","type")) %>%
    separate(window, into=c("window_start", "window_end") ,sep="_to_") %>% distinct() %>%
    mutate(Date=as.Date(date), window_start=as.numeric(window_start), window_end=as.numeric(window_end)) %>%
    group_by(siteA, siteB, window_start, couple, type) %>% mutate(days=as.numeric(Date-min(Date)), ld=abs(ld)) %>% filter(ld>0.0001) %>%
    arrange(days)  %>%
    mutate(siteA=siteA+window_start, siteB=siteB+window_start)%>% mutate(ID=paste(couple, type)) %>% mutate(ID=id) %>% filter(couple==couple_w)
} else{
  pacbio <- read_csv("linkage_recipient_all.csv",
                       col_names = c("date", "siteA", "siteB", "ddash","baseA","baseB", "pA", "pB", "N", "ld","couple", "window","type")) %>%
    separate(window, into=c("window_start", "window_end") ,sep="_to_") %>% distinct() %>%
    mutate(Date=as.Date(date), window_start=as.numeric(window_start), window_end=as.numeric(window_end)) %>%
    group_by(siteA, siteB, window_start, couple, type) %>% mutate(days=as.numeric(Date-min(Date)), ld=abs(ld)) %>% filter(ld>0.0001) %>%
    arrange(days)  %>%
    mutate(siteA=siteA+window_start, siteB=siteB+window_start) %>% mutate(ID=id) %>% filter(couple==couple_w)}


pacbio_dates <- pacbio %>% ungroup() %>% select(Date, ID) %>% distinct()
#if no corresponding date, find the closest date and record the VL result for that (loop through the couples) - add on info on whether this is imputed, and how far away is the imputed
ids_pacbio <- unique(pacbio$ID)
pacbio_vl <- data.frame()
#We need more than one time point
df <- pacbio_dates 
if (length(df$Date) < 2) {
  stop("Program terminated: 2 ")
}

features <- partners_features %>% 
  filter(ID == id) %>% 
  mutate(vl_dt = as.Date(vl_dt)) %>% 
  select(vl_dt, vl_result) %>% mutate(vl_result=ifelse(is.na(vl_result), 0, vl_result))

# Pre-allocate vectors
vl <- numeric(length(df$Date))
date_diffs <- numeric(length(df$Date))

if (max(features$vl_result)==0){
  for (i in seq_along(df$Date)) {
    vl[i] <- 0
    date_diffs[i] <- 0
  }
} else {
  
  for (i in seq_along(df$Date)) {
    dt <- df$Date[i]
    x <- which.min(abs(as.numeric(dt) - as.numeric(features$vl_dt)))
    date_diff <- min(abs(as.numeric(dt) - as.numeric(features$vl_dt)))
    vl[i] <- features$vl_result[x]
    date_diffs[i] <- min(abs(as.numeric(dt) - as.numeric(features$vl_dt)))
  }}
  
  if (is.infinite(abs(max(vl)))) {
    df$vl_result <- NA
    df$date_diffs <- NA
  } else {
    df$vl_result <- vl
    df$date_diffs <- date_diffs
  }
  
  # Now process this ID and write directly to file
  pacbio_vl <- df
  
  pacbio_vl$sample_id <- paste0(couple_w, "_", type, "_", pacbio_vl$Date)  
  
  pacbio <- pacbio %>% left_join(pacbio_vl)
  pacbio_first_timepoint <- pacbio %>% filter(ddash>0.2) %>% group_by(siteA, siteB, couple, window_start) %>% filter(days==min(days)) %>% mutate(first_day = days) %>%
    ungroup() %>% select(siteA, siteB, couple, window_start, first_day, type) %>% right_join(pacbio)
  
  pacbio_filtered <- pacbio_first_timepoint %>% filter(days>=first_day) %>%
    group_by(couple, type, siteA, siteB, window_start) %>%
    arrange(days) %>%
    mutate(last_days = days-lag(days), dist=abs(siteA-siteB), ddash_lag=lag(ddash), lagPa=lag(pA), lagPb=lag(pB), vl_lag=lag(vl_result), date_diffs_lag=lag(date_diffs)) %>%
    filter(last_days>50) %>%
    drop_na() %>% mutate(distt=(last_days/1.8)*dist, ddratio = -1*log(ddash/ddash_lag), viralload=(vl_result+ vl_lag)/2) %>% 
    dplyr:: mutate(pAdiff = abs(pA-lagPa), pBdiff=abs(pB-lagPb)) 
      
  pacbio_filtered$siteID <- paste0(pacbio_filtered$ID, "_", pacbio_filtered$siteA, "_", 
                                     pacbio_filtered$siteB, "_", pacbio_filtered$last_days)
  
  # Calculate frequency changes
  pacbio_filtered$freq_changeA1 <- pacbio_filtered$lagPa * (1 - pacbio_filtered$pA)
  pacbio_filtered$freq_changeA2 <- (1 - pacbio_filtered$lagPa) * (pacbio_filtered$pA)
  pacbio_filtered$freq_changeA <- abs(1/(pacbio_filtered$last_days/1.8) * 
                                          log(pacbio_filtered$freq_changeA1/pacbio_filtered$freq_changeA2))
  
  pacbio_filtered$freq_changeB1 <- pacbio_filtered$lagPb * (1 - pacbio_filtered$pB)
  pacbio_filtered$freq_changeB2 <- (1 - pacbio_filtered$lagPb) * (pacbio_filtered$pB)
  pacbio_filtered$freq_changeB <- abs(1/(pacbio_filtered$last_days/1.8) * 
                                          log(pacbio_filtered$freq_changeB1/pacbio_filtered$freq_changeB2))
   
  
  # Aggregate and select final columns
  bootstrap_pacbio <- pacbio_filtered %>% 
    group_by(siteID) %>% 
    mutate(ddratio = mean(ddratio), 
           freq_changeA = mean(freq_changeA),
           freq_changeB = mean(freq_changeB), 
           N = mean(N), 
           pA = mean(pA), 
           pB = mean(pB), 
           lagPa = mean(lagPa), 
           lagPb = mean(lagPb), 
           ddash_lag = mean(ddash_lag), 
           ddash = mean(ddash)) %>%
    select(siteA, siteB, days, ID, siteID, couple, type, ddratio, Date,
           freq_changeA, freq_changeB, N, pA, pB, last_days, distt, dist, 
           lagPa, lagPb, ddash_lag, ddash, sample_id, viralload, vl_lag, vl_result, date_diffs, date_diffs_lag) %>% 
    distinct()
  
  # Write this ID's data to file
  write_csv(bootstrap_pacbio, 
            file = paste0("modified_linkage_pacbio_final_",type,"_", couple_w,".csv"))
  
  
