## set of scripts for cluster to generate following analyses:
#recombination rate across genome for illumina and pacbio data 
#the effect of selection - using different frequency cut off (change per generation)
#stratifying by viral load

#script 1 prepared the dataset - identifies sites where the initial LD is high enough for us to measure decay
#outputs the modified illumina and pacbio datasets
#incorperates information on HIV subtype, viral load and sex 


args <- commandArgs(trailingOnly=TRUE)
couple_w <- as.numeric(args[1])
type <- (args[2])
#id = paste(couple_w, type)

if (type=="donor"){
id = paste(couple_w, type)}else {id = paste(couple_w, "recip") }

print(couple_w)
print(type)

library(tidyverse)

partners_features <- read_csv("/well/fraser/users/hqh588/pacbio_updated/metadata/partners_features.csv") %>%
  mutate(id_temp = sub("UW-", "", pt_id), couple=sub("\\d{2}$", "", id_temp)) %>%
  mutate(date=vl_dt, couple=as.numeric(couple), type=ifelse(status=="source","donor","recip"), ID=paste(couple, type))

sex_subtype <- partners_features %>% select(ID, sex, subtype_bestref) %>% distinct() %>% drop_na() %>%
  group_by(ID) %>% mutate(subtype=ifelse("A" %in% subtype_bestref | "A1" %in% subtype_bestref |"A2" %in% subtype_bestref | "A3" %in% subtype_bestref | "A6" %in% subtype_bestref , "A",
                                         ifelse("B" %in% subtype_bestref, "B",
                                                ifelse("C" %in% subtype_bestref, "C",
                                                       ifelse("D" %in% subtype_bestref, "D", "other"))))) %>%
  select(-subtype_bestref) %>% distinct() 

# 
if (type=="donor"){
illumina <- read_csv("/well/fraser/users/hqh588/illumina/output/linkage_recombination/linkage_donor_250.csv",
                     col_names = c("date", "siteA", "siteB", "ddash","baseA","baseB", "pA", "pB", "N", "ld","couple", "window","type")) %>%
  separate(window, into=c("window_start", "window_end") ,sep="_to_") %>% distinct() %>%
  mutate(Date=as.Date(date), window_start=as.numeric(window_start), window_end=as.numeric(window_end)) %>%
  group_by(siteA, siteB, window_start, couple, type) %>% mutate(days=as.numeric(Date-min(Date)), ld=abs(ld)) %>% filter(ld>0.0001) %>%
  arrange(days)  %>%
  mutate(siteA=siteA+window_start, siteB=siteB+window_start)%>% mutate(ID=paste(couple, type)) %>% mutate(ID=id) %>% filter(couple==couple_w)
} else{
illumina <- read_csv("/well/fraser/users/hqh588/illumina/output/linkage_recombination/linkage_recipient_250.csv",
                           col_names = c("date", "siteA", "siteB", "ddash","baseA","baseB", "pA", "pB", "N", "ld","couple", "window","type")) %>%
  separate(window, into=c("window_start", "window_end") ,sep="_to_") %>% distinct() %>%
  mutate(Date=as.Date(date), window_start=as.numeric(window_start), window_end=as.numeric(window_end)) %>%
  group_by(siteA, siteB, window_start, couple, type) %>% mutate(days=as.numeric(Date-min(Date)), ld=abs(ld)) %>% filter(ld>0.0001) %>%
  arrange(days)  %>%
  mutate(siteA=siteA+window_start, siteB=siteB+window_start) %>% mutate(ID=id) %>% filter(couple==couple_w)}


illumina_dates <- illumina %>% ungroup() %>% select(Date, ID) %>% distinct()
#if no corresponding date, find the closest date and record the VL result for that (loop through the couples) - add on info on whether this is imputed, and how far away is the imputed
ids_illumina <- unique(illumina$ID)
illumina_vl <- data.frame()
#Dataset becomes too large to 
# Process one ID at a time
  df <- illumina_dates


if (length(df$Date) < 2) {
    stop("Program terminated: only 2 time points ")
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
  illumina_vl <- rbind(illumina_vl, df)
  
  # Now process this ID and write directly to file
  illumina_vl <- df
  
  illumina_vl$sample_id <- paste0(couple_w, "_", type, "_", illumina_vl$Date)  
  
  illumina <- illumina %>% left_join(illumina_vl)
  illumina_first_timepoint <- illumina %>% filter(ddash>0.2) %>% group_by(siteA, siteB, couple, window_start) %>% filter(days==min(days)) %>% mutate(first_day = days) %>%
    ungroup() %>% select(siteA, siteB, couple, window_start, first_day, type) %>% right_join(illumina)
  
  illumina_filtered <- illumina_first_timepoint %>% filter(days>=first_day) %>%
    group_by(couple, type, siteA, siteB, window_start) %>%
    arrange(days) %>%
    mutate(last_days = days-lag(days), dist=abs(siteA-siteB), ddash_lag=lag(ddash), lagPa=lag(pA), lagPb=lag(pB), vl_lag=lag(vl_result), date_diffs_lag=lag(date_diffs)) %>%
    filter(last_days>50) %>%
    drop_na() %>% mutate(distt=(last_days/1.8)*dist, ddratio = -1*log(ddash/ddash_lag), viralload=(vl_result+ vl_lag)/2) %>% 
    dplyr:: mutate(pAdiff = abs(pA-lagPa), pBdiff=abs(pB-lagPb)) %>%
    left_join(sex_subtype %>% filter(ID==id))
  
  # Continue with the rest of your calculations
  illumina_filtered$siteID <- paste0(illumina_filtered$ID, "_", illumina_filtered$siteA, "_", 
                                   illumina_filtered$siteB, "_", illumina_filtered$last_days)
  
  # Calculate frequency changes
  illumina_filtered$freq_changeA1 <- illumina_filtered$lagPa * (1 - illumina_filtered$pA)
  illumina_filtered$freq_changeA2 <- (1 - illumina_filtered$lagPa) * (illumina_filtered$pA)
  illumina_filtered$freq_changeA <- abs(1/(illumina_filtered$last_days/1.8) * 
                                        log(illumina_filtered$freq_changeA1/illumina_filtered$freq_changeA2))
  
  illumina_filtered$freq_changeB1 <- illumina_filtered$lagPb * (1 - illumina_filtered$pB)
  illumina_filtered$freq_changeB2 <- (1 - illumina_filtered$lagPb) * (illumina_filtered$pB)
  illumina_filtered$freq_changeB <- abs(1/(illumina_filtered$last_days/1.8) * 
                                        log(illumina_filtered$freq_changeB1/illumina_filtered$freq_changeB2))
  
  #If vl is taken on the same day for both observations, then record mean viral load, otherwise fill with NA: 
  #  illumina_filtered  <- illumina_filtered %>% ungroup() %>% mutate(viralload=ifelse(date_diffs==0 & lag_date_diff==0 & vl_result>0 & lag_vl>0, (vl_result+lag_vl)/2, NA))
  
  
  # Aggregate and select final columns
  bootstrap_illumina <- illumina_filtered %>% 
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
  write_csv(bootstrap_illumina, 
            file = paste0("/well/fraser/users/hqh588/illumina/output/linkage_recombination/modified_linkage_illumina_final_",type,"_", couple_w,".csv"))
  
  
  