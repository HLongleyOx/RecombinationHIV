library(tidyverse)

source("/well/fraser/users/hqh588/pacbio_updated/scripts/Rscripts/recombination_function.R")


args <- commandArgs(trailingOnly=TRUE)
w <- as.numeric(args[1])
print(w)
v1 <- seq(6615, 6695)
v2 <- seq(6693, 6812)
v3 <- seq(7110, 7217)
v4 <- seq(7377, 7478)
v5 <- seq(7602, 7634)

vloops <- c(v1, v2, v3, v4, v5)

bootstrap_illumina_filtered <- read_csv("/well/fraser/users/hqh588/illumina/output/linkage_recombination/modified_linkage_illumina_final_v.csv") %>%
 filter(!siteA %in% vloops, !siteB %in% vloops)

bootstrap_pacbio_filtered <- read_csv("/well/fraser/users/hqh588/pacbio_updated/output/linkage_recombination/modified_linkage_pacbio_final_v.csv") %>%
mutate(across(
    .cols = -c(siteID, type, ID),  # Exclude columns that should remain non-numeric
    .fns = ~ as.numeric(.),  # Convert all other columns to numeric
    .names = "{.col}"  # Retain original column names
  )) %>% filter(!siteA %in% vloops, !siteB %in% vloops)

  
#Only include pairs where pA, pB, lagpA, lagpB are all greater than 0.05, or N>2 (whichever is larger) 
bootstrap_pacbio_filtered <- bootstrap_pacbio_filtered %>% filter(pA>0.05, pB>0.05, lagPa>0.05, lagPb>0.05)
bootstrap_illumina_filtered <- bootstrap_illumina_filtered %>% filter(pA>0.05, pB>0.05, lagPa>0.05, lagPb>0.05)
dualtypes <- read_csv("/exafs1/well/fraser/users/hqh588/pacbio_updated/metadata/dualtypes.csv")


sex_subtype = read_csv("/well/fraser/users/hqh588/pacbio_updated/output/sex_subtype.csv")
subtypeA <- sex_subtype %>% filter(subtype=="A")
subtypeC <- sex_subtype %>% filter(subtype=="C")
subtypeD <- sex_subtype %>% filter(subtype=="D")

print(head(bootstrap_pacbio_filtered ))
illumina_recom <- recombination_estimate(bootstrap_illumina_filtered %>% filter(siteA>=w, siteB>=w, siteA<(w+750),
                                                                              siteB<(w+750)),
                                       viralload_max = 10^8,
                                       viralload_min = 1000,
                                       sel=1,
                                       timeframe_min = 50,
                                       timeframe_max = 5000,
                                       subtype_choice = "All",
                                       start=7,
                                       end=10000, max_distt = 50000)



pacbio_recom <-  recombination_estimate(bootstrap_pacbio_filtered %>% filter(siteA>=w, siteB>=w, siteA<(w+500),
                                                                                siteB<(w+500)),
viralload_max = 10^8,
viralload_min = 1000,
sel=1,
timeframe_min = 50,
timeframe_max = 5000,
subtype_choice = "All",
start=7,
end=10000, max_distt = 75000)

pacbio_recom750 <-  recombination_estimate(bootstrap_pacbio_filtered %>% filter(siteA>=w, siteB>=w, siteA<(w+750),
                                                                               siteB<(w+750), 
),
viralload_max = 10^8,
viralload_min = 1000,
sel=1,
timeframe_min = 50,
timeframe_max = 5000,
subtype_choice = "All",
start=7,
end=10000, max_distt = 75000)

#subtype


pacbio_recomA <-  recombination_estimate(bootstrap_pacbio_filtered %>% filter(siteA>=w, siteB>=w, siteA<(w+750),
                                                                             siteB<(w+750), !ID %in% dualtypes$ID, ID %in% subtypeA$ID 
),
viralload_max = 10^8,
viralload_min = 1000,
sel=1,
timeframe_min = 50,
timeframe_max = 5000,
subtype_choice = "All",
start=7,
end=10000, max_distt = 75000)


pacbio_recomC <-  recombination_estimate(bootstrap_pacbio_filtered %>% filter(siteA>=w, siteB>=w, siteA<(w+750),
                                                                             siteB<(w+750), !ID %in% dualtypes$ID , ID %in% subtypeC$ID 
),
viralload_max = 10^8,
viralload_min = 1000,
sel=1,
timeframe_min = 50,
timeframe_max = 5000,
subtype_choice = "All",
start=7,
end=10000, max_distt = 75000)


pacbio_recomD <-  recombination_estimate(bootstrap_pacbio_filtered %>% filter(siteA>=w, siteB>=w, siteA<(w+750),
                                                                             siteB<(w+500), !ID %in% dualtypes$ID, ID %in% subtypeD$ID 
),
viralload_max = 10^8,
viralload_min = 1000,
sel=1,
timeframe_min = 50,
timeframe_max = 5000,
subtype_choice = "All",
start=7,
end=10000, max_distt = 75000)


write_csv(pacbio_recom, file=paste0("/well/fraser/users/hqh588/pacbio_updated/output/recombination_output/final_output/pacbio_windows_",w, ".csv"))
write_csv(pacbio_recomA, file=paste0("/well/fraser/users/hqh588/pacbio_updated/output/recombination_output/final_output/pacbio_windowsA_",w, ".csv"))
write_csv(pacbio_recomD, file=paste0("/well/fraser/users/hqh588/pacbio_updated/output/recombination_output/final_output/pacbio_windowsD_",w, ".csv"))
write_csv(pacbio_recomC, file=paste0("/well/fraser/users/hqh588/pacbio_updated/output/recombination_output/final_output/pacbio_windowsC_",w, ".csv"))
write_csv(illumina_recom, file=paste0("/well/fraser/users/hqh588/pacbio_updated/output/recombination_output/final_output/illumina_windows_",w, ".csv"))
write_csv(pacbio_recom750, file=paste0("/well/fraser/users/hqh588/pacbio_updated/output/recombination_output/final_output/pacbio_long_windows_",w, ".csv"))










