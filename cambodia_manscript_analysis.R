########################################################
#_______________________________________________________
# @Organization --  UCSF 
# @Project -- Cambodia multipathogen burden
# @Author -- Samantha Bents, sjbents@stanford.edu, November 30, 2024 
# @Description -- Assess integrated intervention using rao diversity 
#_______________________________________________________
########################################################
library(here)
source(here("0-config.R"))
source(here("0-functions.R"))


#######################################################
# Load disease data. 
cambodia_serology_public = readr::read_csv(file = here("projects/6-multipathogen-burden/data/cambodia", "cambodia_serology_public.csv"))
#ttmb     : Tetanous toxoid
#bm14     : Lymphatic filariasis bm14
#bm33     : Lymphatic filariasis bm33
#wb123    : Lymphatic filariasis wb123
#nie      : Strongyloides stercoralis NIE
#sag2a    : Toxoplasma gondii SAG2A
#t24      : Taenia solium T24
#pfmsp19  : Plasmodium falciparum MSP-1(19)
#pvmsp19  : Plasmodium vivax MSP-1(19)
# Load location data 
gps_cambodia =   readr::read_csv(file = here("projects/6-multipathogen-burden/data/cambodia", "cambodia_ea_dhs.csv")) 

load(file = here("projects/6-multipathogen-burden/data/cambodia/cambodia_serology.Rdata"))
gps_dat = cambodia_serology %>%
  dplyr::select(dhsclust, psuid) %>% distinct()

#######################################################
# Transform into seropositive vs seronegative using: https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0004699#pntd-0004699-t001
# https://journals.asm.org/doi/10.1128/cvi.00052-16#T1 for tetanus
cambodia_seropositivity = cambodia_serology_public %>%
  mutate(tetanus_unprotected = ifelse(ttmb >= 100, 0,1 )) %>% # indicates were vaccinated I think 
  mutate(lymph14 = ifelse(bm14 >= 65, 1, 0 )) %>%
  mutate(lymph123 = ifelse(wb123 >= 115, 1, 0 )) %>%
  mutate(lymph33 = ifelse(bm33 >= 966, 1, 0 )) %>%
  mutate(strong = ifelse(nie  >= 792, 1, 0 )) %>%
  mutate(toxo = ifelse(sag2a  >= 159, 1, 0 )) %>%
  mutate(cyst = ifelse(t24  >= 486, 1, 0 )) %>%
  mutate(malaria_f = ifelse(pfmsp19  >= 343, 1, 0 )) %>%
  mutate(malaria_v = ifelse(pvmsp19  >= 196, 1, 0 )) %>%
  dplyr::select(region, psuid, age, tetanus_unprotected, lymph14, lymph123, 
                lymph33, strong, toxo, cyst, malaria_f, malaria_v) %>%
  mutate(malaria_any = ifelse(malaria_f == 1 | malaria_v == 1, 1, 0)) %>%
  mutate(lymph_any = ifelse(lymph14 == 1 | lymph123 == 1| lymph33 == 1, 1, 0)) %>%
  mutate(psuid = psuid + 1)

head(cambodia_seropositivity)
round(colSums(cambodia_seropositivity[4:14])/length(cambodia_seropositivity), digits = 2) 

psuid_n = length(unique(cambodia_seropositivity$psuid))

#######################################################
# Method 1: Population-level denominator and assess each STH independently. 
exlcude = c("lymph14", "lymph123", "lymph33", "malaria_any", "tetanus_unprotected", "cyst", "toxo")

# Number of non-NA measurements assessed total by pathogen 
parasite_denom = cambodia_seropositivity  %>%
  pivot_longer(cols=c('tetanus_unprotected', 'lymph14', 'lymph123', 'lymph33', 'strong', 'toxo', 'cyst',
                      'malaria_f', 'malaria_v', 'malaria_any',  'lymph_any'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(pathogen) %>% # assess by overall prevalence in population
  summarise(non_na_count = sum(!is.na(presence))) %>%
  filter(!pathogen %in% exlcude)

# Number of parasite instances by block
parasite_numer = cambodia_seropositivity %>%
  pivot_longer(cols=c('tetanus_unprotected', 'lymph14', 'lymph123', 'lymph33', 'strong', 'toxo', 'cyst',
                      'malaria_f', 'malaria_v', 'malaria_any',  'lymph_any'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(psuid, pathogen) %>%
  summarize(across(presence, ~sum(.x, na.rm = TRUE)))  %>%
  filter(!pathogen %in% exlcude)

# Join together and calculate fraction prevalence and fraction if disease was distributed evenly amongst blocks. 
parasite_prevalence_block = left_join(parasite_denom, parasite_numer, by = c("pathogen")) %>%
  mutate(fraction = presence/non_na_count ) %>% # fraction = proportion of total pathogen prevalence accounted for by each block
  ungroup() %>%
  group_by(pathogen) %>%
  mutate(pop_presence = sum(presence), even_prevalence = (pop_presence/non_na_count)/psuid_n) %>%
  dplyr::select(-non_na_count) %>% ungroup()
head(parasite_prevalence_block)


#######################################################
# Calculate Rao diversity at each block. 
product_results = list()


for (i in 1:psuid_n) {
  
  print(i)
  # Subset the data for the current block
  block_subset <- parasite_prevalence_block %>%
    mutate(fraction_mirrored = fraction) %>%
    ungroup() %>%
    filter(psuid == i) %>%
    dplyr::select(fraction_mirrored, fraction)
  
  # Get unique combinations of pathogen prevalences 
  fractions <- block_subset$fraction
  unique_combinations <- t(combn(fractions, 2, simplify = TRUE))  # Get pairs (no duplicates)
  
  # Compute products for each pair
  products <- apply(unique_combinations, 1, prod)
  
  products_df <- data.frame(
    psuid = i, 
    fraction_mirrored = unique_combinations[, 1],
    fraction = unique_combinations[, 2],
    product = products
  )
  # Store the results in the list
  print(products_df)
  product_results[[i]] <- products_df
  
}


# Combine all results into a single data frame
product_results <- do.call(rbind, product_results)

# Sum by block
block_prevalence_final = product_results %>%
  group_by(psuid) %>%
  summarize(across(product, ~sum(.x)))
head(block_prevalence_final)

# Calculate expected rao diversity for a given block.
block_subset = parasite_prevalence_block %>% 
  mutate(even_prevalence_mirrored = even_prevalence) %>%
  ungroup() %>% filter(psuid == 1) %>%
  dplyr::select(even_prevalence_mirrored, even_prevalence)

even_prevalences = block_subset$even_prevalence
unique_combinations <- t(combn(even_prevalences, 2, simplify = TRUE))  # Get pairs (no duplicates)

# Compute products for each pair
products <- apply(unique_combinations, 1, prod)
null_rao = sum(products)

# Compare observed rao to expected rao and add gps data. 
method1 = block_prevalence_final %>%
  mutate(rao = product/null_rao) %>%
  ungroup()
summary(method1$rao)  


#######################################################
#######################################################
# Method 3: Pathogens-specific denominator and assess each STH independently. 

# Number of non-NA measurements assessed total by pathogen 
parasite_denom = cambodia_seropositivity  %>%
  pivot_longer(cols=c('tetanus_unprotected', 'lymph14', 'lymph123', 'lymph33', 'strong', 'toxo', 'cyst',
                      'malaria_f', 'malaria_v', 'malaria_any',  'lymph_any'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(pathogen) %>% # assess by overall prevalence in population
  summarise(non_na_count = sum(presence, na.rm = TRUE))  %>%
  filter(!pathogen %in% exlcude)


# Number of total parasite instances by block
# Number of parasite instances by block
parasite_numer = cambodia_seropositivity %>%
  pivot_longer(cols=c('tetanus_unprotected', 'lymph14', 'lymph123', 'lymph33', 'strong', 'toxo', 'cyst',
                      'malaria_f', 'malaria_v', 'malaria_any',  'lymph_any'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(psuid, pathogen) %>%
  summarize(across(presence, ~sum(.x, na.rm = TRUE)))  %>%
  filter(!pathogen %in% exlcude)


# Join together and calculate fraction unvaccinated and fraction unvaccinated if disease was distributed evenly amongst blocks. 
parasite_prevalence_block = left_join(parasite_denom, parasite_numer, by = c("pathogen")) %>%
  mutate(fraction = presence/non_na_count ) %>% # fraction = proportion of total prevalence accounted for by each block
  ungroup() %>%
  group_by(pathogen) %>%
  mutate(pop_presence = sum(presence), even_prevalence = 1/psuid_n) %>%
  dplyr::select(-non_na_count) 

#######################################################
# Calculate Rao diversity at each block. 
product_results = list()

for (i in 1:psuid_n) {
  
  # Subset the data for the current block
  block_subset <- parasite_prevalence_block %>%
    mutate(fraction_mirrored = fraction) %>%
    ungroup() %>%
    filter(psuid == i) %>%
    dplyr::select(fraction_mirrored, fraction)
  
  # Get unique combinations of pathogen prevalences 
  fractions <- block_subset$fraction
  unique_combinations <- t(combn(fractions, 2, simplify = TRUE))  # Get pairs (no duplicates)
  
  # Compute products for each pair
  products <- apply(unique_combinations, 1, prod)
  
  products_df <- data.frame(
    psuid = i, 
    fraction_mirrored = unique_combinations[, 1],
    fraction = unique_combinations[, 2],
    product = products
  )
  
  # Store the results in the list
  product_results[[i]] <- products_df
}


# Combine all results into a single data frame
product_results <- do.call(rbind, product_results)

# Sum by block
block_prevalence_final = product_results %>%
  group_by(psuid) %>%
  summarize(across(product, ~sum(.x)))
head(block_prevalence_final)

# Calculate expected rao diversity for a given block.
block_subset = parasite_prevalence_block  %>% 
  mutate(even_prevalence_mirrored = even_prevalence) %>%
  ungroup() %>% filter(psuid == 1) %>%
  dplyr::select(even_prevalence_mirrored, even_prevalence)

even_prevalences = block_subset$even_prevalence
unique_combinations <- t(combn(even_prevalences, 2, simplify = TRUE))  # Get pairs (no duplicates)

# Compute products for each pair
products <- apply(unique_combinations, 1, prod)
null_rao = sum(products)

# Compare observed rao to expected rao and add gps data. 
method3 = block_prevalence_final %>%
  mutate(rao = product/null_rao) %>%
  ungroup()
summary(method3$rao)



#######################################################
# Compare method efficiency using # of blocks 
pt = .5 # targeting threshold for disease

block_50 <- list()
efficiency_list <- list()
methods = list(method1 = method1, method3 = method3)

for(i in seq_along(methods)){
  
  method_df = methods[[i]]
  
  compare_strategy = left_join(parasite_prevalence_block, method_df, by = "psuid")
  
  # Malaria f strategy 
  # motivated by itself
  malf_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_f") %>%
    arrange(-fraction) %>% # Arrange from highest prevalence block to lowest 
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  order_malf = print(malf_motivated$psuid)
  
  malf_50 <- which(malf_motivated$above == 2)[1]
  
  # rao motivated 
  rao_malf_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_f") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  # print(rao_malf_motivated$psuid)
  
  # print(head(rao_malf_motivated))
  
  rao_malf_50 <- which(rao_malf_motivated$above == 2)[1]
  
  #Other pathogens - motivated by malv v
  order_malv = compare_strategy %>%
    filter(pathogen == "malaria_v") %>%
    arrange(-fraction) 
  order_malv = print(order_malv$psuid)
  
  malf_malv_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_malv) %>%
    filter(pathogen == "malaria_f") %>%
    arrange(order_malv) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  
  malf_malv_50 <- which(malf_malv_motivated$above == 2)[1]
  
  # Other pathogens = lf
  order_lf = compare_strategy %>%
    filter(pathogen == "lymph_any") %>%
    arrange(-fraction) 
  order_lf = print(order_lf$psuid)
  
  malf_lf_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_lf) %>%
    filter(pathogen == "malaria_f") %>%
    arrange(order_lf) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  
  malf_lf_50 <- which(malf_lf_motivated$above == 2)[1]
  
  malf = ggplot(data = malf_motivated, aes(x = block_label/psuid_n*100, y = target*100)) +
    geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "Malaria falciparum"), cex = 2) +
    geom_line(data = rao_malf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = malf_malv_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Malaria vivax"), cex = 2) + 
    geom_line(data = malf_lf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "LF"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("Malaria falciparum" = "black", "Rao" = "#FF4F00",
                                  "Malaria vivax" = "#F0A20E", "LF" = "#0094C6")) +
    guides(color =guide_legend(title="Strategy")) +
    ggtitle("      Malaria falciparum")+
    theme(axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "black", size = 12.5),
          legend.position = "right",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 11, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          strip.text = element_text(colour = "black", size = 10, hjust = 0),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill=NA))  +
    theme(legend.position = "none") 
  
  # Calculate AUC for malf 
  auc_malf <- trapz(malf_motivated$block_label/psuid_n, malf_motivated$target)
  auc_malf_rao <- trapz(rao_malf_motivated$block_label/psuid_n, rao_malf_motivated$target)
  auc_malf_malv <- trapz(malf_malv_motivated$block_label/psuid_n, malf_malv_motivated$target)
  auc_malf_lf <- trapz(malf_lf_motivated$block_label/psuid_n, malf_lf_motivated$target)
  malf_auc = list(auc_malf,auc_malf_rao, auc_malf_malv, auc_malf_lf )
  print(malf_auc)
  
  # strong strategy
  strong_motivated = compare_strategy %>% 
    filter(pathogen == "strong") %>%
    arrange(-fraction) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  strong_50 <- which(strong_motivated$above == 2)[1]
  
  rao_strong_motivated = compare_strategy %>% 
    filter(pathogen == "strong") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  rao_strong_50 <- which(rao_strong_motivated$above == 2)[1]
  
  #Other pathogens - motivated by malv v
  strong_malv_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_malv) %>%
    filter(pathogen == "strong") %>%
    arrange(order_malv) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  strong_malv_50 <- which(strong_malv_motivated$above == 2)[1]
  
  # Other pathogens = lf
  strong_lf_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_lf) %>%
    filter(pathogen == "strong") %>%
    ungroup() %>%
    arrange(order_lf) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  strong_lf_50 <- which(strong_lf_motivated$above == 2)[1]
  
  # Other pathogens = malf
  strong_malf_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_malf) %>%
    filter(pathogen == "strong") %>%
    arrange(order_malf) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  strong_malf_50 <- which(strong_malf_motivated$above == 2)[1]
  
  strong = ggplot(data = strong_motivated, aes(x = block_label/psuid_n*100, y = target*100)) +
    geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "Strongyloides"), cex = 2) +
    geom_line(data = rao_strong_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = strong_malv_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Malaria vivax"), cex = 2) + 
    geom_line(data = strong_lf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "LF"), cex = 2) + 
    geom_line(data = strong_malf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Malaria falciparum"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("Strongyloides" = "black", "Rao" = "#FF4F00", "Malaria falciparum" = "#004B6F" ,
                                  "Malaria vivax" = "#F0A20E", "LF" = "#0094C6")) +
    guides(color =guide_legend(title="Strategy")) +
    ggtitle("      Strongyloides")+
   # theme(legend.position = "right") +
    theme(axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "black", size = 12.5),
          legend.position = "right",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 11, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          strip.text = element_text(colour = "black", size = 10, hjust = 0),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill=NA)) + 
    theme(legend.position = "none") 
  
  # Calculate AUC for strongy  
  auc_strong <- trapz(strong_motivated$block_label/psuid_n, strong_motivated$target)
  auc_strong_rao <- trapz(rao_strong_motivated$block_label/psuid_n, rao_strong_motivated$target)
  auc_strong_malv <- trapz(strong_malv_motivated$block_label/psuid_n, strong_malv_motivated$target)
  auc_strong_malf <- trapz(strong_malf_motivated$block_label/psuid_n, strong_malf_motivated$target)
  auc_strong_lf <- trapz(strong_lf_motivated$block_label/psuid_n, strong_lf_motivated$target)
  strong_auc = list(auc_strong,  auc_strong_rao, auc_strong_malv, auc_strong_malf,  auc_strong_lf  )
  print(strong_auc)
  
  
  # LF
  lymph_motivated = compare_strategy %>% 
    filter(pathogen == "lymph_any") %>%
    arrange(-fraction) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  lymph_50 <- which(lymph_motivated$above == 2)[1]
  
  rao_lymph_motivated = compare_strategy %>% 
    filter(pathogen == "lymph_any") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  rao_lymph_50 <- which(rao_lymph_motivated$above == 2)[1]
  
  # Other pathogens = malf
  lf_malf_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_malf) %>%
    filter(pathogen == "lymph_any") %>%
    arrange(order_malf) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  lf_malf_50 <- which(lf_malf_motivated$above == 2)[1]
  
  # Other pathogens - mal v 
  lf_malv_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_malv) %>%
    filter(pathogen == "lymph_any") %>%
    arrange(order_malv) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  lf_malv_50 <- which(lf_malv_motivated$above == 2)[1]
  
  lymph = ggplot(data = lymph_motivated, aes(x = block_label/psuid_n*100, y = target*100)) +
    geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "LF"), cex = 2) +
    geom_line(data = rao_lymph_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = lf_malv_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Malaria vivax"), cex = 2) + 
    geom_line(data = lf_malf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Malaria falciparum"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("LF" = "black", "Rao" = "#FF4F00",
                                  "Malaria vivax" = "#F0A20E" , "Malaria falciparum" = "#004B6F" )) +
    guides(color =guide_legend(title="Strategy")) +
    ggtitle("      Lymphatic filariasis")+
    theme(legend.position = "bottom") +
    theme(axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "black", size = 12.5),
          legend.position = "right",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 11, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          strip.text = element_text(colour = "black", size = 10, hjust = 0),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill=NA))  +
    theme(legend.position = "none") 
  
  # Calculate AUC for malf 
  auc_lf <- trapz(lymph_motivated$block_label/psuid_n, lymph_motivated$target)
  auc_lf_rao <- trapz(rao_lymph_motivated$block_label/psuid_n, rao_lymph_motivated$target)
  auc_lf_malv <- trapz(lf_malv_motivated$block_label/psuid_n,lf_malv_motivated$target)
  auc_lf_malf <- trapz(lf_malf_motivated$block_label/psuid_n, lf_malf_motivated$target)
  lf_auc = list(auc_lf, auc_lf_rao, auc_lf_malv, auc_lf_malf )
  print(lf_auc)
  
  
  # malaria vivax 
  malv_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_v") %>%
    arrange(-fraction) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  malv_50 <- which(malv_motivated$above == 2)[1]
  
  rao_malv_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_v") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  rao_malv_50 <- which(rao_malv_motivated$above == 2)[1]
  
  malv_lf_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_lf) %>%
    filter(pathogen == "malaria_v") %>%
    arrange(order_lf) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  malv_lf_50 <- which(malv_lf_motivated$above == 2)[1]
  
  # Other pathogens - mal v 
  malv_malf_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_malf) %>%
    filter(pathogen == "malaria_v") %>%
    arrange(order_malf) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  malv_malf_50 <- which(malv_malf_motivated$above == 2)[1]
  
  malv = ggplot(data = malv_motivated, aes(x = block_label/psuid_n*100, y = target*100)) +
    geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "Malaria vivax"), cex = 2) +
    geom_line(data = rao_malv_motivated , aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = malv_lf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "LF"), cex = 2) + 
    geom_line(data = malv_malf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Malaria falciparum"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("Malaria vivax" = "black", "Rao" = "#FF4F00", "Malaria falciparum" ="#004B6F" ,
                               "LF" = "#0094C6")) +
  #  scale_color_viridis_d(option = "H", direction = 1, end = 1, begin = .6) +
    guides(color =guide_legend(title="Strategy")) +
    ggtitle("      Malaria vivax")+
    theme(legend.position = "right") +
    theme(axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "black", size = 12.5),
          legend.position = "right",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 11, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          strip.text = element_text(colour = "black", size = 10, hjust = 0),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill=NA)) +
    theme(legend.position = "none") 
  
  # Calculate AUC for malv
  auc_malv <- trapz(malv_motivated$block_label/psuid_n, malv_motivated$target)
  auc_malv_rao <- trapz(rao_malv_motivated$block_label/psuid_n, rao_malv_motivated$target)
  auc_malv_malf <- trapz(malv_malf_motivated$block_label/psuid_n, malv_malf_motivated$target)
  auc_malv_lf <- trapz(malv_lf_motivated$block_label/psuid_n, malv_lf_motivated$target)
  malv_auc = list( auc_malv, auc_malv_rao,  auc_malv_malf, auc_malv_lf )
  print(malv_auc)
  
  block_50[[i]] <- list(malf_50 = malf_50, rao_malf_50 = rao_malf_50, malf_malv_50 = malf_malv_50,
                        malf_lf_50 = malf_lf_50, strong_50 = strong_50, rao_strong_50 = rao_strong_50, 
                        strong_malv_50 = strong_malv_50,  strong_lf_50 = strong_lf_50, strong_malf_50 = 
                          strong_malf_50, lymph_50 = lymph_50, rao_lymph_50 = rao_lymph_50,
                        lf_malf_50  =lf_malf_50,lf_malv_50 = lf_malv_50, malv_50 = malv_50, rao_malv_50 = rao_malv_50,
                        malv_lf_50 = malv_lf_50 , malv_malf_50 = malv_malf_50)
  
  
  four_path_plot = plot_grid(malf, malv, strong, lymph)
  efficiency_list[[i]] = four_path_plot 
}


print(block_50[[1]])
print(block_50[[2]])

# Plot efficiency of single-pathogen strategies vs rao strategies. 
method1_efficiency <- plot_grid(efficiency_list = efficiency_list[[1]], ncol = 2)
method1_efficiency 
method2_efficiency <- plot_grid(efficiency_list = efficiency_list[[2]], ncol = 2)
method2_efficiency


# make legend for strategy 
legend1 = ggplot(data = parasite_prevalence_block) +
  geom_line(aes(x = psuid, y = presence, col = "Rao-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = fraction, col = "Malaria falciparum-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = even_prevalence, col = "Malaria vivax-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = pop_presence, col = "LF-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = pop_presence* even_prevalence, col = "Single-pathogen motivated"), cex = 2) + 
  scale_color_manual(values = c("Malaria vivax-motivated" = "#F0A20E" , "Rao-motivated" = "#FF4F00", "Malaria falciparum-motivated" = "#004B6F",
                                "LF-motivated" = "#0094C6", "Single-pathogen motivated" = "black")) +
  guides(color =guide_legend(title="Strategy"))  +
  theme(legend.position = "right") +
  geom_blank() 
legend1


leg = get_legend(legend1)
grid.newpage()
legend_only = grid.draw(leg)

plot_grid(method1_efficiency , legend_only, rel_widths = c(1, .2), labels = "a")
legend_only = grid.draw(leg)




###########################

# Figure 3a
# simulations (compared to method 1 )
sim_50 <- list()

for(i in 1:1000) {
  
  method_df = method1
  compare_strategy = left_join(parasite_prevalence_block, method_df, by = "psuid")
  
  # rao motivated 
  rao_malf_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_f") %>%
    mutate(random = sample(1:100, 100, replace = FALSE)) %>%
    arrange(-random) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  
  rao_malf_50_sim <- which(rao_malf_motivated$above == 2)[1]
  sim_50[[i]] <- rao_malf_50_sim
  
}

malaria_f_vector <- unlist(sim_50)
malaria_f_sim <- data.frame(malaria_f_vector)


# simulations
sim_50 <- list()

for(i in 1:1000) {
  
  method_df = method1
  compare_strategy = left_join(parasite_prevalence_block, method_df, by = "psuid")
  
  # rao motivated 
  rao_malv_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_v") %>%
    mutate(random = sample(1:100, 100, replace = FALSE)) %>%
    arrange(-random) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  
  rao_malv_50_sim <- which( rao_malv_motivated $above == 2)[1]
  sim_50[[i]] <- rao_malv_50_sim
  
}

malaria_v_vector <- unlist(sim_50)
malaria_v_sim <- data.frame(malaria_v_vector)


sim_50 <- list()

for(i in 1:1000) {
  
  method_df = method1
  compare_strategy = left_join(parasite_prevalence_block, method_df, by = "psuid")
  
  # rao motivated 
  rao_strong_motivated = compare_strategy %>% 
    filter(pathogen == "strong") %>%
    mutate(random = sample(1:100, 100, replace = FALSE)) %>%
    arrange(-random) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  
  rao_strong_50_sim <- which( rao_strong_motivated $above == 2)[1]
  sim_50[[i]] <- rao_strong_50_sim
  
}


strong_vector <- unlist(sim_50)
strong_sim <- data.frame(strong_vector)


# LF 
print(block_50[[1]])
sim_50 <- list()

for(i in 1:1000) {
  
  method_df = method1
  compare_strategy = left_join(parasite_prevalence_block, method_df, by = "psuid")
  
  # rao motivated 
  rao_lf_motivated = compare_strategy %>% 
    filter(pathogen == "lymph_any") %>%
    mutate(random = sample(1:100, 100, replace = FALSE)) %>%
    arrange(-random) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(cum < pt, 1, 2))
  
  rao_lf_50_sim <- which( rao_lf_motivated $above == 2)[1]
  sim_50[[i]] <- rao_lf_50_sim
  
}

lf_vector <- unlist(sim_50)
lf_sim <- data.frame(lf_vector)

############################## make a boxplot with them all

head(malaria_f_sim)
mf_boxplot = malaria_f_sim %>%
  mutate(pathogen = "Malaria falciparum") %>%
  mutate(Rao = 8) %>%
  mutate(Optimal = 6) %>%
  mutate(P.vivax = 67) %>%
  mutate(LF = 49)
head(mf_boxplot)

mv_boxplot = malaria_v_sim %>%
  mutate(pathogen = "Malaria vivax") %>%
  mutate(Rao = 14) %>%
  mutate(Optimal = 11) %>%
  mutate(P.falciparum = 49) %>%
  mutate(LF = 47)
head(mv_boxplot)

strong_boxplot = strong_sim %>%
  mutate(pathogen = "Strongyloides") %>%
  mutate(Rao = 39) %>%
  mutate(Optimal = 32) %>%
  mutate(P.falciparum = 54) %>%
  mutate(P.vivax = 50) %>%
  mutate(LF = 50)
head(strong_boxplot)

lf_boxplot = lf_sim %>%
  mutate(pathogen = "LF") %>%
  mutate(Rao = 34) %>%
  mutate(Optimal = 31) %>%
  mutate(P.falciparum = 50) %>%
  mutate(P.vivax = 52) 


target_goal = ggplot(data = mv_boxplot) + 
  geom_violin(aes(x = pathogen, y = malaria_v_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
 # geom_boxplot(aes(x = pathogen, y = malaria_v_vector), fill = "gray95", width=.6) + 
  geom_point(aes(x = pathogen, y = Rao),  col = "#FF4F00", cex = 6)  + #.5, rao  
  geom_point(aes(x = pathogen, y = Optimal), col = "black",  cex = 6)  + 
  geom_point(aes(x = pathogen, y = P.falciparum), col = "#004B6F",  cex = 6)  + 
  geom_point(aes(x = pathogen, y = LF), col = "#0094C6",  cex = 6)  +
  theme_bw() +
  ylab("Clusters to target 50% of disease") +
  geom_violin(data = mf_boxplot, aes(x = pathogen, y = malaria_f_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
 # geom_boxplot(data = mf_boxplot, aes(x = pathogen, y = malaria_f_vector), fill = "gray95", width=.6) + 
  geom_point(data = mf_boxplot, aes(x = pathogen, y = Rao), col = "#FF4F00",  cex = 6)  + #.5, rao  
  geom_point(data = mf_boxplot, aes(x = pathogen, y = Optimal), col = "black",  cex = 6)  + 
  geom_point(data = mf_boxplot, aes(x = pathogen, y = P.vivax), col = "#F0A20E",  cex = 6)  + 
  geom_point(data = mf_boxplot, aes(x = pathogen, y = LF), col = "#0094C6",  cex = 6)  +
  geom_violin(data = strong_boxplot, aes(x = pathogen, y = strong_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
 # geom_boxplot(data = strong_boxplot, aes(x = pathogen, y = strong_vector), fill = "gray95", width=.6) + 
  geom_point(data = strong_boxplot, aes(x = pathogen, y = Rao, col = "Rao-motivated"),  cex = 6)  + #.5, rao  
  geom_point(data = strong_boxplot, aes(x = pathogen, y = Optimal, col = "Single-pathogen motivated"),  cex = 6)  + 
  geom_point(data = strong_boxplot, aes(x = pathogen, y = P.vivax, col = "Malaria vivax-motivated") , cex = 6)  + 
  geom_point(data = strong_boxplot, aes(x = pathogen, y = P.falciparum, col = "Malaria falciparum-motivated"),  cex = 6)  + 
  geom_point(data = strong_boxplot, aes(x = pathogen, y = LF, col = "LF-motivated"),   cex = 6)  +
  geom_violin(data = lf_boxplot, aes(x = pathogen, y = lf_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
 # geom_boxplot(data = lf_boxplot, aes(x = pathogen, y = lf_vector), fill = "gray95", width=.6) + 
  geom_point(data = lf_boxplot, aes(x = pathogen, y = Rao), col = "#FF4F00",  cex = 6)  + #.5, rao  
  geom_point(data = lf_boxplot, aes(x = pathogen, y = Optimal), col = "black", cex = 6)  + 
  geom_point(data = lf_boxplot, aes(x = pathogen, y = P.vivax), col = "#F0A20E",  cex = 6)  + 
  geom_point(data = lf_boxplot, aes(x = pathogen, y = P.falciparum), col = "#004B6F",  cex = 6)   +
  scale_color_manual(values = c("Rao-motivated" = "#FF4F00", "Single-pathogen motivated" = "black",
                                "Malaria vivax-motivated" = "#F0A20E",  
                                "Malaria falciparum-motivated" = "#004B6F", "LF-motivated" = "#0094C6" )) +
  scale_fill_manual(values = c("Rao-motivated" = "#FF4F00", "Single-pathogen motivated" = "black",
                               "Malaria vivax-motivated" = "#F0A20E",  
                               "Malaria falciparum-motivated" = "#004B6F", "LF-motivated" = "#0094C6" )) +
  guides(color =guide_legend(title="Strategy")) +
  xlab("Pathogen") +
  ggtitle("a.") +
  theme(legend.position = "bottom") +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "right",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 16, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))  +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80))  +
  ylim(c(0, 80))
target_goal 

# Optimal vs Rao
##########################
# blocks required results 
# LF: optimal (27, 35), rao (31, 45)
# Strongy: optimal (29, 37), rao (36, 47)
# p vivax: optimal (8, 20), rao (9, 25)
# p falciparum: optimal (5, 13), rao: (5, 17)

# make into line plot 
mean = c(11, 14, 6, 8, 32, 39, 31, 34                          )
lower = c(8, 9, 5, 5, 29, 36, 27, 31)
upper = c(20, 25, 13, 17 , 37, 47, 35, 45)
path = c( "Malaria vivax", "Malaria vivax", "Malaria falciparum", "Malaria falciparum", 
          "Strongyloides", "Strongyloides", "LF", "LF")
strat = c( "Rao-motivated", "Single-pathogen motivated", "Rao-motivated", "Single-pathogen motivated", "Rao-motivated", "Single-pathogen motivated", "Rao-motivated", "Single-pathogen motivated")

paths = data.frame(mean, lower, upper, path, strat)
paths$path = factor(paths$path, levels = c( "LF", "Malaria falciparum", "Malaria vivax", "Strongyloides"))


bootstrap = ggplot(data = paths, aes(x = path, y = mean, col = strat)) + 
  theme_light()+
  geom_point( position=position_dodge(.5), cex = 2) +
  # ylab("Area under efficiency curve")+
  ylab("Clusters to target 50% of disease")+
  theme(axis.text=element_text(size=12)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.4,
                position=position_dodge(.5), lwd= 1) + 
  xlab("Pathogen") +
  ggtitle("b.") +
  theme(axis.text = element_text(size = 13, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        legend.position = "none",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        #  legend.margin=margin(1,1.5,0.5,0.5,unit = "line"),
        strip.text = element_text(colour = "black", size = 15, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  guides(color =guide_legend(title="Strategy")) +
  scale_color_manual(values = c("black", "#FF4F00")) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80))  +
  ylim(c(0, 80))
bootstrap

figure_2_cam = plot_grid(target_goal, bootstrap, rel_widths = c(.7, .45))



########## Map out- Figure 1 #################################################################
#################### Map out rao #################################
# Number of observations per psuid 
cambodia_count = 
  cambodia_seropositivity %>%
  count(psuid) 
count_n = print(cambodia_count$n)

# number of cases by psuid
cambodia_case = 
  cambodia_seropositivity %>%
  group_by(psuid) %>%
  summarize(across(tetanus_unprotected:lymph_any, ~sum(.x, na.rm = TRUE))) %>%
  mutate(across(tetanus_unprotected:lymph_any, ~ .x / count_n))

# gps coordinates with seropositivity 
cambodia_locations = left_join(cambodia_case, gps_dat %>% mutate(psuid = psuid + 1), by = "psuid") %>%
  left_join(gps_cambodia, by = "dhsclust" ) %>%
  dplyr::select(psuid, tetanus_unprotected, strong, malaria_f, malaria_v, lymph_any, dhsclust,
                dhslat, dhslon) %>%
  left_join(method1, by = "psuid") 
head(cambodia_locations)


# Correlation between pathogens
library(wesanderson)
pal <- wes_palette("Darjeeling1", 100, type = "continuous")

corr_plot = cambodia_locations %>%
  dplyr::select(strong, malaria_f, malaria_v, lymph_any)
head(corr_plot)
ggplot(data = corr_plot) +
  geom_histogram(aes(x = lymph_any)) + xlim(c(0,1 )) + theme_bw()


colnames(corr_plot) = c("Strongyloides", "Malaria falciparum", "Malaria vivax", "LF")
Correlation <- round(cor(corr_plot), 2)
cor = ggcorrplot(Correlation,
                 hc.order = TRUE,
                 type = "lower",
                 outline.color = "white", 
                 colors = c("blue4", "white", "greenyellow")) +
  scale_fill_viridis(option = "B", direction = 1, end = .5, begin = 1,
                     limits = c(0, 1),  # Set continuous scale limits
                     guide = guide_colorbar(title = "Correlation", 
                                            title.position = "top",
                                            label.position = "right",
                                            ticks.colour = "white",
                                            ticks.linewidth = 1,
                                            frame.colour = "white")) 
cor 




######################################################################
# Figure 1 with masking now 
######################################################################
# Prepare masking 
pal <- wes_palette("Zissou1", 100, type = "continuous")

### Load in population density data, from 2012 at this link: https://hub.worldpop.org/geodata/summary?id=45540
cambodia_pop_density = read_stars(here::here("projects/6-multipathogen-burden/data/cambodia/khm_pd_2012_1km_UNadj.tif"))

# Make into a shapefile
cam_pd =st_as_sf(cambodia_pop_density)

ggplot(data = cam_pd %>% filter(log(khm_pd_2012_1km_UNadj.tif) > 1)) +
  geom_sf(aes(col = log(khm_pd_2012_1km_UNadj.tif))) +
  theme_minimal() +
  scale_color_viridis_c(begin = .2, guide = guide_colorbar(title = "Log pop density") ) +
  ggtitle("Cambodia population density")


# Create masked shape file 
masked_density_sf = cam_pd %>%
  mutate(log_popdensity = log(khm_pd_2012_1km_UNadj.tif)) %>%
  filter(log_popdensity > log(1))
crs_shapefile <- st_crs(masked_density_sf)


# Convert the shapefile into a raster 
#if (!identical(crs_raster, crs_shapefile)) {
#  sf_cam_pd <- st_transform(cam_pd, crs_raster)
#}

# Convert the shapefile to a spatial object compatible with raster
sf_cam_masked <- as(masked_density_sf, "Spatial")



#######################################
# For individual pathogen serology 
maps_i = c( "malaria_f", "malaria_v", "lymph_any", "strong")
path_names = c( "Malaria falciparum", "Malaria vivax",  "Lymphatic filariasis", "Strongyloides")

plot_list <- list()

for(i in 1:length(maps_i)) {
  
  # Map out 
  # Load country-level boundaries for the world
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  # Filter for Cambodia
  cambodia_map <- countries[countries$name == "Cambodia", ]
  
  # Prepare inputs for mapping
  dbgps_block = cambodia_locations %>% 
    mutate(lat = dhslat, lon = dhslon) %>% 
    dplyr::select(psuid, lon, lat) %>%
    dplyr::distinct()  
  
  dbgps = cambodia_locations %>%
    mutate(lat = dhslat, lon = dhslon) %>%
    ungroup() %>%
    dplyr::select(psuid, !!maps_i[i], lat, lon) %>%
    distinct() %>%
    drop_na()
  
  # Create bounding box from the shapefile. 
  b_boxgrid <- st_make_grid(cambodia_map,
                            n = c(100, 100),
                            what = "centers"  
  )
  
  # Transform cluster centroids into shapefile. 
  dbgps_sf_utm = dbgps_block %>%
    ungroup() %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    dplyr::select(geometry) %>%  
    st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km")
  
  # Identify a buffer of 10 km.
  cl_buff <- st_buffer(dbgps_sf_utm, dist = 200) %>%  
    summarise(geometry = st_union(geometry)) %>%
    st_cast("POLYGON") %>% 
    st_transform(crs = 4326) 
  
  # Fit the model dynamically based on the current map
 # formula_string <- paste(maps_i[i], "~ lat + lon + Matern(1|lat + lon)")
#  b_fit_rao <- spaMM::fitme(
#    as.formula(formula_string),
#    data = dbgps,
#    family = gaussian(link = "identity")
#  )
  
  dbgps = dbgps %>%
    mutate( cases = round(!!sym(maps_i[i]) * 100, digits = 0), 
            noncases = 100 - cases)
  print(head(dbgps))
  
  # Fit the model dynamically based on the current map
  # formula_string <- paste(cbind(cases, noncases), "~ lat + lon + Matern(1|lat + lon)")
  b_fit_rao <- spaMM::fitme(
    #  as.formula(formula_string),
    formula = cbind(cases, noncases) ~  lat + lon + Matern(1|lat + lon), 
    data = dbgps,
    family = binomial(link = "logit"))
  
  # Predictions and rasterization
  b_preds_rao <- get_grid_preds(input_grid = b_boxgrid, spamm_model_fit = b_fit_rao)
  b_preds_rao_coords <- st_coordinates(b_preds_rao)
  
  b_preds_rao_raster <- points_to_raster(
    x = b_preds_rao_coords[,1],
    y = b_preds_rao_coords[,2],
    z = b_preds_rao$pred,
    mask1 = cambodia_map,
    mask2 = cl_buff,
    crop1 = cambodia_map
  )
  
  # Mask the raster with the shapefile to retain only the areas inside the polygons of the shapefile
  masked_raster <- mask(b_preds_rao_raster, sf_cam_masked)
  
  b_preds_rao_raster2 <- raster_to_tibble(masked_raster) %>%
    drop_na(value)  %>%
    mutate(value = value * 100)
  
  print(summary(b_preds_rao_raster2$value))
  
  # Plotting
  map <- ggplot(data = b_preds_rao_raster2 %>%
                  mutate(value = replace(value, value < 0, 0)) %>%
                  drop_na(value) %>%
                  mutate(value = replace(value, value > 100, 100))) +  
    geom_tile(aes(x = x, y = y, fill = value), na.rm = TRUE) +
    geom_tile(aes(x = x, y = y, fill = value), na.rm = TRUE) +
  #  geom_sf(data = cambodia_map, fill = NA, color = "black", size = 15, lwd = 1) +  # Black border for country
    ggtitle(path_names[i]) +  # Now using path_names[i] directly
    coord_sf(crs = 4326) +
    #  scale_fill_viridis(option = "D", direction = 1,
    #                     limits = c(0, 100),
    #                     na.value = NA,
    #                     breaks = seq(0, 100, by = 10), 
    #                     guide = guide_colorbar(
    #                       title = "Seroprevalence (%)",
    #                       title.position = "top",
    #                       title.hjust = 0.5,
    #                       direction = "horizontal",
    #                       label.position = "bottom",
    #                       label.hjust = 0.5,
    #                       barheight = unit(10, "pt"),
    #                       barwidth = unit(200, "pt"),
    #                       ticks.colour = "white",
    #                       ticks.linewidth = 1,
    #                       frame.colour = "white"
    #                     )
    #  ) +
    scale_fill_gradientn(colours = pal,
                         limits = c(0, 100),
                         na.value = NA,
                         breaks = seq(0, 100, by = 10), 
                         guide = guide_colorbar(
                           title = "Seroprevalence (%)",
                           title.position = "top",
                           title.hjust = 0.5,
                           direction = "horizontal",
                           label.position = "bottom",
                           label.hjust = 0.5,
                           barheight = unit(10, "pt"),
                           barwidth = unit(200, "pt"),
                           ticks.colour = "white",
                           ticks.linewidth = 1,
                           frame.colour = "white") ) +
    labs(x = "Latitude", y = "Longitude") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 16)
    )   +
    theme( plot.title = element_text(colour = "black", size = 13.5, face = "bold")) 
  
  plot_list[[i]] = map
}

sero_plot_cam <- plot_grid(plotlist = plot_list, ncol = 2)
sero_plot_cam

######################################################
# for rao 

maps_i = c("rao")
path_names = c("Rao's quadratic entropy")

plot_list <- list()

for(i in 1:length(maps_i)) {
  
  # Map out 
  # Load country-level boundaries for the world
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  # Filter for Cambodia
  cambodia_map <- countries[countries$name == "Cambodia", ]
  
  # Prepare inputs for mapping
  dbgps_block = cambodia_locations %>% 
    mutate(lat = dhslat, lon = dhslon) %>% 
    dplyr::select(psuid, lon, lat) %>%
    dplyr::distinct()  
  
  dbgps = cambodia_locations %>%
    mutate(lat = dhslat, lon = dhslon) %>%
    ungroup() %>%
    dplyr::select(psuid, !!maps_i[i], lat, lon) %>%
    distinct() %>%
    drop_na()
  
  # Create bounding box from the shapefile. 
  # This creates a grid over the Cambodia map, where the grid is size 100 x 100 
  b_boxgrid <- st_make_grid(cambodia_map,
                            n = c(100, 100),
                            what = "centers"  
  )
  
  # Transform cluster centroids into shapefile. 
  dbgps_sf_utm = dbgps_block %>%
    ungroup() %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% # st_as_sf converts lat and lon into coordinates 
    dplyr::select(geometry) %>%  
    st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km") # specifies they are in UTM coordinations, km insteaad of degrees
  
  # Identify a buffer of 10 km.
  cl_buff <- st_buffer(dbgps_sf_utm, dist = 200) %>%  # 200 m buffer around spatial points in dbgps_sf_utm
    summarise(geometry = st_union(geometry)) %>% # combines buffers into unified geometry 
    st_cast("POLYGON") %>% 
    st_transform(crs = 4326) # make into a polygon and concert back to WGS84
  
  dbgps = dbgps %>%
    mutate( log_rao = log(!!sym(maps_i[i]) ))
  print(head(dbgps$log_rao))
  
  
  # Fit the model dynamically based on the current map
 # formula_string <- paste(maps_i[i], "~ lat + lon + Matern(1|lat + lon)")
  formula_string <- paste("log_rao ~ lat + lon + Matern(1|lat + lon)")
  b_fit_rao <- spaMM::fitme(
   as.formula(formula_string),
   data = dbgps,
    family = gaussian(link = "identity")
  )
  
  # Predictions and rasterization
  b_preds_rao <- get_grid_preds(input_grid = b_boxgrid, spamm_model_fit = b_fit_rao) # generates predictions from the fitted model
  b_preds_rao_coords <- st_coordinates(b_preds_rao)
  
  b_preds_rao_raster <- points_to_raster(
    x = b_preds_rao_coords[,1],
    y = b_preds_rao_coords[,2],
    z = b_preds_rao$pred,
    mask1 = cambodia_map,
    mask2 = cl_buff,
    crop1 = cambodia_map
  ) # converts to a raster with the map nad buffer as spatial constriats 
  
  # Mask the raster with the shapefile to retain only the areas inside the polygons of the shapefile
  masked_raster <- mask(b_preds_rao_raster, sf_cam_masked)
  
  b_preds_rao_raster2 <- raster_to_tibble(masked_raster) %>%
    drop_na() %>%
    mutate(value_exp = exp(value))
  
  max_val = print(max( b_preds_rao_raster2$value_exp))
    
  # Plotting
  map <- ggplot(data = b_preds_rao_raster2 %>%
                  drop_na(value_exp)) +  
    geom_tile(aes(x = x, y = y, fill = value_exp), na.rm = TRUE) +
  #  geom_sf(data = cambodia_map, fill = NA, color = "black", size = 15, lwd = 1) +  # Black border for country
    ggtitle(path_names[i]) +  # Now using path_names[i] directly
    coord_sf(crs = 4326) +
    scale_fill_gradientn(colours = pal,
                         limits = c(0, max_val),
                         # na.value = NA,
                         breaks = seq(0, max_val, by = 1), 
                         guide = guide_colorbar(
                           title = "Rao's quadratic entropy",
                           title.position = "top",
                           title.hjust = 0.5,
                           direction = "horizontal",
                           label.position = "bottom",
                           label.hjust = 0.5,
                           barheight = unit(10, "pt"),
                           barwidth = unit(200, "pt"),
                           ticks.colour = "white",
                           ticks.linewidth = 1,
                           frame.colour = "white") ) +
    #  scale_fill_viridis(option = "D", direction = 1,
    #                     limits = c(0, max_val),
    #                     na.value = NA,
    #                     breaks = seq(0, max_val, by = 1), 
    #                     guide = guide_colorbar(
    #                       title = "Rao diversity",
    #                       title.position = "top",
    #                       title.hjust = 0.5,
    #                       direction = "horizontal",
    #                       label.position = "bottom",
    #                       label.hjust = 0.5,
    #                       barheight = unit(10, "pt"),
    #                       barwidth = unit(200, "pt"),
    #                       ticks.colour = "white",
    #                       ticks.linewidth = 1,
    #                       frame.colour = "white"
    #                     )
    #  ) +
    labs(x = "Latitude", y = "Longitude") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 16)
    )  +
    theme( plot.title = element_text(colour = "black", size = 14.5, face = "bold")) 
  
  plot_list[[i]] = map
}


rao_plot_cam <- plot_grid(plotlist = plot_list)
rao_plot_cam


# Figure 
multi_cam =plot_grid(rao_plot_cam, cor, ncol = 1, rel_heights = c(.5, .2), labels  = c("a", "b"))
multi_cam
sero_four_plot_cam = plot_grid(sero_plot_cam, labels = "c")
sero_four_plot_cam
figure_1_cam = plot_grid(multi_cam, sero_four_plot_cam)
figure_1_cam


##################################################################################################
# March 29, this concludes the figures ###########################################################

# Figure 1 - combine locations
plot_grid(figure_1_cam, figure_1_bangl, ncol = 1)
# 1400 x 1500

# Figure 2 


# Figure 3
plot_grid(figure_2_cam, figure_2_bangl, 
          ncol = 1)







# who did you write book with
























































#######################################
# for serology 
maps_i = c("strong", "malaria_f", "malaria_v", "lymph_any")
path_names = c("Strongyloides", "Malaria falciparum", "Malaria vivax", "Lymphatic filariasis")

head(cambodia_locations)

plot_list <- list()

for(i in 1:length(maps_i)) {
  
  # Map out 
  # Load country-level boundaries for the world
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  # Filter for Cambodia
  cambodia_map <- countries[countries$name == "Cambodia", ]
  
  # Prepare inputs for mapping
  dbgps_block = cambodia_locations %>% 
    mutate(lat = dhslat, lon = dhslon) %>% 
    dplyr::select(psuid, lon, lat) %>%
    dplyr::distinct()  
  
  dbgps = cambodia_locations %>%
    mutate(lat = dhslat, lon = dhslon) %>%
    ungroup() %>%
    dplyr::select(psuid, !!maps_i[i], lat, lon) %>%
    distinct() %>%
    drop_na()
  
  print(head(dbgps))
  
  # Create bounding box from the shapefile. 
  b_boxgrid <- st_make_grid(cambodia_map,
                            n = c(100, 100),
                            what = "centers"  
  )
  
  # Transform cluster centroids into shapefile. 
  dbgps_sf_utm = dbgps_block %>%
    ungroup() %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    dplyr::select(geometry) %>%  
    st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km")
  
  # Identify a buffer of 10 km.
  cl_buff <- st_buffer(dbgps_sf_utm, dist = 200) %>%  
    summarise(geometry = st_union(geometry)) %>%
    st_cast("POLYGON") %>% 
    st_transform(crs = 4326) 
  
#  # Fit the model dynamically based on the current map
#  formula_string <- paste(maps_i[i], "~ lat + lon + Matern(1|lat + lon)")
#  b_fit_rao <- spaMM::fitme(
#    as.formula(formula_string),
#    data = dbgps,
#    family = gaussian(link = "identity")
#  )
  
  dbgps = dbgps %>%
    mutate( cases = round(!!sym(maps_i[i]) * 100, digits = 0), 
            noncases = 100 - cases)
  
  print(head(dbgps))
  
  # Fit the model dynamically based on the current map
 # formula_string <- paste(cbind(cases, noncases), "~ lat + lon + Matern(1|lat + lon)")
  b_fit_rao <- spaMM::fitme(
  #  as.formula(formula_string),
    formula = cbind(cases, noncases) ~  lat + lon + Matern(1|lat + lon), 
    data = dbgps,
    family = binomial(link = "logit")
  )
  
  
  
  # Predictions and rasterization
  b_preds_rao <- get_grid_preds(input_grid = b_boxgrid, spamm_model_fit = b_fit_rao)
  b_preds_rao_coords <- st_coordinates(b_preds_rao)
  
  b_preds_rao_raster <- points_to_raster(
    x = b_preds_rao_coords[,1],
    y = b_preds_rao_coords[,2],
    z = b_preds_rao$pred,
    mask1 = cambodia_map,
    mask2 = cl_buff,
    crop1 = cambodia_map
  )
  
  b_preds_rao_raster2 <- raster_to_tibble(b_preds_rao_raster) %>%
    drop_na(value)  %>%
    mutate(value = value * 100)
  
  # Plotting
  map <- ggplot(data = b_preds_rao_raster2 %>%
                  mutate(value = replace(value, value < 0, 0)) %>%
                  drop_na(value) %>%
                  mutate(value = replace(value, value > 100, 100))) +  
    geom_tile(aes(x = x, y = y, fill = value), na.rm = TRUE) +
    ggtitle(path_names[i]) +  # Now using path_names[i] directly
    coord_sf(crs = 4326) +
    scale_fill_viridis(option = "D", direction = 1,
                       limits = c(0, 100),
                       na.value = NA,
                       breaks = seq(0, 100, by = 10), 
                       guide = guide_colorbar(
                         title = "Seroprevalence (%)",
                         title.position = "top",
                         title.hjust = 0.5,
                         direction = "horizontal",
                         label.position = "bottom",
                         label.hjust = 0.5,
                         barheight = unit(10, "pt"),
                         barwidth = unit(200, "pt"),
                         ticks.colour = "white",
                         ticks.linewidth = 1,
                         frame.colour = "white"
                       )
    ) +
    labs(x = "Latitude", y = "Longitude") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 16)
    )  
  
  plot_list[[i]] = map
}

sero_plot <- plot_grid(plotlist = plot_list, ncol = 2)
sero_plot

######################################################
# for rao 
head(cambodia_locations)
mean_rao <- mean(cambodia_locations$rao, na.rm = TRUE)  # Calculate mean
var_rao <- var(cambodia_locations$rao, na.rm = TRUE)    # Calculate variance

ggplot(data = cambodia_locations) +
  geom_histogram(aes(x = log(rao))) +
  geom_text(aes(x = mean_rao + 3, y = 30), 
            label = paste("Mean: ", round(mean_rao, 1)), 
            color = "red", size = 5, vjust = 4) +
  geom_text(aes(x = mean_rao + 3, y = 27), 
            label = paste("Variance: ", round(var_rao, 1)), 
            color = "blue", size = 5, vjust = 4)  +
  ylab("Frequency") +
  xlab("Log rao") + theme_bw() +
  xlim(c(-5, 5))

  
  
install.packages("wesanderson")
library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")


maps_i = c("rao")
path_names = c("Rao diversity")

plot_list <- list()

for(i in 1:length(maps_i)) {
  
  # Map out 
  # Load country-level boundaries for the world
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  # Filter for Cambodia
  cambodia_map <- countries[countries$name == "Cambodia", ]
  
  # Prepare inputs for mapping
  dbgps_block = cambodia_locations %>% 
    mutate(lat = dhslat, lon = dhslon) %>% 
    dplyr::select(psuid, lon, lat) %>%
    dplyr::distinct()  
  
  dbgps = cambodia_locations %>%
    mutate(lat = dhslat, lon = dhslon) %>%
    ungroup() %>%
    dplyr::select(psuid, rao, lat, lon) %>% # !!maps_i[i]
    distinct() %>%
    drop_na() %>%
    mutate(rao_standard = (rao/max(rao)*100) ) %>%
    mutate(rao_whole = round(rao_standard, digits = 0))
  
  print(summary(dbgps$rao_whole))
  print(mean(dbgps$rao_whole))
  print(var(dbgps$rao_whole))
  
  
  # Create bounding box from the shapefile. 
  # This creates a grid over the Cambodia map, where the grid is size 100 x 100 
  b_boxgrid <- st_make_grid(cambodia_map,
                            n = c(100, 100),
                            what = "centers"  
  )
  
  # Transform cluster centroids into shapefile. 
  dbgps_sf_utm = dbgps_block %>%
    ungroup() %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% # st_as_sf converts lat and lon into coordinates 
    dplyr::select(geometry) %>%  
    st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km") # specifies they are in UTM coordinations, km insteaad of degrees
  
  # Identify a buffer of 10 km.
  cl_buff <- st_buffer(dbgps_sf_utm, dist = 200) %>%  # 200 m buffer around spatial points in dbgps_sf_utm
    summarise(geometry = st_union(geometry)) %>% # combines buffers into unified geometry 
    st_cast("POLYGON") %>% 
    st_transform(crs = 4326) # make into a polygon and concert back to WGS84
  
  # Fit the model dynamically based on the current map
 # formula_string <- paste(maps_i[i], "~ lat + lon + Matern(1|lat + lon)")
#  b_fit_rao <- spaMM::fitme(
#    as.formula(formula_string),
#    data = dbgps,
#    family = gaussian(link = "identity")
#  )
  
  # Fit the model dynamically based on the current map
  formula_string <- paste("rao_whole ~ lat + lon + Matern(1|lat + lon)")
  b_fit_rao <- spaMM::fitme(
    as.formula(formula_string),
    data = dbgps,
  #  family = negbin2(link = "log")
    family = gaussian(link = "identity")
  #  family = poisson(link = "log")
  )
  
  # Predictions and rasterization
  b_preds_rao <- get_grid_preds(input_grid = b_boxgrid, spamm_model_fit = b_fit_rao) # generates predictions from the fitted model
  b_preds_rao_coords <- st_coordinates(b_preds_rao)
  
  b_preds_rao_raster <- points_to_raster(
    x = b_preds_rao_coords[,1],
    y = b_preds_rao_coords[,2],
    z = b_preds_rao$pred,
    mask1 = cambodia_map,
    mask2 = cl_buff,
    crop1 = cambodia_map
  ) # converts to a raster with the map nad buffer as spatial constriats 
  
  b_preds_rao_raster2 <- raster_to_tibble(b_preds_rao_raster) %>%
    drop_na(value) %>%
    mutate(value = value/100*max(dbgps$rao))
  
  max_val = print(max(b_preds_rao_raster2$value))
  
  # Plotting
  map <- ggplot(data = b_preds_rao_raster2 %>%
                  drop_na(value)) +  
    geom_tile(aes(x = x, y = y, fill = value), na.rm = TRUE) +
    ggtitle(path_names[i]) +  # Now using path_names[i] directly
    coord_sf(crs = 4326) +
  #  scale_fill_viridis(option = "D", direction = 1,
  #                     limits = c(0, max_val),
  #                     na.value = NA,
  #                     breaks = seq(0, max_val, by = 1), 
  #                     guide = guide_colorbar(
  #                       title = "Rao diversity",
  #                       title.position = "top",
  #                       title.hjust = 0.5,
  #                       direction = "horizontal",
  #                       label.position = "bottom",
  #                       label.hjust = 0.5,
  #                       barheight = unit(10, "pt"),
  #                       barwidth = unit(200, "pt"),
  #                       ticks.colour = "white",
  #                       ticks.linewidth = 1,
  #                       frame.colour = "white"
  #                     )
  #  ) +
    scale_fill_gradientn(colours = pal,
                        # limits = c(0, max_val),
                        # na.value = NA,
                        # breaks = seq(0, max_val, by = 1), 
                         guide = guide_colorbar(
                           title = "Rao diversity",
                           title.position = "top",
                           title.hjust = 0.5,
                           direction = "horizontal",
                           label.position = "bottom",
                           label.hjust = 0.5,
                           barheight = unit(10, "pt"),
                           barwidth = unit(200, "pt"),
                           ticks.colour = "white",
                           ticks.linewidth = 1,
                           frame.colour = "white") ) +
    labs(x = "Latitude", y = "Longitude") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 16)
    )  
  
  plot_list[[i]] = map
}


rao_plot <- plot_grid(plotlist = plot_list)
rao_plot


# Figure 
multi =plot_grid(rao_plot, cor, ncol = 1, rel_heights = c(.5, .2))
figure_1 = plot_grid(multi, sero_plot)
figure_1




#################################
####### try masking #############

countries <- ne_countries(scale = "medium", returnclass = "sf")
# Filter for Cambodia
cambodia_map <- countries[countries$name == "Cambodia", ]

# Prepare inputs for mapping
dbgps_block = cambodia_locations %>% 
  mutate(lat = dhslat, lon = dhslon) %>% 
  dplyr::select(psuid, lon, lat) %>%
  dplyr::distinct()  

dbgps = cambodia_locations %>%
  mutate(lat = dhslat, lon = dhslon) %>%
  ungroup() %>%
  dplyr::select(psuid, rao, lat, lon) %>%
  distinct() %>%
  drop_na()

# Create bounding box from the shapefile. 
# This creates a grid over the Cambodia map, where the grid is size 100 x 100 
b_boxgrid <- st_make_grid(cambodia_map,
                          n = c(100, 100),
                          what = "centers"  
)

# Transform cluster centroids into shapefile. 
dbgps_sf_utm = dbgps_block %>%
  ungroup() %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% # st_as_sf converts lat and lon into coordinates 
  dplyr::select(geometry) %>%  
  st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km") # specifies they are in UTM coordinations, km insteaad of degrees

# Identify a buffer of 10 km.
cl_buff <- st_buffer(dbgps_sf_utm, dist = 200) %>%  # 200 m buffer around spatial points in dbgps_sf_utm
  summarise(geometry = st_union(geometry)) %>% # combines buffers into unified geometry 
  st_cast("POLYGON") %>% 
  st_transform(crs = 4326) # make into a polygon and concert back to WGS84

# Fit the model dynamically based on the current map
formula_string <- paste("rao ~ lat + lon + Matern(1|lat + lon)")
b_fit_rao <- spaMM::fitme(
  as.formula(formula_string),
  data = dbgps,
  family = gaussian(link = "identity")
)

# Predictions and rasterization
b_preds_rao <- get_grid_preds(input_grid = b_boxgrid, spamm_model_fit = b_fit_rao) # generates predictions from the fitted model
b_preds_rao_coords <- st_coordinates(b_preds_rao)

b_preds_rao_raster <- points_to_raster(
  x = b_preds_rao_coords[,1],
  y = b_preds_rao_coords[,2],
  z = b_preds_rao$pred,
  mask1 = cambodia_map,
  mask2 = cl_buff,
  crop1 = cambodia_map
) # converts to a raster with the map nad buffer as spatial constraints 
b_preds_rao_raster
crs_raster = crs(b_preds_rao_raster)

### Load in population density data, from 2012 at this link: https://hub.worldpop.org/geodata/summary?id=45540
cambodia_pop_density = read_stars(here::here("projects/6-multipathogen-burden/data/cambodia/khm_pd_2012_1km_UNadj.tif"))

# Make into a shapefile
cam_pd =st_as_sf(cambodia_pop_density)

ggplot(data = cam_pd %>% filter(log(khm_pd_2012_1km_UNadj.tif) > 1)) +
  geom_sf(aes(col = log(khm_pd_2012_1km_UNadj.tif))) +
  theme_minimal() +
  scale_color_viridis_c(begin = .2, guide = guide_colorbar(title = "Log pop density") ) +
  ggtitle("Cambodia population density")


# Create masked shape file 
masked_density_sf = cam_pd %>%
  mutate(log_popdensity = log(khm_pd_2012_1km_UNadj.tif)) %>%
  filter(log_popdensity > log(5))
crs_shapefile <- st_crs(masked_density_sf)
  
  
# Convert the shapefile into a raster 
if (!identical(crs_raster, crs_shapefile)) {
  sf_cam_pd <- st_transform(cam_pd, crs_raster)
}

# Convert the shapefile to a spatial object compatible with raster
sf_cam_masked <- as(masked_density_sf, "Spatial")

# Mask the raster with the shapefile to retain only the areas inside the polygons of the shapefile
masked_raster <- mask(b_preds_rao_raster, sf_cam_masked)

b_preds_rao_raster2 <- raster_to_tibble(masked_raster) %>%
  drop_na()

max_val = print(max(b_preds_rao_raster2$value))

# Plotting
map <- ggplot(data = b_preds_rao_raster2 ) +  
  geom_sf(data = cambodia_map, fill = NA, color = "black", size = 15, lwd = 1.2) +  # Black border for country
  geom_tile(aes(x = x, y = y, fill = value), na.rm = TRUE) +
  ggtitle(path_names[i]) +  # Now using path_names[i] directly
  coord_sf(crs = 4326) +
  scale_fill_viridis(option = "D", direction = 1,
                     limits = c(0, max_val),
                     na.value = NA,
                     breaks = seq(0, max_val, by = 1), 
                     guide = guide_colorbar(
                       title = "Rao diversity",
                       title.position = "top",
                       title.hjust = 0.5,
                       direction = "horizontal",
                       label.position = "bottom",
                       label.hjust = 0.5,
                       barheight = unit(10, "pt"),
                       barwidth = unit(200, "pt"),
                       ticks.colour = "white",
                       ticks.linewidth = 1,
                       frame.colour = "white"
                     )
  )  + 
  theme_bw() +
  ylab("Longtitude") + xlab("Latitude") +
  theme(
    legend.position = "bottom",
    plot.tag = element_text(face = "bold", size = 16) 
  )   


map




######################################################################
# Figure 1 with masking now 
######################################################################
# Prepare masking 
pal <- wes_palette("Zissou1", 100, type = "continuous")

### Load in population density data, from 2012 at this link: https://hub.worldpop.org/geodata/summary?id=45540
cambodia_pop_density = read_stars(here::here("projects/6-multipathogen-burden/data/cambodia/khm_pd_2012_1km_UNadj.tif"))

# Make into a shapefile
cam_pd =st_as_sf(cambodia_pop_density)

ggplot(data = cam_pd %>% filter(log(khm_pd_2012_1km_UNadj.tif) > 1)) +
  geom_sf(aes(col = log(khm_pd_2012_1km_UNadj.tif))) +
  theme_minimal() +
  scale_color_viridis_c(begin = .2, guide = guide_colorbar(title = "Log pop density") ) +
  ggtitle("Cambodia population density")


# Create masked shape file 
masked_density_sf = cam_pd %>%
  mutate(log_popdensity = log(khm_pd_2012_1km_UNadj.tif)) %>%
  filter(log_popdensity > log(1))
crs_shapefile <- st_crs(masked_density_sf)


# Convert the shapefile into a raster 
if (!identical(crs_raster, crs_shapefile)) {
  sf_cam_pd <- st_transform(cam_pd, crs_raster)
}

# Convert the shapefile to a spatial object compatible with raster
sf_cam_masked <- as(masked_density_sf, "Spatial")



#######################################
# for serology 
maps_i = c("strong", "malaria_f", "malaria_v", "lymph_any")
path_names = c("Strongyloides", "Malaria falciparum", "Malaria vivax", "Lymphatic filariasis")

plot_list <- list()

for(i in 1:length(maps_i)) {
  
  # Map out 
  # Load country-level boundaries for the world
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  # Filter for Cambodia
  cambodia_map <- countries[countries$name == "Cambodia", ]
  
  # Prepare inputs for mapping
  dbgps_block = cambodia_locations %>% 
    mutate(lat = dhslat, lon = dhslon) %>% 
    dplyr::select(psuid, lon, lat) %>%
    dplyr::distinct()  
  
  dbgps = cambodia_locations %>%
    mutate(lat = dhslat, lon = dhslon) %>%
    ungroup() %>%
    dplyr::select(psuid, !!maps_i[i], lat, lon) %>%
    distinct() %>%
    drop_na()
  
  # Create bounding box from the shapefile. 
  b_boxgrid <- st_make_grid(cambodia_map,
                            n = c(100, 100),
                            what = "centers"  
  )
  
  # Transform cluster centroids into shapefile. 
  dbgps_sf_utm = dbgps_block %>%
    ungroup() %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    dplyr::select(geometry) %>%  
    st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km")
  
  # Identify a buffer of 10 km.
  cl_buff <- st_buffer(dbgps_sf_utm, dist = 200) %>%  
    summarise(geometry = st_union(geometry)) %>%
    st_cast("POLYGON") %>% 
    st_transform(crs = 4326) 
  
  # Fit the model dynamically based on the current map
  formula_string <- paste(maps_i[i], "~ lat + lon + Matern(1|lat + lon)")
  b_fit_rao <- spaMM::fitme(
    as.formula(formula_string),
    data = dbgps,
    family = gaussian(link = "identity")
  )
  
  # Predictions and rasterization
  b_preds_rao <- get_grid_preds(input_grid = b_boxgrid, spamm_model_fit = b_fit_rao)
  b_preds_rao_coords <- st_coordinates(b_preds_rao)
  
  b_preds_rao_raster <- points_to_raster(
    x = b_preds_rao_coords[,1],
    y = b_preds_rao_coords[,2],
    z = b_preds_rao$pred,
    mask1 = cambodia_map,
    mask2 = cl_buff,
    crop1 = cambodia_map
  )
  
  # Mask the raster with the shapefile to retain only the areas inside the polygons of the shapefile
  masked_raster <- mask(b_preds_rao_raster, sf_cam_masked)
  
  b_preds_rao_raster2 <- raster_to_tibble(masked_raster) %>%
    drop_na(value)  %>%
    mutate(value = value * 100)
  
  # Plotting
  map <- ggplot(data = b_preds_rao_raster2 %>%
                  mutate(value = replace(value, value < 0, 0)) %>%
                  drop_na(value) %>%
                  mutate(value = replace(value, value > 100, 100))) +  
    geom_tile(aes(x = x, y = y, fill = value), na.rm = TRUE) +
    geom_tile(aes(x = x, y = y, fill = value), na.rm = TRUE) +
   # geom_sf(data = cambodia_map, fill = NA, color = "black", size = 15, lwd = 1) +  # Black border for country
    ggtitle(path_names[i]) +  # Now using path_names[i] directly
    coord_sf(crs = 4326) +
  #  scale_fill_viridis(option = "D", direction = 1,
  #                     limits = c(0, 100),
  #                     na.value = NA,
  #                     breaks = seq(0, 100, by = 10), 
  #                     guide = guide_colorbar(
  #                       title = "Seroprevalence (%)",
  #                       title.position = "top",
  #                       title.hjust = 0.5,
  #                       direction = "horizontal",
  #                       label.position = "bottom",
  #                       label.hjust = 0.5,
  #                       barheight = unit(10, "pt"),
  #                       barwidth = unit(200, "pt"),
  #                       ticks.colour = "white",
  #                       ticks.linewidth = 1,
  #                       frame.colour = "white"
  #                     )
  #  ) +
    scale_fill_gradientn(colours = pal,
                         limits = c(0, 100),
                         na.value = NA,
                          breaks = seq(0, 100, by = 10), 
                         guide = guide_colorbar(
                           title = "Seroprevalence (%)",
                           title.position = "top",
                           title.hjust = 0.5,
                           direction = "horizontal",
                           label.position = "bottom",
                           label.hjust = 0.5,
                           barheight = unit(10, "pt"),
                           barwidth = unit(200, "pt"),
                           ticks.colour = "white",
                           ticks.linewidth = 1,
                           frame.colour = "white") ) +
    labs(x = "Latitude", y = "Longitude") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 16)
    )   +
    theme( plot.title = element_text(colour = "black", size = 13.5, face = "bold")) 
  
  plot_list[[i]] = map
}

sero_plot <- plot_grid(plotlist = plot_list, ncol = 2)
sero_plot

######################################################
# for rao 

maps_i = c("rao")
path_names = c("Rao diversity")

plot_list <- list()

for(i in 1:length(maps_i)) {
  
  # Map out 
  # Load country-level boundaries for the world
  countries <- ne_countries(scale = "medium", returnclass = "sf")
  # Filter for Cambodia
  cambodia_map <- countries[countries$name == "Cambodia", ]
  
  # Prepare inputs for mapping
  dbgps_block = cambodia_locations %>% 
    mutate(lat = dhslat, lon = dhslon) %>% 
    dplyr::select(psuid, lon, lat) %>%
    dplyr::distinct()  
  
  dbgps = cambodia_locations %>%
    mutate(lat = dhslat, lon = dhslon) %>%
    ungroup() %>%
    dplyr::select(psuid, !!maps_i[i], lat, lon) %>%
    distinct() %>%
    drop_na()
  
  # Create bounding box from the shapefile. 
  # This creates a grid over the Cambodia map, where the grid is size 100 x 100 
  b_boxgrid <- st_make_grid(cambodia_map,
                            n = c(100, 100),
                            what = "centers"  
  )
  
  # Transform cluster centroids into shapefile. 
  dbgps_sf_utm = dbgps_block %>%
    ungroup() %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% # st_as_sf converts lat and lon into coordinates 
    dplyr::select(geometry) %>%  
    st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km") # specifies they are in UTM coordinations, km insteaad of degrees
  
  # Identify a buffer of 10 km.
  cl_buff <- st_buffer(dbgps_sf_utm, dist = 200) %>%  # 200 m buffer around spatial points in dbgps_sf_utm
    summarise(geometry = st_union(geometry)) %>% # combines buffers into unified geometry 
    st_cast("POLYGON") %>% 
    st_transform(crs = 4326) # make into a polygon and concert back to WGS84
  
  # Fit the model dynamically based on the current map
  formula_string <- paste(maps_i[i], "~ lat + lon + Matern(1|lat + lon)")
  b_fit_rao <- spaMM::fitme(
    as.formula(formula_string),
    data = dbgps,
    family = gaussian(link = "identity")
  )
  
  # Predictions and rasterization
  b_preds_rao <- get_grid_preds(input_grid = b_boxgrid, spamm_model_fit = b_fit_rao) # generates predictions from the fitted model
  b_preds_rao_coords <- st_coordinates(b_preds_rao)
  
  b_preds_rao_raster <- points_to_raster(
    x = b_preds_rao_coords[,1],
    y = b_preds_rao_coords[,2],
    z = b_preds_rao$pred,
    mask1 = cambodia_map,
    mask2 = cl_buff,
    crop1 = cambodia_map
  ) # converts to a raster with the map nad buffer as spatial constriats 
  
  # Mask the raster with the shapefile to retain only the areas inside the polygons of the shapefile
  masked_raster <- mask(b_preds_rao_raster, sf_cam_masked)
  
  b_preds_rao_raster2 <- raster_to_tibble(masked_raster) %>%
    drop_na()
  
  max_val = print(max(b_preds_rao_raster2$value))

  
  # Plotting
  map <- ggplot(data = b_preds_rao_raster2 %>%
                  drop_na(value)) +  
    geom_tile(aes(x = x, y = y, fill = value), na.rm = TRUE) +
  #  geom_sf(data = cambodia_map, fill = NA, color = "black", size = 15, lwd = 1) +  # Black border for country
    ggtitle(path_names[i]) +  # Now using path_names[i] directly
    coord_sf(crs = 4326) +
    scale_fill_gradientn(colours = pal,
                          limits = c(0, max_val),
                         # na.value = NA,
                         breaks = seq(0, max_val, by = 1), 
                         guide = guide_colorbar(
                           title = "Rao diversity",
                           title.position = "top",
                           title.hjust = 0.5,
                           direction = "horizontal",
                           label.position = "bottom",
                           label.hjust = 0.5,
                           barheight = unit(10, "pt"),
                           barwidth = unit(200, "pt"),
                           ticks.colour = "white",
                           ticks.linewidth = 1,
                           frame.colour = "white") ) +
  #  scale_fill_viridis(option = "D", direction = 1,
  #                     limits = c(0, max_val),
  #                     na.value = NA,
  #                     breaks = seq(0, max_val, by = 1), 
  #                     guide = guide_colorbar(
  #                       title = "Rao diversity",
  #                       title.position = "top",
  #                       title.hjust = 0.5,
  #                       direction = "horizontal",
  #                       label.position = "bottom",
  #                       label.hjust = 0.5,
  #                       barheight = unit(10, "pt"),
  #                       barwidth = unit(200, "pt"),
  #                       ticks.colour = "white",
  #                       ticks.linewidth = 1,
  #                       frame.colour = "white"
  #                     )
  #  ) +
    labs(x = "Latitude", y = "Longitude") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 16)
    )  +
    theme( plot.title = element_text(colour = "black", size = 14.5, face = "bold")) 
  
  plot_list[[i]] = map
}


rao_plot <- plot_grid(plotlist = plot_list)
rao_plot


# Figure 
multi =plot_grid(rao_plot, cor, ncol = 1, rel_heights = c(.5, .2), labels  = c("a", "b"))
multi
sero_four_plot = plot_grid(sero_plot, labels = "c")
sero_four_plot
figure_1 = plot_grid(multi, sero_four_plot)
figure_1

"#A8D8FF"

(6.3 + 6.9 + 8.4 + 22.9 + 8.8 + 6.8 + 7.1 + 7.1)/ 5
(6.3 + 6.9 + 8.4 + 22.9 + 8.8 + 6.8 + 7.1 + 7.1)/ (21.8 + 23.4 + 40.7 + 26.3 + 24.8)
  
21.50*.907 + 104.50*.093
21.50*.991 + 104.50*.009

.0017*172084 + (1 - .0017)*40.66
115*.66 + 40.66*.33

25000000*.005 * (2150/15600000)

14900000* (2150/15600000)

40.66*.5 + .06*.5
