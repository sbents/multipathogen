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
  mutate(tetanus_unprotected = ifelse(ttmb >= 100, 0,1 )) %>% 
  mutate(lymph14 = ifelse(bm14 >= 65, 1, 0 )) %>%
  mutate(lymph123 = ifelse(wb123 >= 115, 1, 0 )) %>%
  mutate(lymph33 = ifelse(bm33 >= 966, 1, 0 )) %>% #Pat suggested dropping this 
  mutate(strong = ifelse(nie  >= 792, 1, 0 )) %>%
  mutate(toxo = ifelse(sag2a  >= 159, 1, 0 )) %>%
  mutate(cyst = ifelse(t24  >= 486, 1, 0 )) %>%
  mutate(malaria_f = ifelse(pfmsp19  >= 343, 1, 0 )) %>%
  mutate(malaria_v = ifelse(pvmsp19  >= 196, 1, 0 )) %>%
  dplyr::select(region, psuid, age, tetanus_unprotected, lymph14, lymph123, 
                lymph33, strong, toxo, cyst, malaria_f, malaria_v) %>%
  mutate(malaria_any = ifelse(malaria_f == 1 | malaria_v == 1, 1, 0)) %>%
 # mutate(lymph_any = ifelse(lymph14 == 1 | lymph123 == 1| lymph33 == 1, 1, 0)) %>% 9/2/2025
  mutate(lymph_any = ifelse(lymph14 == 1 & lymph123 == 1, 1, 0)) %>%
 # mutate(lymph_any = ifelse(lymph123 == 1, 1, 0)) %>%
  mutate(psuid = psuid + 1)

head(cambodia_seropositivity)
round(colSums(cambodia_seropositivity[4:14])/length(cambodia_seropositivity), digits = 2) 

psuid_n = length(unique(cambodia_seropositivity$psuid))

# check peeople per psuid 
people_per_cluster = cambodia_seropositivity %>%
  mutate(person = 1) %>%
  group_by(psuid) %>%
  count(person)
summary(people_per_cluster$n)

# Age distribution
summary(cambodia_seropositivity$age)


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

# plot cluster level distribution of prevalences 
histo = parasite_prevalence_block %>%
  mutate(pathogen = replace(pathogen, pathogen == "lymph_any", "Lymphatic filariasis")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "malaria_f", "P. falciparum")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "malaria_v", "P. vivax")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "strong", "Strongyloides")) %>%
  mutate(pathogen = factor(pathogen, levels = c("P. falciparum", "P. vivax", 
                                                "Strongyloides", "Lymphatic filariasis")))
head(histo)
ggplot(data = histo ) +
  geom_histogram(aes(x = fraction), fill = "blue4", bins = 15) +
  facet_wrap(vars(pathogen)) + theme_bw()+
  theme(axis.text = element_text(size = 11, color = "black"),
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
        strip.text = element_text(colour = "black", size = 16, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  ylab("Count") + xlab("Prevalence") +
  ggtitle("Cambodia")




#######
prev_overall = parasite_prevalence_block %>%
  mutate(prev = even_prevalence*100) %>%
  distinct(pathogen, prev)
head(prev_overall)



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
head(method1)
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
pt = .75 # targeting threshold for disease, was .5

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
    mutate(above = ifelse(target < pt, 1, 2))
  order_malf = print(malf_motivated$psuid)
  
  print(head(malf_motivated))
  
  malf_50 <- which(malf_motivated$above == 2)[1]
  
  # rao motivated 
  rao_malf_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_f") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
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
    mutate(above = ifelse(target < pt, 1, 2))
  
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
    mutate(above = ifelse(target < pt, 1, 2))
  
  malf_lf_50 <- which(malf_lf_motivated$above == 2)[1]
  
  # Wealth 
 malf_wealth_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_f") %>%
    left_join(rao_wealth_cambodia, by = "psuid") %>%
    arrange(hv271) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  malf = ggplot(data = malf_motivated, aes(x = block_label/psuid_n*100, y = target*100)) +
    geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "P. falciparum"), cex = 2) +
    geom_line(data = rao_malf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = malf_malv_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "P. vivax"), cex = 2) + 
    geom_line(data = malf_lf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "LF"), cex = 2) + 
  #  geom_line(data =  malf_wealth_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "wealth"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("P. falciparum" = "black", "Rao" = "#FF4F00",
           #                       "P. vivax" = "#F0A20E", "LF" = "#0094C6")) + # "wealth" = "darkblue"
           "P. vivax" = "#F0A20E", "LF" = "darkslategray2")) + # "wealth" = "darkblue"
    guides(color =guide_legend(title="Strategy")) +
    ggtitle("      P. falciparum")+
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
    mutate(above = ifelse(target < pt, 1, 2))
  strong_50 <- which(strong_motivated$above == 2)[1]
  
  rao_strong_motivated = compare_strategy %>% 
    filter(pathogen == "strong") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
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
    mutate(above = ifelse(target < pt, 1, 2))
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
    mutate(above = ifelse(target < pt, 1, 2))
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
    mutate(above = ifelse(target < pt, 1, 2))
  strong_malf_50 <- which(strong_malf_motivated$above == 2)[1]
  
  # Wealth 
  strong_wealth_motivated = compare_strategy %>% 
    filter(pathogen == "strong") %>%
    left_join(rao_wealth_cambodia, by = "psuid") %>%
    arrange(hv271) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  strong = ggplot(data = strong_motivated, aes(x = block_label/psuid_n*100, y = target*100)) +
    geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "Strongyloides"), cex = 2) +
    geom_line(data = rao_strong_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = strong_malv_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "P. vivax"), cex = 2) + 
    geom_line(data = strong_lf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "LF"), cex = 2) + 
    geom_line(data = strong_malf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "P. falciparum"), cex = 2) + 
  #  geom_line(data = strong_wealth_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "wealth"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("Strongyloides" = "black", "Rao" = "#FF4F00", "P. falciparum" = "deepskyblue3" ,
                                  "P. vivax" = "#F0A20E", "LF" = "darkslategray2" )) + #"wealth" = "darkblue"
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
    mutate(above = ifelse(target < pt, 1, 2))
  lymph_50 <- which(lymph_motivated$above == 2)[1]
  
  rao_lymph_motivated = compare_strategy %>% 
    filter(pathogen == "lymph_any") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(target< pt, 1, 2))
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
    mutate(above = ifelse(target < pt, 1, 2))
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
    mutate(above = ifelse(target < pt, 1, 2))
  lf_malv_50 <- which(lf_malv_motivated$above == 2)[1]
  
  # wealth
  lf_wealth_motivated = compare_strategy %>% 
    filter(pathogen == "lymph_any") %>%
    left_join(rao_wealth_cambodia, by = "psuid") %>%
    arrange(hv271) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  
  lymph = ggplot(data = lymph_motivated, aes(x = block_label/psuid_n*100, y = target*100)) +
    geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "LF"), cex = 2) +
    geom_line(data = rao_lymph_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = lf_malv_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "P. vivax"), cex = 2) + 
    geom_line(data = lf_malf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "P. falciparum"), cex = 2) + 
  #  geom_line(data = lf_wealth_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "wealth"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("LF" = "black", "Rao" = "#FF4F00",
                             #     "P. vivax" = "#F0A20E" , "P. falciparum" = "#004B6F")) + # 
                             "P. vivax" = "#F0A20E" , "P. falciparum" = "deepskyblue3")) + # 
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
    mutate(above = ifelse(target < pt, 1, 2))
  malv_50 <- which(malv_motivated$above == 2)[1]
  
  rao_malv_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_v") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(target< pt, 1, 2))
  rao_malv_50 <- which(rao_malv_motivated$above == 2)[1]
  
  malv_lf_motivated = compare_strategy %>% 
    mutate(customorder = factor(psuid), levels = order_lf) %>%
    filter(pathogen == "malaria_v") %>%
    arrange(order_lf) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
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
    mutate(above = ifelse(target < pt, 1, 2))
  malv_malf_50 <- which(malv_malf_motivated$above == 2)[1]
  
  # wealth
  malv_wealth_motivated = compare_strategy %>% 
    filter(pathogen == "malaria_v") %>%
    left_join(rao_wealth_cambodia, by = "psuid") %>%
    arrange(hv271) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, psuid_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  
  malv = ggplot(data = malv_motivated, aes(x = block_label/psuid_n*100, y = target*100)) +
    geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "P. vivax"), cex = 2) +
    geom_line(data = rao_malv_motivated , aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = malv_lf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "LF"), cex = 2) + 
    geom_line(data = malv_malf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "P. falciparum"), cex = 2) + 
  #  geom_line(data = malv_wealth_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "wealth"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("P. vivax" = "black", "Rao" = "#FF4F00", "P. falciparum" = "deepskyblue3" ,
                               "LF" =  "darkslategray2" )) + # "wealth" = "darkblue"
  #  scale_color_viridis_d(option = "H", direction = 1, end = 1, begin = .6) +
    guides(color =guide_legend(title="Strategy")) +
    ggtitle("      P. vivax")+
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


print(block_50[[1]]) # Method 1 
print(block_50[[2]]) # Method 3 

# Plot efficiency of single-pathogen strategies vs rao strategies. 
#method1_efficiency <- plot_grid(efficiency_list = efficiency_list[[1]], ncol = 2)
method1_efficiency <- efficiency_list[[1]]
method1_efficiency 

#method2_efficiency <- plot_grid(efficiency_list = efficiency_list[[2]], ncol = 2)
#method2_efficiency


# make legend for strategy 
legend1 = ggplot(data = parasite_prevalence_block) +
  geom_line(aes(x = psuid, y = presence, col = "Rao-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = fraction, col = "P. falciparum-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = even_prevalence, col = "P. vivax-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = pop_presence, col = "LF-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = pop_presence* even_prevalence, col = "Single-pathogen motivated"), cex = 2) + 
  scale_color_manual(values = c("P. vivax-motivated" = "#F0A20E" , "Rao-motivated" = "#FF4F00", "P. falciparum-motivated" =  "deepskyblue3",
                                "LF-motivated" =  "darkslategray2", "Single-pathogen motivated" = "black")) +
  guides(color =guide_legend(title="Strategy"))  +
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, color = "black"))+ 
  geom_blank() 
legend1
legend_only <- gtable::gtable_filter(ggplotGrob(legend1), "guide-box")

fig2_cam_eff <- cowplot::plot_grid(
  method1_efficiency,   # Main plot
  legend_only,           # Extracted legend
  ncol = 1,              # Layout with one column
  rel_heights = c(1, 0.1)  # Adjust relative heights to make space for the legend
)

fig2_cam_eff

# Display the final plot

#leg = get_legend(legend1)
#legend_only <- ggpubr::as_ggplot(leg) + theme(plot.margin = margin(0,0,0,0))
#legend_only <- cowplot::get_legend(legend1)
#legend_only
#grid.newpage()
#legend_only = grid.draw(leg)
#plot_grid(method1_efficiency , legend_only, rel_widths = c(1, .2), labels = "a")
#legend_only = grid.draw(leg)

# was patchwork
#fig2_cam_eff <- method1_efficiency + legend_only + 
#  plot_layout(widths = c(.03, .1)) & 
#  theme(legend.position = "bottom")




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
    mutate(above = ifelse(target < pt, 1, 2))
  
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
    mutate(above = ifelse(target < pt, 1, 2))
  
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
    mutate(above = ifelse(target < pt, 1, 2))
  
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
    mutate(above = ifelse(target < pt, 1, 2))
  
  rao_lf_50_sim <- which( rao_lf_motivated $above == 2)[1]
  sim_50[[i]] <- rao_lf_50_sim
  
}

lf_vector <- unlist(sim_50)
lf_sim <- data.frame(lf_vector)

############################## make a boxplot with them all
print(block_50[[1]])

head(malaria_f_sim)
mf_boxplot = malaria_f_sim %>%
  mutate(pathogen = "P. falciparum") %>%
#  mutate(Rao = 8) %>%       # 50% targeting 
#  mutate(Optimal = 6) %>%   # 50% targeting 
#  mutate(P.vivax = 67) %>%  # 50% targeting 
#  mutate(LF = 49)           # 50% targeting 
#  mutate(Rao = 19) %>%       # 75% targeting 
#  mutate(Optimal = 11) %>%   # 75% targeting 
#  mutate(P.vivax = 75) %>%   # 75% targeting 
#  mutate(LF = 73)            # 75% targeting 
  mutate(Rao = 13) %>%       # 75% targeting, new LF
  mutate(Optimal = 11) %>%   # 75% targeting, new LF
  mutate(P.vivax = 75) %>%   # 75% targeting, new LF 
  mutate(LF = 80)            # 75% targeting, new LF 
head(mf_boxplot)

mv_boxplot = malaria_v_sim %>%
  mutate(pathogen = "P. vivax") %>%
 # mutate(Rao = 14) %>%            # 50% targeting 
 # mutate(Optimal = 11) %>%        # 50% targeting 
 # mutate(P.falciparum = 49) %>%   # 50% targeting 
 # mutate(LF = 47)                 # 50% targeting 
 #  mutate(Rao = 36) %>%            # 75% targeting 
#   mutate(Optimal = 23) %>%        # 75% targeting 
#   mutate(P.falciparum = 74) %>%   # 75% targeting 
#   mutate(LF = 74)                 # 75% targeting 
  mutate(Rao = 30) %>%            # 75% targeting, new LF 
  mutate(Optimal = 23) %>%        # 75% targeting, new LF 
  mutate(P.falciparum = 74) %>%   # 75% targeting, new LF 
  mutate(LF = 78)                 # 75% targeting, new LF 
head(mv_boxplot)

strong_boxplot = strong_sim %>%
  mutate(pathogen = "Strongyloides") %>%
 # mutate(Rao = 39) %>%               # 50% targeting 
 #  mutate(Optimal = 32) %>%          # 50% targeting 
 #  mutate(P.falciparum = 54) %>%     # 50% targeting 
 #  mutate(P.vivax = 50) %>%          # 50% targeting 
 #  mutate(LF = 50)                   # 50% targeting 
#  mutate(Rao = 64) %>%                # 75% targeting 
#  mutate(Optimal = 56) %>%            # 75% targeting 
#  mutate(P.falciparum = 75) %>%       # 75% targeting 
#  mutate(P.vivax = 72) %>%            # 75% targeting 
#  mutate(LF = 81)                     # 75% targeting 
  mutate(Rao = 74) %>%                # 75% targeting, new LF 
  mutate(Optimal = 56) %>%            # 75% targeting, new LF 
  mutate(P.falciparum = 75) %>%       # 75% targeting, new LF 
  mutate(P.vivax = 72) %>%            # 75% targeting, new LF 
  mutate(LF = 79)                     # 75% targeting, new LF 
head(strong_boxplot)

lf_boxplot = lf_sim %>%
  mutate(pathogen = "LF") %>%
#  mutate(Rao = 34) %>%               # 50% targeting
#  mutate(Optimal = 31) %>%           # 50% targeting
#  mutate(P.falciparum = 50) %>%      # 50% targeting
#  mutate(P.vivax = 52)               # 50% targeting
#  mutate(Rao = 59) %>%                # 75% targeting
#  mutate(Optimal = 55) %>%            # 75% targeting
#  mutate(P.falciparum = 77) %>%       # 75% targeting
#   mutate(P.vivax = 72)              # 75% targeting
  mutate(Rao = 31) %>%                # 75% targeting, new LF
  mutate(Optimal = 18) %>%            # 75% targeting, new LF
  mutate(P.falciparum = 77) %>%       # 75% targeting, new LF
  mutate(P.vivax = 72)               # 75% targeting, new LF

# efficient
# lf 
(72- 31)/72
# strong 
(72-64)/64
# mv
(74 - 30)/74
# mf
(75 - 13)/75


target_goal = ggplot(data = mv_boxplot) + 
  geom_violin(aes(x = pathogen, y = malaria_v_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
 # geom_boxplot(aes(x = pathogen, y = malaria_v_vector), fill = "gray95", width=.6) + 
  geom_point(aes(x = pathogen, y = Rao),  col = "#FF4F00", cex = 6)  + #.5, rao  
  geom_point(aes(x = pathogen, y = Optimal), col = "black",  cex = 6)  + 
  geom_point(aes(x = pathogen, y = P.falciparum), col = "deepskyblue3",  cex = 6)  + 
  geom_point(aes(x = pathogen, y = LF), col = "darkslategray2",  cex = 6)  +
  theme_bw() +
  #ylab("Clusters to target 50% of disease") +
  ylab("Percent of clusters to target 75% of disease (%)") +
  geom_violin(data = mf_boxplot, aes(x = pathogen, y = malaria_f_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
 # geom_boxplot(data = mf_boxplot, aes(x = pathogen, y = malaria_f_vector), fill = "gray95", width=.6) + 
  geom_point(data = mf_boxplot, aes(x = pathogen, y = Rao), col = "#FF4F00",  cex = 6)  + #.5, rao  
  geom_point(data = mf_boxplot, aes(x = pathogen, y = Optimal), col = "black",  cex = 6)  + 
  geom_point(data = mf_boxplot, aes(x = pathogen, y = P.vivax), col = "#F0A20E",  cex = 6)  + 
  geom_point(data = mf_boxplot, aes(x = pathogen, y = LF), col = "darkslategray2",  cex = 6)  +
  geom_violin(data = strong_boxplot, aes(x = pathogen, y = strong_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
 # geom_boxplot(data = strong_boxplot, aes(x = pathogen, y = strong_vector), fill = "gray95", width=.6) + 
  geom_point(data = strong_boxplot, aes(x = pathogen, y = Optimal, col = "Single-pathogen motivated"),  cex = 6)  + 
  geom_point(data = strong_boxplot, aes(x = pathogen, y = P.vivax, col = "P. vivax-motivated") , cex = 6)  + 
  geom_point(data = strong_boxplot, aes(x = pathogen, y = P.falciparum, col = "P. falciparum-motivated"),  cex = 6)  + 
  geom_point(data = strong_boxplot, aes(x = pathogen, y = LF, col = "LF-motivated"),   cex = 6)  +
  geom_point(data = strong_boxplot, aes(x = pathogen, y = Rao, col = "Rao-motivated"),  cex = 6)  + #.5, rao  
  geom_violin(data = lf_boxplot, aes(x = pathogen, y = lf_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
 # geom_boxplot(data = lf_boxplot, aes(x = pathogen, y = lf_vector), fill = "gray95", width=.6) + 
  geom_point(data = lf_boxplot, aes(x = pathogen, y = Rao), col = "#FF4F00",  cex = 6)  + #.5, rao  
  geom_point(data = lf_boxplot, aes(x = pathogen, y = Optimal), col = "black", cex = 6)  + 
  geom_point(data = lf_boxplot, aes(x = pathogen, y = P.vivax), col = "#F0A20E",  cex = 6)  + 
  geom_point(data = lf_boxplot, aes(x = pathogen, y = P.falciparum), col = "deepskyblue3",  cex = 6)   +
  scale_color_manual(values = c("Rao-motivated" = "#FF4F00", "Single-pathogen motivated" = "black",
                                "P. vivax-motivated" = "#F0A20E",  
                                "P. falciparum-motivated" = "deepskyblue3", "LF-motivated" = "darkslategray2" )) +
  scale_fill_manual(values = c("Rao-motivated" = "#FF4F00", "Single-pathogen motivated" = "black",
                               "P. vivax-motivated" = "#F0A20E",  
                               "P. falciparum-motivated" = "deepskyblue3", "LF-motivated" = "darkslategray2" )) +
  guides(color =guide_legend(title="Strategy")) +
  xlab("Pathogen") +
  ggtitle("b") +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 11, color = "black"),
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
        strip.text = element_text(colour = "black", size = 16, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA))  +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80))  +
  ylim(c(0, 100))
target_goal 


# Optimal vs Rao
# Bootstrapping 

print(unique(parasite_prevalence_block$pathogen))

store_statistic_opt = list()

for(j in 1:200){
  
  sample_replace = parasite_prevalence_block %>%
    filter(pathogen == "lymph_any") %>% # arbtriary to sample from
    group_modify(~ .x[sample(nrow(.x), replace = TRUE), ]) %>%
    left_join(parasite_prevalence_block, sample_replace, by = "psuid" ,  relationship =
                "many-to-many") %>% # join so we obtain random sample of psuid across all pathogens
    mutate(pathogen = pathogen.y, presence = presence.y,  even_prevalence = even_prevalence.y,
           pop_presence = pop_presence.y, fraction = fraction.y) %>%
    dplyr::select(pathogen, psuid, presence, pathogen, even_prevalence, pop_presence, fraction) %>%
    arrange(pathogen) %>%
    mutate(psuid_replace = rep(seq(1, 100, 1),  4))
  
  head(sample_replace)
  
  # Calculate Rao diversity at each block. 
  product_results = list()
  
  for (i in 1:psuid_n) {
    
    print(i)
    # Subset the data for the current block
    block_subset <- sample_replace %>%
      mutate(fraction_mirrored = fraction) %>%
      ungroup() %>%
      filter(psuid_replace== i) %>%
      dplyr::select(fraction_mirrored, fraction)
    
    # Get unique combinations of pathogen prevalences 
    fractions <- block_subset$fraction
    unique_combinations <- t(combn(fractions, 2, simplify = TRUE))  # Get pairs (no duplicates)
    
    # Compute products for each pair
    products <- apply(unique_combinations, 1, prod)
    
    products_df <- data.frame(
      psuid_replace = i, 
      fraction_mirrored = unique_combinations[, 1],
      fraction = unique_combinations[, 2],
      product = products
    )
    # Store the results in the list
    print(products_df)
    product_results[[i]] <- products_df
    
  }
  
  # Combine all results into a single data frame
  product_results_replace <- do.call(rbind, product_results)
  
  # Sum by block
  block_prevalence_final_replace = product_results_replace %>%
    group_by(psuid_replace) %>%
    summarize(across(product, ~sum(.x)))
  head(block_prevalence_final_replace)
  
  # Calculate expected rao diversity for a given block.
  block_subset_replace = parasite_prevalence_block  %>% 
    mutate(even_prevalence_mirrored = even_prevalence) %>%
    ungroup() %>% filter(psuid == 1) %>%
    dplyr::select(even_prevalence_mirrored, even_prevalence)
  
  even_prevalences = block_subset_replace$even_prevalence
  unique_combinations <- t(combn(even_prevalences, 2, simplify = TRUE))  # Get pairs (no duplicates)
  
  # Compute products for each pair
  products <- apply(unique_combinations, 1, prod)
  null_rao = sum(products)
  
  # Compare observed rao to expected rao and add gps data. 
  recode_psuid = sample_replace %>%
    dplyr::select(psuid, psuid_replace) %>%
    distinct()
  
  method1_replace = block_prevalence_final_replace %>%
    mutate(rao = product/null_rao) %>%
    distinct() %>%
    left_join(recode_psuid, by = "psuid_replace")
  head(method1_replace)
  
  compare_strategy_replace = left_join(sample_replace, method1_replace, by = "psuid_replace")
  
  path = print(unique(compare_strategy_replace))
  
  # Malaria f strategy 
  # motivated by itself
  # TO DO: change disease here 
  
 # malf_motivated = compare_strategy_replace %>% 
#    filter(pathogen == "lymph_any")  %>%
#    arrange(-fraction) %>% # Arrange from highest prevalence block to lowest 
#    mutate(total_fraction = sum(fraction)) %>%
#    mutate(cum = cumsum(fraction)) %>%
#    mutate(target = cum/total_fraction) %>%
#    mutate(block_label = seq(1, psuid_n, 1)) %>%
 #   mutate(above = ifelse(target < pt, 1, 2))
#  auc_malf <- trapz(malf_motivated$block_label/psuid_n, malf_motivated$target)
  
#  blocks_needed <- which( malf_motivated$above == 2)[1]
# rao motivated 
   rao_malf_motivated =  compare_strategy_replace %>% 
    filter(pathogen == "malaria_v") %>%
      arrange(-rao) %>%
      mutate(total_fraction = sum(fraction)) %>%
      mutate(cum = cumsum(fraction)) %>%
      mutate(target = cum/total_fraction) %>%
      mutate(block_label = seq(1, psuid_n, 1)) %>%
#   #  # mutate(above = ifelse(cum < pt, 1, 2))
      mutate(above = ifelse(target < pt, 1, 2))
  #  auc_malf_rao <- trapz(rao_malf_motivated$block_label/psuid_n, rao_malf_motivated$target)
   
  blocks_needed <- which( rao_malf_motivated$above == 2)[1]
  
  # statistic = auc_malf  
  statistic = blocks_needed
  print(statistic)
  
  store_statistic_opt[[j]]= statistic
  
} 


auc_vector_opt <- unlist(store_statistic_opt)
auc_list_opt <- data.frame(auc_vector_opt) %>%
  drop_na()
head(auc_list_opt)


quantiles = c(.025, .5, .975)
quant_025 = quantile(auc_list_opt$auc_vector_opt, probs = quantiles[1])
quant_50 = quantile(auc_list_opt$auc_vector_opt, probs = quantiles[2])
quant_975 = quantile(auc_list_opt$auc_vector_opt, probs = quantiles[3])


# it should be within that reality yeah?? 
# so you shoudl always be able to get to 75% of disease

# LF new definiton, targeting goal of 75% 
# LF:        optimal:  11 (18, 24)   rao: 27 (18, 33)
# strongy:   optimal:  56 (54, 59)   rao:  73 (65, 79)
# malaria f: optimal:  11 (8, 15)    rao:  17 (11, 27)
# malaria v: optimal:  24 (19, 29)   rao:  29 (21, 40)


# manual output for 75% 
# LF:        optimal:  55 (51, 59)   rao: 59 (54, 63)
# strongy:   optimal:  56 (54, 59)   rao:  64 (60, 68)
# malaria f: optimal:  11 (8, 15)    rao:  14 (8, 19)
# malaria v: optimal:  24 (19, 29)   rao:  35 (24, 54)

##########################
# blocks required results 
# LF: optimal (27, 35), rao (31, 45)
# Strongy: optimal (29, 37), rao (36, 47)
# p vivax: optimal (8, 20), rao (9, 25)
# p falciparum: optimal (5, 13), rao: (5, 17)

# make into line plot 
# these results are for 50% targeting goal 
#mean = c(11, 14, 6, 8, 32, 39, 31, 34                          )
#lower = c(8, 9, 5, 5, 29, 36, 27, 31)
#upper = c(20, 25, 13, 17 , 37, 47, 35, 45)

## 75% targeting goal 
mean = c(  29, 24, 17,  11,  73, 56, 27, 18)
lower = c( 21, 19, 11,   8,  65,  54, 18, 11)
upper = c( 40, 29, 27,  15,   79, 59, 33, 24)
path = c( "P. vivax", "P. vivax", "P. falciparum", "P. falciparum", 
          "Strongyloides", "Strongyloides", "LF", "LF")
strat = c( "Rao-motivated", "Single-pathogen motivated", "Rao-motivated", "Single-pathogen motivated", "Rao-motivated", "Single-pathogen motivated", "Rao-motivated", "Single-pathogen motivated")


paths = data.frame(mean, lower, upper, path, strat)
paths$path = factor(paths$path, levels = c( "LF", "P. falciparum", "P. vivax", "Strongyloides"))


bootstrap = ggplot(data = paths, aes(x = path, y = (mean/psuid_n)*100, col = strat)) + 
  theme_light()+
  geom_point( position=position_dodge(.5), cex = 2) +
  # ylab("Area under efficiency curve")+
  ylab("Percent of clusters to target 75% of disease (%)")+
  theme(axis.text=element_text(size=12)) +
  geom_errorbar(aes(ymin= (lower/psuid_n)*100, ymax= (upper/psuid_n)*100), width=.4,
                position=position_dodge(.5), lwd= 1) + 
  xlab("Pathogen") +
  ggtitle("c") +
  theme(axis.text = element_text(size = 11, color = "black"),
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
  scale_color_manual(values = c(   "#FF4F00", "black")) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80))  +
  ylim(c(0, 100))
bootstrap

figure_2_cam = plot_grid(target_goal, bootstrap, rel_widths = c(.75, .65))
figure_2_cam 


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

head(gps_cambodia)

# gps coordinates with seropositivity 
cambodia_locations = left_join(cambodia_case, gps_dat %>% mutate(psuid = psuid + 1), by = "psuid") %>%
  left_join(gps_cambodia, by = "dhsclust" ) %>%
  dplyr::select(psuid, tetanus_unprotected, strong, malaria_f, malaria_v, lymph_any, dhsclust,
                dhslat, dhslon) %>%
  left_join(method1, by = "psuid") 
head(cambodia_locations)

# Mean radii 
head(gps_cambodia)
coords <- as.matrix(cambodia_locations[, c("dhslon", "dhslat")])
# Calculate full distance matrix (in meters)
dist_matrix <- distm(coords, fun = distHaversine)  # Returns meters

# Set diagonal to NA to ignore self-distance
diag(dist_matrix) <- NA

# Calculate mean distance to k nearest neighbors (e.g., 3)
k <- 3
mean_knn_distances <- apply(dist_matrix, 1, function(x) {
  mean(sort(x, na.last = NA)[1:k])
})
# Add to original data frame (converted to km)
mean_knn_radius_km <- mean_knn_distances / 1000
summary(mean_knn_radius_km )


ggplot(data = cambodia_locations) +
  geom_point(aes(x =dhslat, y = dhslon, col =  malaria_f), cex = 4)
ggplot(data = cambodia_locations) +
  geom_point(aes(x =dhslat, y = dhslon, col =  malaria_v), cex = 4)
ggplot(data = cambodia_locations) +
  geom_point(aes(x =dhslat, y = dhslon, col =  lymph_any), cex = 4)


# Correlation between pathogens
pal <- wes_palette("Darjeeling1", 100, type = "continuous")

corr_plot = cambodia_locations %>%
  dplyr::select(strong, malaria_f, malaria_v, lymph_any) 
head(corr_plot)
ggplot(data = corr_plot) +
  geom_histogram(aes(x = lymph_any)) + xlim(c(0,1 )) + theme_bw()


colnames(corr_plot) = c("Strongyloides", "P. falciparum", "P. vivax", "LF")
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
path_names = c( "P. falciparum", "P. vivax",  "Lymphatic filariasis", "Strongyloides")

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
  
  plot_materncor <- plot_matern_corr(b_fit_rao) 
  print(plot_materncor)
  
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
    geom_sf(data = cambodia_map, fill = NA, color = "black", size = 15, lwd = 1) +  # Black border for country
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
path_names = c("Rao's quadratic index")

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
  cl_buff <- st_buffer(dbgps_sf_utm, dist =100) %>%  # 200 m buffer around spatial points in dbgps_sf_utm
    summarise(geometry = st_union(geometry)) %>% # combines buffers into unified geometry 
    st_cast("POLYGON") %>% 
    st_transform(crs = 4326) # make into a polygon and concert back to WGS84
  
  dbgps = dbgps %>%
    mutate(rao = ifelse(rao == 0, .001, rao)) %>%
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
  
  plot_materncor <- plot_matern_corr(b_fit_rao) 
  print(plot_materncor)
  
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
    geom_sf(data = cambodia_map, fill = NA, color = "black", size = 15, lwd = 1) +  # Black border for country
    ggtitle(path_names[i]) +  # Now using path_names[i] directly
    coord_sf(crs = 4326) +
    scale_fill_gradientn(colours = pal,
                         limits = c(0, max_val),
                         # na.value = NA,
                         breaks = seq(0, max_val, by = 1), 
                         guide = guide_colorbar(
                           title = "Rao's quadratic index",
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
fig2_cam_eff
plot_grid(fig2_cam_eff, figure_2_cam, labels = c("a"), 
          nrow = 2, rel_heights =c(.5, .3))
# 1200 x 1300



# Figure 3
figure_2_cam
plot_grid(figure_2_cam, figure_2_bangl, 
          ncol = 1)


# actual feasible path
# gps coordinates with seropositivity 
cambodia_locations_order = left_join(cambodia_case, gps_dat %>% mutate(psuid = psuid + 1), by = "psuid") %>%
  left_join(gps_cambodia, by = "dhsclust" ) %>%
  dplyr::select(psuid, tetanus_unprotected, strong, malaria_f, malaria_v, lymph_any, dhsclust,
                dhslat, dhslon) %>%
  left_join(method1, by = "psuid") %>%
  arrange(-rao) %>% 
  mutate(order = row_number()) %>%
  mutate(rao = ifelse(rao == 0, .00001, rao))
head(cambodia_locations_order)

ggplot(data = cambodia_locations_order, aes(x = dhslat, y = dhslon, col = rao)) +
  geom_point()

ggplot(cambodia_locations_order %>% filter(order >  74), aes(x = dhslon, y = dhslat)) +
  geom_path(aes(group = 1), color = "blue", size = 1) +  # Connect points
  geom_point(color = "red", size = 3) +  # Plot points
  geom_text(aes(label = round(rao, 2)), vjust = -1) +  # Label with rao
  labs(title = "Path from Highest to Lowest rao: Quartile 4",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

ggplot(cambodia_locations_order %>% filter(rao > 0), aes(x = dhslon, y = dhslat)) +
  geom_path(aes(group = 1), color = "blue", size = 1) +  # Connect points
  geom_point(color = "red", size = 3) +  # Plot points
  geom_text(aes(label = round(rao, 2)), vjust = -1) +  # Label with rao
  labs(title = "Path from Highest to Lowest rao",
       x = "Longitude", y = "Latitude") +
  theme_minimal()



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

################ add sampling sites 
gps_cambodia =   readr::read_csv(file = here("projects/6-multipathogen-burden/data/cambodia", "cambodia_ea_dhs.csv")) 
head(gps_cambodia)

countries <- ne_countries(scale = "medium", returnclass = "sf")
# Filter for Cambodia
cambodia_map <- countries[countries$name == "Cambodia", ]

cam_sampling = ggplot() + 
  geom_sf(data = cambodia_map, fill = NA, color = "black", size = 15, lwd = 1) +
  theme_bw() +
  geom_point(data = gps_cambodia, aes(x = dhslon, y = dhslat), col = "darkred", fill = "red", cex = 3, pch = 23) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("Cambodia")
cam_sampling

# Download Cambodia GADM data at level 2
cam_districts_sp <- getData("GADM", country = "KHM", level = 2)

# Convert to sf
sf_cam_districts <- st_as_sf(cam_districts_sp)



################## Check cluster level prevalence 

# Number of non-NA measurements assessed total by pathogen 
parasite_denom = cambodia_seropositivity  %>%
  pivot_longer(cols=c('tetanus_unprotected', 'lymph14', 'lymph123', 'lymph33', 'strong', 'toxo', 'cyst',
                      'malaria_f', 'malaria_v', 'malaria_any',  'lymph_any'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(psuid, pathogen) %>% # assess by overall prevalence in population
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
head(parasite_numer)

prev_psuid = left_join(parasite_denom, parasite_numer, by =  c("psuid", "pathogen")) %>%
  mutate(prev_psuid = presence/non_na_count) %>%
  filter(pathogen == "lymph_any")
summary(prev_psuid$prev_psuid)
head(prev_psuid)

