# Compare diversity metrics: Rao vs Shannon vs Gini Simpson

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
  mutate(lymph_any = ifelse(lymph14 == 1 & lymph123 == 1, 1, 0)) %>%
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
  ungroup() %>%
  arrange(-rao)
head(method1)
summary(method1$rao)  

# rao order of locations 
rao_order = print(method1$psuid)

################ July 13 ################################################################
######################################

# Gini Simpson
gini_simpson = parasite_prevalence_block %>%
  group_by(psuid ) %>%
  mutate(psuid_presence = sum(presence)) %>%
  ungroup() %>%
  mutate(fraction_gini = presence/psuid_presence) %>%
  mutate(gs_component = fraction_gini^2) %>%
  group_by(psuid) %>%
  mutate(gs_entropy = sum(gs_component)) %>%
  ungroup() %>%
  mutate(gs_entropy_complement = 1 - gs_entropy ) %>%
  distinct(psuid, gs_entropy, gs_entropy_complement) %>%
  arrange(-gs_entropy_complement)
head(gini_simpson)
gs_order = print(gini_simpson$psuid) 

#######################################
# Shannon diversity 
shannon = parasite_prevalence_block %>%
  group_by(psuid ) %>%
  mutate(psuid_presence = sum(presence)) %>%
  ungroup() %>%
  mutate(fraction_shannon = presence/psuid_presence) %>%
  mutate(fraction_shannon = fraction_shannon + .0000000000001) %>%
  mutate(shannon_component = fraction_shannon*log2(fraction_shannon)) %>%
  group_by(psuid) %>%
  mutate(shannon_entropy = - sum(shannon_component)) %>%
  distinct(psuid, shannon_entropy) %>%
  arrange(-shannon_entropy)
shannon_order = print(shannon$psuid) 
print(shannon_order)

# Calculate the percent of disease targeted for each pathogen in Cambodia 
rao_ordered <- parasite_prevalence_block %>%
  group_by(pathogen) %>%
  mutate(customorder = factor(psuid, levels = rao_order)) %>%
  arrange(customorder) %>%  # arrange by the new factor
  mutate(
    cum = cumsum(fraction),
    cum_presence = cumsum(presence),
    total_fraction = sum(fraction),
    target = cum / total_fraction,
    block_label = seq(1, 100, 1)
  ) %>%
  ungroup() %>%
  mutate(Method = "Rao") %>%
  group_by(psuid) %>%
  mutate(people_treated = sum(cum_presence))

gs_ordered = parasite_prevalence_block %>%
  group_by(pathogen) %>%
  mutate(customorder = factor(psuid, levels = gs_order)) %>%
  arrange(customorder) %>%  # arrange by the new factor
  mutate(
    cum = cumsum(fraction),
    cum_presence = cumsum(presence),
    total_fraction = sum(fraction),
    target = cum / total_fraction,
    block_label = seq(1, 100, 1)
  ) %>%
  ungroup() %>%
  mutate(Method = "Gini-Simpson") %>%
  group_by(psuid) %>%
  mutate(people_treated = sum(cum_presence))

shannon_ordered = parasite_prevalence_block %>%
  group_by(pathogen) %>%
  mutate(customorder = factor(psuid, levels = shannon_order)) %>%
  arrange(customorder) %>%  # arrange by the new factor
  mutate(
    cum = cumsum(fraction),
    cum_presence = cumsum(presence),
    total_fraction = sum(fraction),
    target = cum / total_fraction,
    block_label = seq(1, 100, 1)
  ) %>%
  ungroup() %>%
  mutate(Method = "Shannon") %>%
  group_by(psuid) %>%
  mutate(people_treated = sum(cum_presence))

all_methods = rbind(rao_ordered, gs_ordered, shannon_ordered) %>%
  mutate(pathogen = replace(pathogen, pathogen == "lymph_any", "Lymphatic filariasis")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "malaria_f", "P. falciparum")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "malaria_v", "P. vivax")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "strong", "Strongyloides")) %>%
  mutate(pathogen = factor(pathogen, levels = c("P. falciparum", "P. vivax", 
                                                "Strongyloides", "Lymphatic filariasis")))
head(all_methods)

compare_a = ggplot(data = all_methods ) + 
  geom_line(aes(x = block_label, y = people_treated, col = Method), lwd = 1.9) +
  theme_bw() +
  xlab("Spatial cluster") +
  ylab("Number of people treated") +
  ggtitle("a.") +
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
  scale_color_manual(values = c("royalblue2", "#FF4F00", "gray70")) 


compare_b = ggplot(data = all_methods ) + 
  geom_line(aes(x = block_label, y = target*100, col = Method), lwd = 1.4) +
  facet_wrap(vars(pathogen)) + theme_bw() +
  ylab("Percent of disease targeted (%)") +
  xlab("Spatial cluster") +
 # theme(legend.position = "none") +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
       # legend.position = "none",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 16, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  scale_color_manual(values = c("royalblue2", "#FF4F00", "gray70")) + 
  ggtitle("b.")

plot_grid(compare_a, compare_b)


##############################################################
# Run diversity metric comparison for Bangladesh
#######################################################
# Load disease data. 
sth = read.csv(file = here("data/bangl/parasites/untouched", "bangl_analysis_parasite.csv")) %>%
  dplyr::select(dataid, clusterid, block, personid, tr, al, tt, giar, hw, sth) 
# sth = readr::read_csv
sth1 = read.csv(file = here("data/bangl/parasites/untouched", "bangl_analysis_parasite.csv")) %>%
  filter(tr == "Control" | tr == "Nutrition") 
head(sth1)
summary(sth$agey)

prot = read.csv(file = here("data/bangl/parasites/untouched", "washb-bangladesh-protozoa-public.csv")) %>%
  dplyr::select(dataid, clusterid, block, personid, tr, poseh , poscr) %>% #entamoeba, crypto
  mutate(poscr = replace(poscr, poscr == 9, NA)) %>% # 9 indicates the sample was missing 
  mutate(poseh = replace(poseh, poseh == 9, NA)) 

# Load measles data from final vacciantion records. 
vax = read.csv(file = here("data/bangl/vaccination/final", "washb_bangl_vax_records_midline.csv"))
head(vax)
# age in days is in anthropogegy  

# load this in to determine age, age needs to be at midline 
diarrhea = read.csv(file = here("data/bangl/diarrhea", "washb-bangladesh-diar.csv")) %>%
  filter(svy == 1) %>%
  dplyr::select(dataid, ageyrs)

# Load public id and gps data
public_ids = read.csv(file = here("data/bangl/public_ids", "public-ids.csv"))

gps_dat = read_dta(file = here("data/bangl/gps/untouched", "6. WASHB_Baseline_gps.dta")) %>%
  mutate(dataid = as.numeric(dataid)) %>% # had to add later? 
  left_join(public_ids, by = "dataid") %>% 
  dplyr::select(block, block_r, qgpslong, qgpslat) %>%
  group_by(block) %>%
  mutate(med_qgpslong = median(qgpslong), med_qgpslat = median(qgpslat)) %>%
  distinct(block, block_r, med_qgpslong, med_qgpslat)


#######################################################
# Join disease data and filter for control clusters. 
disease_dat = left_join(sth, prot, by = c("dataid", "clusterid", "block", "tr", "personid"))  %>%# %>%
  filter(tr == "Control" | tr == "Nutrition") 
block_n = length(unique(disease_dat$block))
head(disease_dat)

#######################################################
# Method 1: Population-level denominator and assess each STH independently. 

# Number of non-NA measurements assessed total by pathogen 
parasite_denom = disease_dat %>%
  pivot_longer(cols=c('al', 'tt', 'hw', 'giar', 'poscr', 'poseh', 'sth'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(pathogen) %>% # assess by overall prevalence in population
  summarise(non_na_count = sum(!is.na(presence))) %>%
  filter(pathogen != "sth" & pathogen != "giar" & pathogen != "poscr" & pathogen != "poseh")

# Number of parasite instances by block
parasite_numer = disease_dat %>%
  pivot_longer(cols=c('al', 'tt', 'hw', 'giar', 'poscr', 'poseh', 'sth'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(block, pathogen) %>%
  summarize(across(presence, ~sum(.x, na.rm = TRUE))) %>%
  filter(pathogen != "sth" & pathogen != "giar" & pathogen != "poscr" & pathogen != "poseh")
head(parasite_numer)

# Join together and calculate fraction prevalence and fraction if disease was distributed evenly amongst blocks. 
parasite_prevalence_block = left_join(parasite_denom, parasite_numer, by = c("pathogen")) %>%
  mutate(fraction = presence/non_na_count ) %>% # fraction = proportion of total pathogen prevalence accounted for by each block
  mutate(block_r = block) %>% dplyr::select(-block) %>% #specify that the blocks here are public
  ungroup() %>%
  group_by(pathogen) %>%
  mutate(pop_presence = sum(presence), even_prevalence = (pop_presence/non_na_count)/block_n) %>%
  dplyr::select(-non_na_count)
head(parasite_prevalence_block)
sum(parasite_prevalence_block$presence)

# Measles vaccination. 
measles_denom_block = vax %>%
  group_by(block_r) %>%
  summarise(non_na_count = sum(!is.na(measles))) %>%
  mutate(pathogen = "measles")
sum(measles_denom_block$non_na_count)

measles_numer_block = vax %>%
  group_by(block_r) %>%
  summarize(across(measles, ~sum(.x, na.rm = TRUE))) %>%
  mutate(presence = measles, pathogen = "measles") %>%
  dplyr::select(-measles)
sum(measles_numer_block$presence)

#Join together and calculate fraction unvaccinated and fraction unvaccinated if disease was distributed evenly amongst blocks. 
measles_prevalence_block = left_join(measles_denom_block, measles_numer_block, by = c("pathogen", "block_r")) %>%
  mutate(unvaccinated = non_na_count - presence) %>%
  dplyr::select(-presence) %>%
  mutate(total_assessed = sum(non_na_count), pop_presence = sum(unvaccinated)) %>%
  mutate(presence = unvaccinated) %>% #presence indicates number of unvaccinated individuals 
  mutate(fraction = unvaccinated/total_assessed, even_prevalence = (pop_presence /total_assessed)/block_n) %>%
  dplyr::select(pathogen, block_r, fraction, even_prevalence, pop_presence, presence)
head(measles_prevalence_block)

# Join measles vaccination and STH data. 
measles_path1 = rbind(parasite_prevalence_block, measles_prevalence_block ) %>%
  arrange(block_r) %>%
  ungroup()
head(measles_path1)

# plot cluster level distribution of prevalences 
histo = measles_path1 %>%
  mutate(pathogen = replace(pathogen, pathogen == "al", "Ascaris lumbricoides")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "hw", "Hookworm")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "measles", "Measles")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "tt", "Trichuris trichiura"))
ggplot(data = histo ) +
  geom_histogram(aes(x = fraction), fill = "blue4") +
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
  ggtitle("Bangladesh")

#### Pathogen-specific prevalence 
prev_overall = measles_path1 %>%
  mutate(prev = even_prevalence*90) %>%
  distinct(pathogen, prev)
head(prev_overall)

#######################################################
# Calculate Rao diversity at each block. 
product_results = list()

for (i in 1:block_n) {
  
  # Subset the data for the current block
  block_subset <- measles_path1 %>%
    mutate(fraction_mirrored = fraction) %>%
    ungroup() %>%
    filter(block_r == i) %>%
    dplyr::select(fraction_mirrored, fraction)
  
  # Get unique combinations of pathogen prevalences 
  fractions <- block_subset$fraction
  unique_combinations <- t(combn(fractions, 2, simplify = TRUE))  # Get pairs (no duplicates)
  
  # Compute products for each pair
  products <- apply(unique_combinations, 1, prod)
  
  products_df <- data.frame(
    block_r = i, 
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
  group_by(block_r) %>%
  summarize(across(product, ~sum(.x)))
head(block_prevalence_final)

# Calculate expected rao diversity for a given block.
block_subset = measles_path1 %>% 
  mutate(even_prevalence_mirrored = even_prevalence) %>%
  ungroup() %>% filter(block_r == 1) %>%
  dplyr::select(even_prevalence_mirrored, even_prevalence)


even_prevalences = block_subset$even_prevalence
unique_combinations <- t(combn(even_prevalences, 2, simplify = TRUE))  # Get pairs (no duplicates)

# Compute products for each pair
products <- apply(unique_combinations, 1, prod)
null_rao = sum(products)

# Compare observed rao to expected rao and add gps data. 
method1 = block_prevalence_final %>%
  mutate(rao = product/null_rao) %>%
  left_join(gps_dat, by = "block_r")%>%
  ungroup() %>%
  arrange(-rao)
summary(method1$rao)

# rao order of locations 
rao_order = print(method1$block_r)

# Calculate diversity metrics 
#####################################################
# Gini-simpson
gini_simpson = measles_path1 %>%
  group_by(block_r ) %>%
  mutate(block_presence = sum(presence)) %>%
  ungroup() %>%
  mutate(fraction_gini = presence/block_presence) %>%
  mutate(gs_component = fraction_gini^2) %>%
  group_by(block_r) %>%
  mutate(gs_entropy = sum(gs_component)) %>%
  ungroup() %>%
  mutate(gs_entropy_complement = 1 - gs_entropy ) %>%
  distinct(block_r, gs_entropy, gs_entropy_complement) %>%
  arrange(-gs_entropy_complement)
head(gini_simpson)
gs_order = print(gini_simpson$block_r) 

#####################################################
# Shannon
shannon = measles_path1 %>%
  group_by(block_r ) %>%
  mutate(block_presence = sum(presence)) %>%
  ungroup() %>%
  mutate(fraction_shannon = presence/block_presence) %>%
  mutate(fraction_shannon = fraction_shannon + .0000000000001) %>%
  mutate(shannon_component = fraction_shannon*log2(fraction_shannon)) %>%
  group_by(block_r) %>%
  mutate(shannon_entropy = - sum(shannon_component)) %>%
  distinct(block_r, shannon_entropy) %>%
  arrange(-shannon_entropy)
head(shannon)
shannon_order = print(shannon$block_r) 
print(shannon_order)

# Calculate the percent of disease targeted for each pathogen in Cambodia 
rao_ordered <- measles_path1 %>%
  group_by(pathogen) %>%
  mutate(customorder = factor(block_r, levels = rao_order)) %>%
  arrange(customorder) %>%  # arrange by the new factor
  mutate(
    cum = cumsum(fraction),
    cum_presence = cumsum(presence),
    total_fraction = sum(fraction),
    target = cum / total_fraction,
    block_label = seq(1, 90, 1)
  ) %>%
  ungroup() %>%
  mutate(Method = "Rao") %>%
  group_by(block_r) %>%
  mutate(people_treated = sum(cum_presence))

gs_ordered = measles_path1 %>%
  group_by(pathogen) %>%
  mutate(customorder = factor(block_r, levels = gs_order)) %>%
  arrange(customorder) %>%  # arrange by the new factor
  mutate(
    cum = cumsum(fraction),
    cum_presence = cumsum(presence),
    total_fraction = sum(fraction),
    target = cum / total_fraction,
    block_label = seq(1, 90, 1)
  ) %>%
  ungroup() %>%
  mutate(Method = "Gini-Simpson") %>%
  group_by(block_r) %>%
  mutate(people_treated = sum(cum_presence))

shannon_ordered = measles_path1 %>%
  group_by(pathogen) %>%
  mutate(customorder = factor(block_r, levels = shannon_order)) %>%
  arrange(customorder) %>%  # arrange by the new factor
  mutate(
    cum = cumsum(fraction),
    cum_presence = cumsum(presence),
    total_fraction = sum(fraction),
    target = cum / total_fraction,
    block_label = seq(1, 90, 1)
  ) %>%
  ungroup() %>%
  mutate(Method = "Shannon") %>%
  group_by(block_r) %>%
  mutate(people_treated = sum(cum_presence))

all_methods = rbind(rao_ordered, gs_ordered, shannon_ordered) %>%
  mutate(pathogen = replace(pathogen, pathogen == "al", "Ascaris lumbricoides")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "hw", "Hookworm")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "measles", "Measles")) %>%
  mutate(pathogen = replace(pathogen, pathogen == "tt", "Trichurus trichuria")) %>%
  mutate(pathogen = factor(pathogen, levels = c("Ascaris lumbricoides", "Hookworm", 
                                                "Measles", "Trichurus trichuria")))

compare_a = ggplot(data = all_methods ) + 
  geom_line(aes(x = block_label, y = people_treated, col = Method), lwd = 1.9) +
  theme_bw() +
  xlab("Spatial cluster") +
  ylab("Number of people treated") +
  ggtitle("a.") +
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
  scale_color_manual(values = c("royalblue2", "#FF4F00", "gray70")) 


compare_b = ggplot(data = all_methods ) + 
  geom_line(aes(x = block_label, y = target*100, col = Method), lwd = 1.4) +
  facet_wrap(vars(pathogen)) + theme_bw() +
  ylab("Percent of disease targeted (%)") +
  xlab("Spatial cluster") +
  # theme(legend.position = "none") +
  theme(axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        axis.line = element_blank(),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
        plot.title.position = "plot",
        plot.subtitle = element_text(colour = "black", size = 12.5),
        # legend.position = "none",
        legend.key.width = unit(0.5, "cm"),
        legend.text = element_text(size = 13, color = "black"),
        legend.title = element_text(size = 14, color = "black"),
        strip.text = element_text(colour = "black", size = 16, hjust = 0),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  scale_color_manual(values = c("royalblue2", "#FF4F00", "gray70")) + 
  ggtitle("b.")

plot_grid(compare_a, compare_b)




