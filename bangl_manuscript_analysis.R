########################################################
#~/Library/CloudStorage/Box-Box
#finder --> settings 
#_______________________________________________________
# @Organization --  UCSF 
# @Project -- WASH Benefits Bangladesh multipathogen burden
# @Author -- Samantha Bents, sjbents@stanford.edu, November 30, 2024 
# @Description -- Assess integrated intervention for measles vaccination and soil-transmitted helminths. 
# (hookworm, trichuris, ascaris) using rao diversity 
#_______________________________________________________
########################################################
library(here)
source(here("0-config.R"))
source(here("0-functions.R"))

#######################################################
# Load disease data. 
sth = readr::read_csv(file = here("data/bangl/parasites/untouched", "bangl_analysis_parasite.csv"))%>%
  dplyr::select(dataid, clusterid, block, personid, tr, al, tt, giar, hw, sth) 

prot = readr::read_csv(file = here("data/bangl/parasites/untouched", "washb-bangladesh-protozoa-public.csv")) %>%
  dplyr::select(dataid, clusterid, block, personid, tr, poseh , poscr) %>% #entamoeba, crypto
  mutate(poscr = replace(poscr, poscr == 9, NA)) %>% # 9 indicates the sample was missing 
  mutate(poseh = replace(poseh, poseh == 9, NA)) 

# Load measles data from final vacciantion records. 
vax = readr::read_csv(file = here("data/bangl/vaccination/final", "washb_bangl_vax_records_midline.csv"))

# Load public id and gps data
public_ids = readr::read_csv(file = here("data/bangl/public_ids", "public-ids.csv"))
gps_dat = read_dta(file = here("data/bangl/gps/untouched", "6. WASHB_Baseline_gps.dta")) %>%
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

# Measles vaccination. 
measles_denom_block = vax %>%
  group_by(block_r) %>%
  summarise(non_na_count = sum(!is.na(measles))) %>%
  mutate(pathogen = "measles")

measles_numer_block = vax %>%
  group_by(block_r) %>%
  summarize(across(measles, ~sum(.x, na.rm = TRUE))) %>%
  mutate(presence = measles, pathogen = "measles") %>%
  dplyr::select(-measles)

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
  arrange(block_r)

head(measles_path1)

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
  left_join(gps_dat, by = "block_r")
summary(method1$rao)

#######################################################
#######################################################
# Method 2: Population-level denominator and assess STH overall
parasite_numer = disease_dat %>%
  pivot_longer(cols=c('al', 'tt', 'hw', 'giar', 'poscr', 'poseh', 'sth'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(block, pathogen) %>%
  summarize(across(presence, ~sum(.x, na.rm = TRUE))) %>%
  filter(pathogen == "sth")

# Number of non-NA measurements assessed total by pathogen 
parasite_denom = disease_dat %>%
  pivot_longer(cols=c('al', 'tt', 'hw', 'giar', 'poscr', 'poseh', 'sth'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(pathogen) %>% # assess by overall prevalence in population
  summarise(non_na_count = sum(!is.na(presence))) %>%
  filter(pathogen == "sth")

# Join together and calculate fraction prevalence and fraction if disease was distributed evenly amongst blocks. 
parasite_prevalence_block = left_join(parasite_denom, parasite_numer, by = c("pathogen")) %>%
  mutate(fraction = presence/non_na_count ) %>% # fraction = proportion of total pathogen prevalence accounted for by each block
  mutate(block_r = block) %>% dplyr::select(-block) %>% #specify that the blocks here are public
  ungroup() %>%
  group_by(pathogen) %>%
  mutate(pop_presence = sum(presence), even_prevalence = (pop_presence/non_na_count)/block_n) %>%
  dplyr::select(-non_na_count)
head(parasite_prevalence_block)


# Join measles vaccination and STH data. 
measles_path2 = rbind(parasite_prevalence_block, measles_prevalence_block ) %>%
  arrange(block_r)

#######################################################
# Calculate Rao diversity at each block. 
product_results = list()

for (i in 1:block_n) {
  
  # Subset the data for the current block
  block_subset <- measles_path2 %>%
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
block_subset = measles_path2 %>% 
  mutate(even_prevalence_mirrored = even_prevalence) %>%
  ungroup() %>% filter(block_r == 1) %>%
  dplyr::select(even_prevalence_mirrored, even_prevalence)

even_prevalences = block_subset$even_prevalence
unique_combinations <- t(combn(even_prevalences, 2, simplify = TRUE))  # Get pairs (no duplicates)

# Compute products for each pair
products <- apply(unique_combinations, 1, prod)
null_rao = sum(products)

# Compare observed rao to expected rao and add gps data. 
method2 = block_prevalence_final %>%
  mutate(rao = product/null_rao) %>%
  left_join(gps_dat, by = "block_r")
summary(method2$rao)


#######################################################
#######################################################
# Method 3: Pathogens-specific denominator and assess each STH independently. 

# Number of non-NA measurements assessed total by pathogen 
parasite_denom = disease_dat  %>%
  pivot_longer(cols=c('al', 'tt', 'hw', 'giar', 'poscr', 'poseh', 'sth'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(pathogen) %>% # assess by overall prevalence in population
  summarise(non_na_count = sum(presence, na.rm = TRUE)) %>%
  filter(pathogen != "sth" & pathogen != "giar" & pathogen != "poscr" & pathogen != "poseh")
head(parasite_denom)

# Number of total parasite instances by block
parasite_numer = disease_dat %>%
 # filter(tr == "Control") %>%
  pivot_longer(cols=c('al', 'tt', 'hw', 'giar', 'poscr', 'poseh', 'sth'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(block, pathogen) %>%
  summarize(across(presence, ~sum(.x, na.rm = TRUE))) %>%
  filter(pathogen != "sth" & pathogen != "giar" & pathogen != "poscr" & pathogen != "poseh")

# Join together and calculate fraction unvaccinated and fraction unvaccinated if disease was distributed evenly amongst blocks. 
parasite_prevalence_block = left_join(parasite_denom, parasite_numer, by = c("pathogen")) %>%
  mutate(fraction = presence/non_na_count ) %>% # fraction = proportion of total prevalence accounted for by each block
  mutate(block_r = block) %>% dplyr::select(-block) %>%
  ungroup() %>%
  group_by(pathogen) %>%
  mutate(pop_presence = sum(presence), even_prevalence = 1/block_n) %>%
  dplyr::select(-non_na_count) 

# Measles unvaccinated. 
measles_denom_block = vax %>%
  group_by(block_r) %>%
  summarise(non_na_count = sum(!is.na(measles))) %>%
  mutate(pathogen = "measles")

measles_numer_block = vax %>%
  group_by(block_r) %>%
  summarize(across(measles, ~sum(.x, na.rm = TRUE))) %>%
  mutate(presence = measles, pathogen = "measles") %>%
  dplyr::select(-measles)

measles_prevalence_block = left_join(measles_denom_block, measles_numer_block, by = c("pathogen", "block_r")) %>%
  mutate(unvaccinated = non_na_count - presence) %>%
  dplyr::select(-presence) %>%
  mutate(total_assessed = sum(non_na_count), pop_presence = sum(unvaccinated)) %>%
  mutate(presence = unvaccinated) %>%
  mutate(fraction = unvaccinated/pop_presence, even_prevalence = 1/block_n) %>%
  dplyr::select(pathogen, block_r, fraction, even_prevalence, pop_presence, presence)
head(measles_prevalence_block)

### Join measles and STH 
measles_path3 = rbind(parasite_prevalence_block, measles_prevalence_block ) %>%
  arrange(block_r)
head(measles_path3)

#######################################################
# Calculate Rao diversity at each block. 
product_results = list()

for (i in 1:block_n) {
  
  # Subset the data for the current block
  block_subset <- measles_path3 %>%
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
block_subset = measles_path3 %>% 
  mutate(even_prevalence_mirrored = even_prevalence) %>%
  ungroup() %>% filter(block_r == 1) %>%
  dplyr::select(even_prevalence_mirrored, even_prevalence)


even_prevalences = block_subset$even_prevalence
unique_combinations <- t(combn(even_prevalences, 2, simplify = TRUE))  # Get pairs (no duplicates)

# Compute products for each pair
products <- apply(unique_combinations, 1, prod)
null_rao = sum(products)

# Compare observed rao to expected rao and add gps data. 
method3 = block_prevalence_final %>%
  mutate(rao = product/null_rao) %>%
  left_join(gps_dat, by = "block_r")
summary(method3$rao)

#######################################################
#######################################################
# Method 4: Pathogens-specific denominator and assess STH overall. 
parasite_denom = disease_dat %>%
  pivot_longer(cols=c('al', 'tt', 'hw', 'giar', 'poscr', 'poseh', 'sth'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(pathogen) %>% # assess by overall prevalence in population
  summarise(non_na_count = sum(presence, na.rm = TRUE)) %>%
  filter(pathogen == "sth")

# Number of total parasite instances by block
parasite_numer = disease_dat %>%
  pivot_longer(cols=c('al', 'tt', 'hw', 'giar', 'poscr', 'poseh', 'sth'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(block, pathogen) %>%
  summarize(across(presence, ~sum(.x, na.rm = TRUE))) %>%
  filter(pathogen == "sth")

# Join together and calculate fraction unvaccinated and fraction unvaccinated if disease was distributed evenly amongst blocks. 
parasite_prevalence_block = left_join(parasite_denom, parasite_numer, by = c("pathogen")) %>%
  mutate(fraction = presence/non_na_count ) %>% # fraction = proportion of total prevalence accounted for by each block
  mutate(block_r = block) %>% dplyr::select(-block) %>%
  ungroup() %>%
  group_by(pathogen) %>%
  mutate(pop_presence = sum(presence), even_prevalence = 1/block_n) %>%
  dplyr::select(-non_na_count) 

### Join measles and STH 
measles_path4 = rbind(parasite_prevalence_block, measles_prevalence_block ) %>%
  arrange(block_r)

#######################################################
# Calculate Rao diversity at each block. 
product_results = list()

for (i in 1:block_n) {
  
  # Subset the data for the current block
  block_subset <- measles_path4 %>%
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
block_subset = measles_path4 %>% 
  mutate(even_prevalence_mirrored = even_prevalence) %>%
  ungroup() %>% filter(block_r == 1) %>%
  dplyr::select(even_prevalence_mirrored, even_prevalence)

even_prevalences = block_subset$even_prevalence
unique_combinations <- t(combn(even_prevalences, 2, simplify = TRUE))  # Get pairs (no duplicates)

# Compute products for each pair
products <- apply(unique_combinations, 1, prod)
null_rao = sum(products)

# Compare observed rao to expected rao and add gps data. 
method4 = block_prevalence_final %>%
  mutate(rao = product/null_rao) %>%
  left_join(gps_dat, by = "block_r")
summary(method4$rao)

#######################################################
# Map methods visually.
plot_list <- list()
methods = list(method1 = method1, method2 = method2, method3 = method3, method4 = method4)

for(i in seq_along(methods)){
  
  method_df = methods[[i]]
  
  #Retriece map
  admin_b_study <- get_map("district")
  
  dbgps_block = method_df %>% 
    mutate(lat = med_qgpslat, lon = med_qgpslong) %>%
    dplyr::select(block_r, lon, lat) %>%
    dplyr::distinct()  %>% mutate(block = block_r) %>%
    dplyr::select(-block_r)
  
  dbgps = method_df  %>%
    mutate(lat = med_qgpslat, lon = med_qgpslong) %>%
    ungroup() %>%
    dplyr::select(block_r, rao, lat, lon) %>%
    distinct() %>% mutate(block = block_r)  %>%
    dplyr::select(-block_r)
  
  # Create bonding box from the shapefile. 
  b_boxgrid <- st_make_grid(admin_b_study,
                            n=c(100,100), #this is the number of rows and col to split the shapefile into. 
                            what = "centers" #this means that the function will return the center points of each grid, alternatively the shape file would be returned. 
  )
  
  # Transform cluster centroids into shape file. 
  dbgps_sf_utm = dbgps_block %>%
    ungroup() %>%
    st_as_sf(coords = c("lon", "lat"), crs=4326) %>% # crs=4326 specifies that lat and long are coordinates
    dplyr::select(geometry) %>%  
    st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km") # transforms the spatial data from WGS84 coordinates (longitude and latitude) to UTM coordinates
  
  # Identify a buffer of 10 km (based on Arnold et al. 2024).
  cl_buff <- st_buffer(dbgps_sf_utm, dist=10) %>%  # creates a buffer of 10km around each center 
    summarise(geometry = st_union(geometry)) %>%  # identify their union, then convert to a single polygon that merges overlapping areas of the buffers into one continuous polygon.
    st_cast("POLYGON") %>% # set the CRS back to WGS84 lat and long coordinates 
    st_transform(crs=4326) # this transforms the geometry back to the WGS84 coordinate reference system (EPSG:4326), which uses degrees of latitude and longitude. 
  
  # Fit the model 
  b_fit_rao <- spaMM::fitme(rao ~  lat + lon + Matern(1|lat + lon), # Matern spatial random effect which is fixed 
                            data = dbgps,
                            family = gaussian(link = "identity")) 
  # if you do not specify Matern parameters, spaMM defaults to v = .5 (exponential decay of smoothness)
  # ϕ is the scale parameter, which controls the distance over which the correlation decays
  # Bessel function ensures covariance declines with distance between points, gamma function ensures normalization
  print(summary(b_fit_rao))
  
  #outputs nu and rho parameters, nu indiciates spatial smoothness, rho indicates spatial influence for the correlations modeled by the Matern function.
  plot_materncor <- plot_matern_corr(b_fit_rao) 
  print(plot_materncor)
  
  # this function predicts the response variable = pred and the variance = pred_var over the grid delineated earlier for each shapefile center. 
  b_preds_rao <- get_grid_preds(input_grid=b_boxgrid, 
                                spamm_model_fit = b_fit_rao)
  b_preds_rao_coords <- st_coordinates(b_preds_rao)
  head(b_preds_rao_coords ) # converts shape file centers back to coordinates.
  
  b_preds_rao_raster <- points_to_raster( #for each of the 100 x 100 gridlines, 
    x = b_preds_rao_coords[,1],
    y = b_preds_rao_coords[,2],
    z = b_preds_rao$pred, # was multiplied by 100 before 
    mask1 = admin_b_study,
    mask2 = cl_buff,
    crop1 = admin_b_study
  )
  
  b_preds_rao_raster2 <- raster_to_tibble(b_preds_rao_raster)
  print(head(b_preds_rao_raster2 ))
  
  map_rao <- ggplot(data= b_preds_rao_raster2 %>% filter(value > 1)) + 
    geom_tile(aes(x=x,y=y,fill=value), na.rm=TRUE) +
    ggtitle("Method", i) +
    coord_sf(crs = 4326) +
    scale_fill_viridis(option = "D", direction = 1,
                       limits = c(0, max(b_preds_rao_raster2$value)),
                       na.value = NA,
                       breaks = seq(0, 3,by= .5),
                       guide=guide_colorbar(title="Rao diversity", 
                                            title.position = "top",
                                            title.hjust = 0.5,
                                            direction="horizontal",
                                            label.position = "bottom",
                                            label.hjust = 0.5,
                                            barheight=unit(10,"pt"),
                                            barwidth=unit(200,"pt"),
                                            ticks.colour = "white",
                                            ticks.linewidth = 1,
                                            frame.colour = "white"
                       )
    ) +
    labs(x="Latitude",y="Longitude") +
    xlim(c(89.65, 91)) +
    ylim(c(23.7, 25)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.tag = element_text(face="bold",size=16)
    )  
  
  print(map_rao)
  
  plot_list[[i]] = map_rao
} 

method_plot <- plot_grid(plotlist = plot_list, ncol = 2)
method_plot

#######################################################


# Compare method efficiency. 
#######################################################
pt = .5
block_50 <- list()
efficiency_list <- list()
#methods = list(method1 = method1, method2 = method2, method3 = method3, method4 = method4)
methods = list(method1 = method1)
head(method1)


for(i in seq_along(methods)){
  
  method_df = methods[[i]]
  compare_strategy = left_join(measles_path1, method_df, by = "block_r")
  
  # Measles strategy 
  measles_motivated = compare_strategy %>% 
    filter(pathogen == "measles") %>%
    arrange(-fraction) %>% # Arrange from highest prevalence block to lowest 
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  # order_measles = print(measles_motivated$block_r)
  
  measles_50 <- which(measles_motivated$above == 2)[1]
  
  rao_measles_motivated = compare_strategy %>% 
    filter(pathogen == "measles") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  print(rao_measles_motivated$block)
  
  measles_rao_50 <- which(rao_measles_motivated$above == 2)[1]
  
  measles = ggplot(data = measles_motivated, aes(x = block_label/block_n*90, y = target*100)) +
    geom_line(aes(col = "Measles"), cex = 2) +
    geom_line(data = rao_measles_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Rao"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent vaccinated (%)") +
    scale_color_manual(values = c("Measles" = "black", "Rao" = "#FF4F00")) +
    guides(color =guide_legend(title="Block prioritization strategy")) +
    geom_abline(slope=100/90, intercept = 0, lwd =1, lty = 2,  col = "gray72") +
    ggtitle("      Measles")+
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90)) +
    theme(axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "black", size = 12.5),
         # legend.position = "right",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 11, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          strip.text = element_text(colour = "black", size = 10, hjust = 0),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill=NA)) 
  
  # Calculate AUC for measles  
  # auc_measles <- trapz(measles_motivated$block_label/block_n, measles_motivated$target)
  #  auc_rao <- trapz(rao_measles_motivated$block_label/block_n, rao_measles_motivated$target)
  #  measles_auc = list(auc_measles,  auc_rao)
  #  print(measles_auc)
  
  # Ascaris strategy
  ascaris_motivated = compare_strategy %>% 
    filter(pathogen == "al") %>%
    arrange(-fraction) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1))%>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  print(head(ascaris_motivated))
  
  ascaris_50 <- which(ascaris_motivated$above == 2)[1]
  
  rao_ascaris_motivated = compare_strategy %>% 
    filter(pathogen == "al") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1))%>%
    mutate(above = ifelse(target < pt, 1, 2))
  print(rao_ascaris_motivated$block)
  
  ascaris_rao_50 <- which(rao_ascaris_motivated$above == 2)[1]
  
  ascaris_measles_motivated = compare_strategy %>% 
    filter(pathogen == "al") %>%
    mutate(custom_order = factor(block_r), levels = order_measles) %>%
    arrange(-order_measles) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1))%>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  ascaris_measles_50 <- which(ascaris_measles_motivated$above == 2)[1]
  
  ascaris = ggplot(data = ascaris_motivated, aes(x = block_label/block_n*90, y = target*100)) +
    geom_line(aes(col = "Ascaris"), cex = 2) +
    geom_line(data = rao_ascaris_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = ascaris_measles_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Measles"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("Ascaris" = "black", "Rao" = "#FF4F00", "Measles" = "#0094C6")) +
    guides(color =guide_legend(title="Block prioritization strategy")) +
    geom_abline(slope=100/90, intercept = 0, lwd =1, lty = 2,  col = "gray72") +
    ggtitle("      Ascaris lumbricoides")+
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90)) +
    theme(axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "black", size = 12.5),
     #     legend.position = "right",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 11, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          strip.text = element_text(colour = "black", size = 10, hjust = 0),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill=NA)) 
  
  # Calculate AUC for ascaris
  # auc_ascaris <- trapz(ascaris_motivated$block_label/block_n, ascaris_motivated$target)
  #  auc_rao <- trapz(rao_ascaris_motivated$block_label/block_n, rao_ascaris_motivated$target)
  #  auc_ascaris_measles <- trapz(  ascaris_measles_motivated$block_label/block_n,   ascaris_measles_motivated$target)
  #  ascaris_auc = list(auc_ascaris, auc_rao, auc_ascaris_measles)
  #  print(ascaris_auc)
  
  # T. trichuris strategy
  trichurus_motivated = compare_strategy %>% 
    filter(pathogen == "tt") %>%
    arrange(-fraction) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1))%>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  trichurus_50 <- which(trichurus_motivated$above == 2)[1]
  
  rao_trichurus_motivated = compare_strategy %>% 
    filter(pathogen == "tt") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1))%>%
    mutate(above = ifelse(target < pt, 1, 2))
  print( rao_trichurus_motivated$block)
  
  trichurus_rao_50 <- which(rao_trichurus_motivated$above == 2)[1]
  
  trichurus_measles_motivated = compare_strategy %>% 
    filter(pathogen == "tt") %>%
    mutate(custom_order = factor(block_r), levels = order_measles) %>%
    arrange(-order_measles) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1))%>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  trichurus_measles_50 <- which(trichurus_measles_motivated$above == 2)[1]
  
  trichurus = ggplot(data = trichurus_motivated, aes(x = block_label/block_n*90, y = target*100)) +
    geom_line(aes(col = "T. trichuris"), cex = 2) +
    geom_line(data = rao_trichurus_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = trichurus_measles_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Measles"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("T. trichuris" = "black", "Rao" = "#FF4F00", "Measles" = "#0094C6")) +
    guides(color =guide_legend(title="Block prioritization strategy")) +
    geom_abline(slope=100/90, intercept = 0, lwd =1, lty = 2,  col = "gray72") +
    ggtitle("      Trichuris trichiura")+
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90)) +
    theme(axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "black", size = 12.5),
     #     legend.position = "right",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 11, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          strip.text = element_text(colour = "black", size = 10, hjust = 0),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill=NA)) 
  
  # Calculate AUC for trichurs 
  # auc_trichurus <- trapz(trichurus_motivated$block_label/block_n, trichurus_motivated$target)
  #  auc_rao <- trapz(rao_trichurus_motivated$block_label/block_n, rao_trichurus_motivated$target)
  #  auc_trichurus_measles <- trapz(trichurus_measles_motivated$block_label/block_n,  trichurus_measles_motivated$target)
  #  trichurus_auc = list(auc_trichurus, auc_rao, auc_trichurus_measles)
  #  print(trichurus_auc)
  
  # Hookworm strategy
  hookworm_motivated = compare_strategy %>% 
    filter(pathogen == "hw") %>%
    arrange(-fraction) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1))%>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  hookworm_50 <- which(hookworm_motivated$above == 2)[1]
  
  rao_hookworm_motivated = compare_strategy %>% 
    filter(pathogen == "hw") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1))%>%
    mutate(above = ifelse(target < pt, 1, 2))
  print(rao_hookworm_motivated$block)
  
  hookworm_rao_50 <- which(rao_hookworm_motivated$above == 2)[1]
  
  hookworm_measles_motivated = compare_strategy %>% 
    filter(pathogen == "hw") %>%
    mutate(custom_order = factor(block_r), levels = order_measles) %>%
    arrange(-order_measles) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1))%>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  hookworm_measles_50 <- which(hookworm_measles_motivated$above == 2)[1]
  
  hookworm = ggplot(data = hookworm_motivated, aes(x = block_label/block_n*90, y = target*100)) +
    geom_line(aes(col = "Hookworm"), cex = 2) +
    geom_line(data = rao_hookworm_motivated , aes(x = block_label/block_n*90, y = target*100, col = "Rao"), cex = 2) + 
    geom_line(data = hookworm_measles_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Measles"), cex = 2) + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("Hookworm" = "black", "Rao" = "#FF4F00", "Measles" = "#0094C6")) +
    guides(color =guide_legend(title="Block prioritization strategy")) +
    geom_abline(slope=100/90, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    ggtitle("      Hookworm")+
    theme(legend.position = "none") +
    scale_x_continuous(breaks = c(0, 15, 30, 45, 60, 75, 90)) +
    theme(axis.text = element_text(size = 11, color = "black"),
          axis.title = element_text(size = 12, color = "black"),
          axis.line = element_blank(),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(colour = "black", size = 13.5, face = "bold"),
          plot.title.position = "plot",
          plot.subtitle = element_text(colour = "black", size = 12.5),
       #   legend.position = "right",
          legend.key.width = unit(0.5, "cm"),
          legend.text = element_text(size = 11, color = "black"),
          legend.title = element_text(size = 12, color = "black"),
          strip.text = element_text(colour = "black", size = 10, hjust = 0),
          strip.background = element_rect(colour="white", fill="white"),
          panel.border = element_rect(colour = "black", fill=NA)) 
  
  # Calculate AUC for hw
  #auc_hookworm <- trapz(hookworm_motivated$block_label/block_n, hookworm_motivated$target)
  #auc_rao <- trapz(rao_hookworm_motivated$block_label/block_n, rao_hookworm_motivated$target)
  #auc_hookworm_measles <- trapz(hookworm_measles_motivated$block_label/block_n,  hookworm_measles_motivated$target)
  #hookworm_auc = list(auc_hookworm, auc_rao, auc_hookworm_measles)
  #print(hookworm_auc)
  
  
  four_path_plot = plot_grid(measles, ascaris, trichurus, hookworm)
  efficiency_list[[i]] = four_path_plot 
  
  
  block_50[[i]] <- list(measles_50 = measles_50, measles_rao_50 = measles_rao_50,
                        ascaris_50 = ascaris_50, ascaris_rao_50 = ascaris_rao_50,
                        ascaris_measles_50 = ascaris_measles_50, 
                        hookworm_50 = hookworm_50, hookworm_rao_50 = hookworm_rao_50,
                        hookworm_measles_50 = hookworm_measles_50,
                        trichurus_50 = trichurus_50, trichurus_rao_50 = trichurus_rao_50,
                        trichurus_measles_50 = trichurus_measles_50)
}

# Plot efficiency of single-pathogen strategies vs rao strategies. 
method1_efficiency <- plot_grid(efficiency_list = efficiency_list[[1]], ncol = 2)
method1_efficiency 
method2_efficiency <- plot_grid(efficiency_list = efficiency_list[[2]], ncol = 2)
method2_efficiency
method3_efficiency <- plot_grid(efficiency_list = efficiency_list[[3]], ncol = 2)
method3_efficiency
method4_efficiency <- plot_grid(efficiency_list = efficiency_list[[4]], ncol = 2)
method4_efficiency


print(block_50[[1]])
print(block_50[[2]])

method1_efficiency <- efficiency_list[[1]]
method1_efficiency 

#method2_efficiency <- plot_grid(efficiency_list = efficiency_list[[2]], ncol = 2)
#method2_efficiency


# make legend for strategy 
legend2 = ggplot(data = parasite_prevalence_block) +
  geom_line(aes(x = block_r, y = presence, col = "Measles-motivated"), cex = 2) + 
  geom_line(aes(x = block_r, y = fraction, col = "Single-pathogen motivated"), cex = 2) + 
  geom_line(aes(x = block_r, y = even_prevalence, col = "Rao-motivated"), cex = 2) + 
  geom_line(aes(x = block_r, y = pop_presence* even_prevalence, col = "Single-pathogen motivated"), cex = 2) + 
  scale_color_manual(values = c( "Rao-motivated" = "#FF4F00", 
                                 "Measles-motivated" = "#0094C6", "Single-pathogen motivated" = "black")) +
  guides(color =guide_legend(title="Strategy"))  +
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 14, color = "black"))+ 
  geom_blank() 
legend2
legend_only2 <- gtable::gtable_filter(ggplotGrob(legend2), "guide-box")


fig2_bangl_eff <- cowplot::plot_grid(
  method1_efficiency,   # Main plot
  legend_only2,           # Extracted legend
  ncol = 1,              # Layout with one column
  rel_heights = c(1, 0.1)  # Adjust relative heights to make space for the legend
)

fig2_bangl_eff


#### legend 
#head(parasite_prevalence_block)

#legend2 = ggplot(data = parasite_prevalence_block) +
#  geom_line(aes(x = block_r, y = presence, col = "Measles-motivated"), cex = 2) + 
#  geom_line(aes(x = block_r, y = fraction, col = "Single-pathogen motivated"), cex = 2) + 
#  geom_line(aes(x = block_r, y = even_prevalence, col = "Rao-motivated"), cex = 2) + 
#  geom_line(aes(x = block_r, y = pop_presence* even_prevalence, col = "Single-pathogen motivated"), cex = 2) + 
#  scale_color_manual(values = c( "Rao-motivated" = "#FF4F00", 
#                                "Measles-motivated" = "#0094C6", "Single-pathogen motivated" = "black")) +
#  guides(color =guide_legend(title="Strategy"))  +
#  theme(legend.position = "right") +
#  geom_blank() 
#legend2

#leg = get_legend(legend2)
#grid.newpage()
#legend_only = grid.draw(leg)

#plot_grid(method1_efficiency , legend_only, rel_widths = c(1, .2), labels = "b")
#legend_only = grid.draw(leg)





# simulations (compared to method 1 )
pt = .5
sim_50 <- list()

for(i in 1:1000) {
  
  method_df = method1
  compare_strategy = left_join(measles_path1, method_df, by = "block_r")
  
  # rao motivated 
  measles_rao_motivated = compare_strategy %>% 
    filter(pathogen == "measles") %>%
    mutate(random = sample(1:90, 90, replace = FALSE)) %>%
    arrange(-random) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  rao_measles_50_sim <- which(measles_rao_motivated$above == 2)[1]
  sim_50[[i]] <- rao_measles_50_sim
  
}

measles_vector <- unlist(sim_50)
measles_sim <- data.frame(measles_vector)

# simulations

sim_50 <- list()

for(i in 1:1000) {
  
  method_df = method1
  compare_strategy = left_join(measles_path1, method_df, by = "block_r")
  
  
  # rao motivated 
  rao_ascaris_motivated = compare_strategy %>% 
    filter(pathogen == "al") %>%
    mutate(random = sample(1:90, 90, replace = FALSE)) %>%
    arrange(-random) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  rao_ascaris_50_sim <- which(rao_ascaris_motivated$above == 2)[1]
  sim_50[[i]] <- rao_ascaris_50_sim 
  
}

ascaris_vector <- unlist(sim_50)
ascaris_sim <- data.frame(ascaris_vector)


# tt
sim_50 <- list()

for(i in 1:1000) {
  
  method_df = method1
  compare_strategy = left_join(measles_path1, method_df, by = "block_r")
  
  
  #rao motivated
  rao_trichurus_motivated = compare_strategy %>% 
    filter(pathogen == "tt") %>%
    mutate(random = sample(1:90, 90, replace = FALSE)) %>%
    arrange(-random) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  rao_trichurus_50_sim <- which(rao_trichurus_motivated$above == 2)[1]
  sim_50[[i]] <- rao_trichurus_50_sim 
  
}

trichurus_vector <- unlist(sim_50)
trichurus_sim <- data.frame(trichurus_vector)


sim_50 <- list()

for(i in 1:1000) {
  
  method_df = method1
  compare_strategy = left_join(measles_path1, method_df, by = "block_r")
  
  
  #rao motivated
  rao_hookworm_motivated = compare_strategy %>% 
    filter(pathogen == "hw") %>%
    mutate(random = sample(1:90, 90, replace = FALSE)) %>%
    arrange(-random) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  rao_hookworm_50_sim <- which(rao_hookworm_motivated$above == 2)[1]
  sim_50[[i]] <- rao_hookworm_50_sim 
  
}

hookworm_vector <- unlist(sim_50)
hookworm_sim <- data.frame(hookworm_vector)


#### boxplot with them all
print(block_50[[1]])

me_boxplot = measles_sim %>%
  mutate(pathogen = "Measles") %>%
  mutate(Rao = 42) %>%
  mutate(Optimal = 38)
head(me_boxplot)

as_boxplot = ascaris_sim %>%
  mutate(pathogen = "Ascaris lumbricoides") %>%
  mutate(Rao = 36) %>%
  mutate(Optimal = 31) %>%
  mutate(Measles = 46) 
head(as_boxplot )

tt_boxplot = trichurus_sim %>%
  mutate(pathogen = "Trichuris trichiura") %>%
  mutate(Rao = 19) %>%
  mutate(Optimal = 14) %>%
  mutate(Measles = 40) 
head(tt_boxplot)

hw_boxplot = hookworm_sim %>%
  mutate(pathogen = "Hookworm") %>%
  mutate(Rao = 23) %>%
  mutate(Optimal = 19) %>%
  mutate(Measles = 45) 

target_goal = ggplot(data = me_boxplot) + 
  geom_violin(aes(x = pathogen, y = measles_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
  # geom_boxplot(aes(x = pathogen, y = malaria_v_vector), fill = "gray95", width=.6) + 
  geom_point(aes(x = pathogen, y = Rao),  col = "#FF4F00", cex = 6)  + #.5, rao  
  geom_point(aes(x = pathogen, y = Optimal), col = "black",  cex = 6)  + 
  theme_bw() +
  ylab("Clusters to target 50% of disease") +
  geom_violin(data = as_boxplot, aes(x = pathogen, y = malaria_f_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
  # geom_boxplot(data = mf_boxplot, aes(x = pathogen, y = malaria_f_vector), fill = "gray95", width=.6) + 
  geom_point(data =as_boxplot, aes(x = pathogen, y = Rao), col = "#FF4F00",  cex = 6)  + #.5, rao  
  geom_point(data = as_boxplot, aes(x = pathogen, y = Optimal), col = "black",  cex = 6)  + 
  geom_point(data = as_boxplot, aes(x = pathogen, y = Measles), col = "#0094C6",  cex = 6)  + 
  geom_violin(data = tt_boxplot, aes(x = pathogen, y = strong_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
  # geom_boxplot(data = strong_boxplot, aes(x = pathogen, y = strong_vector), fill = "gray95", width=.6) + 
  geom_point(data = tt_boxplot, aes(x = pathogen, y = Rao, col = "Rao-motivated"),  cex = 6)  + #.5, rao  
  geom_point(data = tt_boxplot, aes(x = pathogen, y = Optimal, col = "Single-pathogen motivated"),  cex = 6)  + 
  geom_point(data = tt_boxplot, aes(x = pathogen, y = Measles, col = "Measles-motivated") , cex = 6)  + 
  geom_violin(data = hw_boxplot, aes(x = pathogen, y = lf_vector), draw_quantiles = c(0.25, 0.5, 0.75), fill = "gray95") +
  # geom_boxplot(data = lf_boxplot, aes(x = pathogen, y = lf_vector), fill = "gray95", width=.6) + 
  geom_point(data = hw_boxplot, aes(x = pathogen, y = Rao), col = "#FF4F00",  cex = 6)  + #.5, rao  
  geom_point(data = hw_boxplot, aes(x = pathogen, y = Optimal), col = "black", cex = 6)  + 
  geom_point(data = hw_boxplot, aes(x = pathogen, y = Measles), col = "#0094C6",  cex = 6)  + 
  scale_color_manual(values = c("Rao-motivated" = "#FF4F00", "Single-pathogen motivated" = "black",
                             "Measles-motivated" = "#0094C6" )) +
  scale_fill_manual(values = c("Rao-motivated" = "#FF4F00", "Single-pathogen motivated" = "black",
                        "Measles-motivated" = "#0094C6" )) +
  guides(color =guide_legend(title="Strategy")) +
  xlab("Pathogen") +
  ggtitle("b.") +
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
  ylim(c(0, 80))
target_goal 


## Bootstrapping figure d 
# Measles 
head(measles_path1)
print(unique(measles_path1$pathogen))

pt = .5
store_statistic = list()

for(j in 1:200){
  
  sample_replace = measles_path1  %>%
    filter(pathogen == "measles") %>% # arbtriary to sample from
    group_modify(~ .x[sample(nrow(.x), replace = TRUE), ]) %>%
    left_join(measles_path1 , sample_replace, by = "block_r" ,  relationship =
                "many-to-many") %>% # join so we obtain random sample of block across all pathogens
    mutate(pathogen = pathogen.y, presence = presence.y,  even_prevalence = even_prevalence.y,
           pop_presence = pop_presence.y, fraction = fraction.y) %>%
    dplyr::select(pathogen, block_r, presence, pathogen, even_prevalence, pop_presence, fraction) %>%
    arrange(pathogen) %>%
    mutate(block_replace = rep(seq(1, 90, 1),  4))
  
  head(sample_replace)
  
  # Calculate Rao diversity at each block. 
  product_results = list()
  
  for (i in 1:block_n) {
    
    print(i)
    # Subset the data for the current block
    block_subset <- sample_replace %>%
      mutate(fraction_mirrored = fraction) %>%
      ungroup() %>%
      filter(block_replace== i) %>%
      dplyr::select(fraction_mirrored, fraction)
    
    # Get unique combinations of pathogen prevalences 
    fractions <- block_subset$fraction
    unique_combinations <- t(combn(fractions, 2, simplify = TRUE))  # Get pairs (no duplicates)
    
    # Compute products for each pair
    products <- apply(unique_combinations, 1, prod)
    
    products_df <- data.frame(
      block_replace = i, 
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
    group_by(block_replace) %>%
    summarize(across(product, ~sum(.x)))
  head(block_prevalence_final_replace)
  
  # Calculate expected rao diversity for a given block.
  block_subset_replace = measles_path1  %>% 
    mutate(even_prevalence_mirrored = even_prevalence) %>%
    ungroup() %>% filter(block_r == 1) %>%
    dplyr::select(even_prevalence_mirrored, even_prevalence)
  
  even_prevalences = block_subset_replace$even_prevalence
  unique_combinations <- t(combn(even_prevalences, 2, simplify = TRUE))  # Get pairs (no duplicates)
  
  # Compute products for each pair
  products <- apply(unique_combinations, 1, prod)
  null_rao = sum(products)
  
  # Compare observed rao to expected rao and add gps data. 
  recode_block = sample_replace %>%
    dplyr::select(block_r, block_replace) %>%
    distinct()
  
  method1_replace = block_prevalence_final_replace %>%
    mutate(rao = product/null_rao) %>%
    distinct() %>%
    left_join(recode_block, by = "block_replace")
  head(method1_replace)
  
  compare_strategy_replace = left_join(sample_replace, method1_replace, by = "block_replace") 
  
  
  path = print(unique(compare_strategy_replace$pathogen))
  
  # motivated by itself
  single_motivated = compare_strategy_replace %>% 
    filter(pathogen == "tt")  %>%
    arrange(-fraction) %>% # Arrange from highest prevalence block to lowest 
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2)) # this should do it
  
  print(head(single_motivated))
  print(max(single_motivated$cum))
 # auc_malf <- trapz(malf_motivated$block_label/block_n, malf_motivated$target)
  
  # rao motivated 
  rao_motivated =  compare_strategy_replace %>% 
    filter(pathogen == "tt") %>%
    arrange(-rao) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
 # auc_malf_rao <- trapz(rao_malf_motivated$block_label/psuid_n, rao_malf_motivated$target)
  
  # for AUC 
  # statistic = auc_malf - auc_malf_rao  
  #statistic = auc_malf 
  #  statistic = auc_malf_rao   
  #  print(statistic)
  
  #blocks needed
 # blocks_needed <- which(single_motivated$above == 2)[1]  # optimal
  blocks_needed <- which(rao_motivated$above == 2)[1] # rao
  
  print(blocks_needed)
  
  statistic = blocks_needed
  store_statistic[[j]]= statistic
  
} 


block_vector <- unlist(store_statistic)
block_list <- data.frame(block_vector) %>%
  drop_na()

quantiles = c(.025, .50, .975)
quant_025 = quantile(block_list$block_vector, probs = quantiles[1])
quant_50 = quantile(block_list$block_vector, probs = quantiles[2])
quant_975 = quantile(block_list$block_vector, probs = quantiles[3])

#################################################################

# Save outputs manually 
# Measles 
# optimal: 38 (36, 39)
# rao:     41 (39, 43)

# Ascaris 
# optimal: 31 (28, 33)
# rao:     35 (32, 39)

# Hookworm
# Optimal: 19 (16, 22)
# Rao:     23 (19, 27)

# T. trichurus
# Optimal: 15 (12, 17)
# Rao:     19 (14, 25)


# make into line plot 
mean = c( 35, 31, 23, 19, 41, 38, 19, 15 )
lower = c(32, 28, 19, 16, 39, 36, 14, 12 )
upper = c(39, 33, 27, 22, 43, 39, 25, 17 )
path = c( "Ascaris lumbricoides", "Ascaris lumbricoides", "Hookworm", "Hookworm", 
          "Measles", "Measles", "Trichuris trichiura", "Trichuris trichiura")
strat = c( "Rao-motivated", "Single-pathogen motivated", "Rao-motivated", "Single-pathogen motivated", "Rao-motivated", "Single-pathogen motivated", "Rao-motivated", "Single-pathogen motivated")


paths = data.frame(mean, lower, upper, path, strat)
paths$path = factor(paths$path, levels = c( "Ascaris lumbricoides", "Hookworm", "Measles", "Trichuris trichiura"))
paths$strat = factor(paths$strat, levels = c("Single-pathogen motivated", "Rao-motivated"))
bootstrap = ggplot(data = paths, aes(x = path, y = mean, col = strat)) + 
  theme_light()+
  geom_point( position=position_dodge(.5), cex = 2) +
  # ylab("Area under efficiency curve")+
  ylab("Clusters to target 50% of disease")+
  theme(axis.text=element_text(size=12)) +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.4,
                position=position_dodge(.5), lwd= 1) + 
  xlab("Pathogen") +
  ggtitle("c.") +
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
  scale_color_manual(values = c( "black", "#FF4F00")) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80))  +
  ylim(c(0, 80))
bootstrap

figure_2_bangl = plot_grid(target_goal, bootstrap, rel_widths = c(.75, .75))
figure_2_bangl


plot_grid(fig2_bangl_eff, figure_2_bangl, labels = c("a"), 
          nrow = 2, rel_heights =c(.5, .3))
# 1200 x 1300

###########################################################
# Map it out now
head(measles_path1)
corr_plot = measles_path1 %>%
  dplyr::select(pathogen, fraction, block_r) %>%
  pivot_wider(names_from = pathogen, values_from = fraction) %>%
  dplyr::select(-block_r)
head(corr_plot)

# map out density function for fun
ggplot(data = corr_plot) +
  geom_histogram(aes(x = al)) + xlim(c(0,.1 )) + theme_bw()

colnames(corr_plot) = c("Ascaris", "Hookworm", "Trichuris", "Measles")
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


# Map methods visually.
# Prepare mapping aesthetics 
pal <- wes_palette("Zissou1", 100, type = "continuous")

countries <- ne_countries(scale = "medium", returnclass = "sf")
# Filter for Cambodia
bangl_map <- countries[countries$name == "Bangladesh", ]


# ######################################
# Rao 
########################################
plot_list <- list()
methods = list(method1 = method1)
head(method1)

for(i in seq_along(methods)){
  
  method_df = methods[[i]]
  
  #Retriece map
  admin_b_study <- get_map("district")
  
  dbgps_block = method_df %>% 
    mutate(lat = med_qgpslat, lon = med_qgpslong) %>%
    dplyr::select(block_r, lon, lat) %>%
    dplyr::distinct()  %>% mutate(block = block_r) %>%
    dplyr::select(-block_r)
  
  dbgps = method_df  %>%
    mutate(lat = med_qgpslat, lon = med_qgpslong) %>%
    ungroup() %>%
    dplyr::select(block_r, rao, lat, lon) %>%
    distinct() %>% mutate(block = block_r)  %>%
    dplyr::select(-block_r)
  
  # Create bonding box from the shapefile. 
  b_boxgrid <- st_make_grid(admin_b_study,
                            n=c(300,300), #this is the number of rows and col to split the shapefile into. 
                            what = "centers" #this means that the function will return the center points of each grid, alternatively the shape file would be returned. 
  )
  
  # Transform cluster centroids into shape file. 
  dbgps_sf_utm = dbgps_block %>%
    ungroup() %>%
    st_as_sf(coords = c("lon", "lat"), crs=4326) %>% # crs=4326 specifies that lat and long are coordinates
    dplyr::select(geometry) %>%  
    st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km") # transforms the spatial data from WGS84 coordinates (longitude and latitude) to UTM coordinates
  
  # Identify a buffer of 10 km (based on Arnold et al. 2024).
  cl_buff <- st_buffer(dbgps_sf_utm, dist= 10) %>%  # creates a buffer of 10km around each center 
    summarise(geometry = st_union(geometry)) %>%  # identify their union, then convert to a single polygon that merges overlapping areas of the buffers into one continuous polygon.
    st_cast("POLYGON") %>% # set the CRS back to WGS84 lat and long coordinates 
    st_transform(crs=4326) # this transforms the geometry back to the WGS84 coordinate reference system (EPSG:4326), which uses degrees of latitude and longitude. 
  
  dbgps = dbgps %>%
    mutate( log_rao = log(rao ))
  
  print(head(dbgps))
  print(summary(dbgps$log_rao))
  
  
  # Fit the model dynamically based on the current map
  # formula_string <- paste(maps_i[i], "~ lat + lon + Matern(1|lat + lon)")
  formula_string <- paste("rao ~ lat + lon + Matern(1|lat + lon)")
  b_fit_rao <- spaMM::fitme(
    as.formula(formula_string),
    data = dbgps,
    family = gaussian(link = "identity")
  )
  
  # Fit the model 
 # b_fit_rao <- spaMM::fitme(rao ~  lat + lon + Matern(1|lat + lon), # Matern spatial random effect which is fixed 
  #                          data = dbgps,
   #                         family = gaussian(link = "identity")) 
  # if you do not specify Matern parameters, spaMM defaults to v = .5 (exponential decay of smoothness)
  # ϕ is the scale parameter, which controls the distance over which the correlation decays
  # Bessel function ensures covariance declines with distance between points, gamma function ensures normalization
 # print(summary(b_fit_rao))
  
  #outputs nu and rho parameters, nu indiciates spatial smoothness, rho indicates spatial influence for the correlations modeled by the Matern function.
  plot_materncor <- plot_matern_corr(b_fit_rao) 
  print(plot_materncor)
  
  # this function predicts the response variable = pred and the variance = pred_var over the grid delineated earlier for each shapefile center. 
  b_preds_rao <- get_grid_preds(input_grid=b_boxgrid, 
                                spamm_model_fit = b_fit_rao)
  b_preds_rao_coords <- st_coordinates(b_preds_rao)
  head(b_preds_rao_coords ) # converts shape file centers back to coordinates.
  
  b_preds_rao_raster <- points_to_raster( #for each of the 100 x 100 gridlines, 
    x = b_preds_rao_coords[,1],
    y = b_preds_rao_coords[,2],
    z = b_preds_rao$pred, # was multiplied by 100 before 
    mask1 = admin_b_study,
    mask2 = cl_buff,
    crop1 = admin_b_study
  )
  
  b_preds_rao_raster2 <- raster_to_tibble(b_preds_rao_raster) %>%
    drop_na() %>%
    mutate(value_exp = exp(value)) %>%
    mutate(value = replace(value, value < 0, 0))
  print(head(b_preds_rao_raster2 ))
  
  max_val = max(b_preds_rao_raster2$value)
  
  map_rao <- ggplot(data= b_preds_rao_raster2) + 
   # geom_sf(data = bangl_map, fill = NA, color = "black", size = 15, lwd = 1) +  # Black border for country
    geom_tile(aes(x=x,y=y,fill=value), na.rm=TRUE) +
   # ggtitle("Method", i) +
    coord_sf(crs = 4326)  +
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
                         frame.colour = "white") )  + # +
    labs(x = "Latitude", y = "Longitude") +
    theme_minimal() +
    ggtitle("Rao's quadratic entropy") +
    theme(
     legend.position = "bottom",
      plot.tag = element_text(face = "bold", size = 16)
    )  +
    theme( plot.title = element_text(colour = "black", size = 14.5, face = "bold")) 
  
  
  
  print(map_rao)
  
  plot_list[[i]] = map_rao
} 

rao_plot <- plot_grid(plotlist = plot_list, ncol = 1)
rao_plot



### make left panel of plot
multi_bangl =plot_grid(rao_plot, cor, ncol = 1, rel_heights = c(.5, .2), labels  = c("a", "b"))
multi_bangl


###################################################
# Plot each pathogen's raw smoothed map 
###################################################
# Number of non-NA measurements assessed total by pathogen 
head(disease_dat)
print(unique(disease_dat$tr))

parasite = disease_dat %>%
  pivot_longer(cols=c('al', 'tt', 'hw', 'giar', 'poscr', 'poseh', 'sth'),
               names_to='pathogen',
               values_to='presence') %>%
  group_by(pathogen, block) %>% # assess by overall prevalence in population
  summarise(non_na_count = sum(presence, na.rm = TRUE), 
            denom = n()) %>%
  filter(pathogen != "sth" & pathogen != "giar" & pathogen != "poscr" & pathogen != "poseh")  %>%
  mutate(block_r = block) %>% dplyr::select(-block) %>% #specify that the blocks here are public
  mutate(cases = non_na_count, noncases = denom - non_na_count)%>%
  dplyr::select(block_r, cases, noncases, pathogen)

# Measles unvaccinated. 
head(vax)
measles = vax %>%
  group_by(block_r) %>%
  summarise(vaccinated = sum(measles, na.rm = TRUE), 
            denom = n()) %>%
  mutate(pathogen = "measles") %>%
  mutate(cases = denom - vaccinated, noncases = vaccinated) %>%
  dplyr::select(block_r, cases, noncases, pathogen)
head(measles)

# Join together 
bangladesh_locations = rbind(parasite, measles) %>%
  left_join(gps_dat, by = "block_r") %>%
  ungroup() %>%
  mutate(pathogen = as.character(pathogen))
head(bangladesh_locations)

maps_i = c( "al", "hw", "tt", "measles")
path_names = c( "Ascaris lumbricoides", "Hookworm",  "Trichuris trichiura", "Measles")
legend_names = c("Proportion infected (%)", "Proportion infected (%)", "Proportion infected (%)", "Fraction unvaccinated (%)")

plot_list <- list()

for(i in 1:length(maps_i)) {
  
  print(paste("Processing", maps_i[i]))
  
  # Prepare inputs for mapping
  dbgps_block = bangladesh_locations %>% 
    mutate(lat = med_qgpslat, lon = med_qgpslong) %>% 
    dplyr::select(block_r, lon, lat) %>%
    dplyr::distinct() 
  
  dbgps = bangladesh_locations %>%
    mutate(lat = med_qgpslat, lon = med_qgpslong) %>%
    ungroup() %>%
   # dplyr::select(block_r, !!maps_i[i], lat, lon) %>%
    dplyr::select(block_r, lat, lon, cases, noncases, pathogen) %>%
    distinct() %>%
    drop_na()  %>% 
    filter(as.character(pathogen) == maps_i[i])
  
  # Create bonding box from the shapefile. 
  b_boxgrid <- st_make_grid(admin_b_study,
                            n=c(300,300), #this is the number of rows and col to split the shapefile into. 
                            what = "centers" #this means that the function will return the center points of each grid, alternatively the shape file would be returned. 
  )
  
  # Transform cluster centroids into shape file. 
  dbgps_sf_utm = dbgps_block %>%
    ungroup() %>%
    st_as_sf(coords = c("lon", "lat"), crs=4326) %>% # crs=4326 specifies that lat and long are coordinates
    dplyr::select(geometry) %>%  
    st_transform("+proj=utm +zone=46 +datum=WGS84 +units=km") # transforms the spatial data from WGS84 coordinates (longitude and latitude) to UTM coordinates
  
  # Identify a buffer of 10 km (based on Arnold et al. 2024).
  cl_buff <- st_buffer(dbgps_sf_utm, dist=10) %>%  # creates a buffer of 10km around each center 
    summarise(geometry = st_union(geometry)) %>%  # identify their union, then convert to a single polygon that merges overlapping areas of the buffers into one continuous polygon.
    st_cast("POLYGON") %>% # set the CRS back to WGS84 lat and long coordinates 
    st_transform(crs=4326) # this transforms the geometry back to the WGS84 coordinate reference system (EPSG:4326), which uses degrees of latitude and longitude. 
  
  print(head(cl_buff))
  
  # Fit the model dynamically based on the current map
  # formula_string <- paste(cbind(cases, noncases), "~ lat + lon + Matern(1|lat + lon)")
  b_fit_rao <- spaMM::fitme(
    #  as.formula(formula_string),
    formula = cbind(cases, noncases) ~  lat + lon + Matern(1|lat + lon), 
    data = dbgps,
    family = binomial(link = "logit"))
  
  #outputs nu and rho parameters, nu indiciates spatial smoothness, rho indicates spatial influence for the correlations modeled by the Matern function.
  plot_materncor <- plot_matern_corr(b_fit_rao) 
  print(plot_materncor)
  
  # this function predicts the response variable = pred and the variance = pred_var over the grid delineated earlier for each shapefile center. 
  b_preds_rao <- get_grid_preds(input_grid=b_boxgrid, 
                                spamm_model_fit = b_fit_rao)
  b_preds_rao_coords <- st_coordinates(b_preds_rao)
  head(b_preds_rao_coords ) # converts shape file centers back to coordinates.
  
  b_preds_rao_raster <- points_to_raster( #for each of the 100 x 100 gridlines, 
    x = b_preds_rao_coords[,1],
    y = b_preds_rao_coords[,2],
    z = b_preds_rao$pred, # was multiplied by 100 before 
    mask1 = admin_b_study,
    mask2 = cl_buff,
    crop1 = admin_b_study
  )
  
  b_preds_rao_raster2 <- raster_to_tibble(b_preds_rao_raster) %>%
    drop_na() 
  print(head(b_preds_rao_raster2 ))
  
  max_val = max(b_preds_rao_raster2$value)
  
  # Plotting
  map <- ggplot(data = b_preds_rao_raster2 %>%
                 # mutate(value = replace(value, value < 0, 0)) %>%
                  drop_na(value) %>%
                  mutate(value = value*100)) +  
    geom_tile(aes(x = x, y = y, fill = value), na.rm = TRUE) +
  #  geom_sf(data = bangl_map, fill = NA, color = "black", size = 15, lwd = 1) +  # Black border for country
    ggtitle(path_names[i]) +  # Now using path_names[i] directly
    coord_sf(crs = 4326) +
    scale_fill_gradientn(colours = pal,
                         limits = c(0, 100),
                         na.value = NA,
                         breaks = seq(0, 100, by = 10), 
                         guide = guide_colorbar(
                           title = legend_names[i],
                         #  title = "Proportion infected (%)",
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


inf_plot <- plot_grid(plotlist = plot_list, ncol = 2, labels = "c")
inf_plot

figure_1_bangl = plot_grid(multi_bangl, inf_plot)
figure_1_bangl

# 1400 x 1000








