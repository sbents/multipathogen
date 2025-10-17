########################################################
#~/Library/CloudStorage/Box-Box
#finder --> settings 
#_______________________________________________________
# @Organization --  UCSF 
# @Project -- WASH Benefits Bangladesh multipathogen burden
# @Author -- Samantha Bents, sjbents@stanford.edu, November 30, 2024 
# @Description -- Make wealth index for Bangladesh 
#_______________________________________________________
########################################################
library(here)
source(here("0-config.R"))
source(here("0-functions.R"))

#######################################################

df <- read.csv(file = here("data/bangl/enrollment_wealth/untouched", "washb-bangladesh-enrol-public.csv"))
head(df)

df <- subset(df, select = -c(svyweek, svyyear))

# Vector of variables needed for the function asset_PCA
# removed "asset_phone" because not enough variation
# removed asset_tvbw and asset_tvcol, because they're redundant with asset_tv and cement which is the same with floor

# no WASH-related vars
varlist = c("landacre","roof","walls","floor",
            "elec","asset_radio","asset_refrig","asset_bike","asset_moto",
            "asset_sewmach","asset_tv","asset_wardrobe","asset_table","asset_chair","asset_clock",
            "asset_khat","asset_chouki","asset_mobile")
varlist<-c("dataid","clusterid","hhid","block", varlist)

varlist_fac = c("roof","walls","floor",
                "elec","asset_radio","asset_refrig","asset_bike","asset_moto",
                "asset_sewmach","asset_tv","asset_wardrobe","asset_table","asset_chair","asset_clock",
                "asset_khat","asset_chouki","asset_mobile")

stat_vars = c("landacre","walls","floor",
              "elec","asset_refrig","asset_bike","asset_moto",
              "asset_sewmach","asset_tv","asset_wardrobe","asset_table","asset_chair","asset_khat",
              "asset_chouki","asset_mobile","wealthscore","wealth_tertile")


##############################
# Section 2 ##################
##############################

# Running the function assetPCA that computes the wealth score and divides them
# into tertiles. This function outputs a dataframe with the dataid variables +
# wealth scores and tertiles

assetPCA<-function(df, varlist, varlist_fac, stat_vars, reorder=F ){
  
  #Subset to only needed variables for subgroup analysis
  ret <- df %>%
    subset(select=c(varlist)) ## 22 vars including 4 ids
  
  #Select assets
  ret<-as.data.frame(ret) 
  id<-subset(df, select=c("dataid","clusterid","hhid","block")) ###only ID
  ret_assets_comp<-ret[,which(!(colnames(ret) %in% c("dataid","clusterid","hhid","block")))] ###only asset vars
  ret_assets_comp[,c(2:18)] <- lapply(ret_assets_comp[,c(2:18)], factor)
  
  print(head(ret_assets_comp))
  #Replace character blank with NA
  for(i in 2:ncol(ret_assets_comp)){
    ret_assets_comp[,i]<-ifelse(ret_assets_comp[,i]=="",NA,ret_assets_comp[,i])
  }
  
  #drop rows with no asset data
  id<-id[rowSums(is.na(ret_assets_comp[,5:ncol(ret_assets_comp)])) != ncol(ret_assets_comp)-4,]  
  ret_assets_comp<-ret_assets_comp[rowSums(is.na(ret_assets_comp[,5:ncol(ret_assets_comp)])) != ncol(ret_assets_comp)-4,]  
  
  
  #Drop assets with great missingness
  for(i in 1:ncol(ret_assets_comp)){
    cat(colnames(ret_assets_comp)[i],"\n")
    print(table(is.na(ret_assets_comp[,i])))
    print(class((ret_assets_comp[,i])))
  }
  
  #### asset_clock has 2859 (vs 2692) NAs
  
  ### removing asset_clock
  #cols.dont.want <- c("asset_clock","latseal")
  cols.dont.want <- c("asset_clock")
  ret_assets_comp <- ret_assets_comp[, ! names(ret_assets_comp) %in% cols.dont.want, drop = F]
  
  #create level for missing factor levels
  table(is.na(ret_assets_comp))
  for(i in 2:ncol(ret_assets_comp)){
    ret_assets_comp[,i]<-as.character(ret_assets_comp[,i])
    ret_assets_comp[is.na(ret_assets_comp[,i]),i]<-"miss"
    ret_assets_comp[,i]<-as.factor(ret_assets_comp[,i])
    
  }
  
  ret_assets_comp$landacre[is.na(ret_assets_comp$landacre)] <- mean(ret_assets_comp$landacre, na.rm = T) ###repLace NAs with the mean
  table(is.na(ret_assets_comp))
  
  #Convert factors into indicators
  ret_assets_comp[,2:length(ret_assets_comp)]<-droplevels(ret_assets_comp[,2:length(ret_assets_comp)])
  #ret_assets_compm<-design_matrix(ret_assets_comp[,2:23])
  
  #ret_assets_compm <- data.frame(ret_assets_comp[, ! colnames(ret_assets_comp) %in% c("tubewell","storewat","latown","latslab","roof","walls","floor",
  #     "elec","asset_radio","asset_refrig","asset_bike","asset_moto",
  #    "asset_sewmach","asset_tv","asset_wardrobe","asset_table","asset_chair","asset_khat",
  #  "asset_chouki","asset_mobile")],
  #model.matrix(~tubewell+storewat+latown+latslab+roof+walls+floor+elec+asset_radio+asset_refrig+asset_bike+asset_moto+
  #asset_sewmach+asset_tv+asset_wardrobe+asset_table+asset_chair+asset_khat+asset_chouki+asset_mobile, ret_assets_comp))
  
  Formula <- as.character(paste("~", "ret_assets_comp[,",2,"]", sep=''))
  for (i in 3:length(varlist_fac)) {
    Formula<-as.character(paste(Formula,"+","ret_assets_comp[,",i,"]",sep=""))
  }
  Formula<-as.formula(Formula)
  
  ###removing WASH-related vars
  ret_assets_compm <- data.frame(ret_assets_comp[, ! colnames(ret_assets_comp) %in% varlist_fac],
                                 model.matrix(Formula, ret_assets_comp))
  
  #Remove columns with almost no variance
  if(length(nearZeroVar(ret_assets_compm))>0){
    ret_assets_compm <-ret_assets_compm[,-nearZeroVar(ret_assets_compm)]
  } ## asset_radio and "roof" got removed
  
  ## Convert the data into matrix ##
  ret_assets_compm <-as.matrix(ret_assets_compm)
  
  for(y in 1:ncol(ret_assets_compm)) {
    for(x in 1:nrow(ret_assets_compm)) {
      if (is.na(ret_assets_compm[x,y])) ret_assets_compm[x,y] = 0
    }}
  
  ##Computing the principal component using eigenvalue decomposition ##
  princ.return <- princomp(ret_assets_compm) 
  print(summary(princ.return))
  screeplot(princ.return)
  
#  print(loadings(princ.return))
  
  ## To get the first principal component in a variable ##
  load <- loadings(princ.return)[,1] 
  
  pr.cp <- ret_assets_compm %*% load  ## Matrix multiplication of the input data with the loading for the 1st PC gives us the 1st PC in matrix form. 
  HHwealth <- as.numeric(pr.cp) ## Gives us the 1st PC in numeric form in pr.
  
  #HHwealth=princ.return$scores[,1] ##factor scores from PC 1
  
  #ret_assets_compm<-as.data.frame(ret_assets_compm)
  #ret_assets_compm$HHwealth_nowash<-HHwealth
  
  #Create 3-level household weath index
  tertiles<-quantile(HHwealth, probs=seq(0, 1, 1/3))
  print(tertiles)
  ret_assets_compm<-as.data.frame(ret_assets_compm)
  ret_assets_compm$HHwealth<-HHwealth
  ret_assets_compm$HHwealth_ter<-rep(1, nrow(ret_assets_compm))
  ret_assets_compm$HHwealth_ter[HHwealth>=tertiles[1]]<-1
  ret_assets_compm$HHwealth_ter[HHwealth>=tertiles[2]]<-2
  ret_assets_compm$HHwealth_ter[HHwealth>=tertiles[3]]<-3
  table(ret_assets_compm$HHwealth_ter)
  
  if(reorder==T){
    ret_assets_compm$HHwealth_ter<-factor(ret_assets_compm$HHwealth_ter, levels=c("1", "2","3"))
  }else{
    levels(ret_assets_compm$HHwealth_ter)<-c("1", "2","3")
  }
  
  #Table assets by pca quintile to identify wealth/poverty levels
  d<-data.frame(id, ret_assets_compm)
  #stat_vars = c("landacre","tubewell","storewat","latown","latslab","walls","floor",
  # "elec","asset_refrig","asset_bike","asset_moto",
  # "asset_sewmach","asset_tv","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_mobile","wealthscore_wash","wealthquin")
  
  ###removing WASH-related vars
  #stat_vars = c("landacre","walls","floor",
  #"elec","asset_refrig","asset_bike","asset_moto",
  #"asset_sewmach","asset_tv","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_mobile","wealth_tertile","wealthscore")
  
  colnames(d) <- c("dataid","clusterid","hhid","block", stat_vars)
  
  #Save just the wealth data
  pca.wealth<-d %>% subset(select=c(dataid,clusterid,hhid,block,wealthscore,wealth_tertile))
  pca.wealth$dataid<-as.numeric(as.character(pca.wealth$dataid))
  
  d <-df %>% subset(., select=c("dataid","clusterid","hhid","block"))
  d$dataid<-as.numeric(as.character(d$dataid))
  d<-left_join(d, pca.wealth, by=c("dataid","clusterid","hhid","block"))
  
  return(as.data.frame(d))
  
}

d <- assetPCA(df, varlist = varlist, varlist_fac = varlist_fac, stat_vars = stat_vars, reorder = F)



## 
block_wealth = d %>%
  group_by(block) %>%
  mutate(mean_wealthscore = mean(wealthscore)) %>%
  distinct(block, mean_wealthscore)  %>%
  mutate(block_r = 91 - block) 
head(block_wealth)
  
# Pull in method 
#######################################################
#######################################################
# Method 1: Population-level denominator and assess each STH independently. 
#######################################################
# Load disease data. 
sth = read.csv(file = here("data/bangl/parasites/untouched", "bangl_analysis_parasite.csv"))%>%
  dplyr::select(dataid, clusterid, block, personid, tr, al, tt, giar, hw, sth) 
# sth = readr::read_csv

prot = read.csv(file = here("data/bangl/parasites/untouched", "washb-bangladesh-protozoa-public.csv")) %>%
  dplyr::select(dataid, clusterid, block, personid, tr, poseh , poscr) %>% #entamoeba, crypto
  mutate(poscr = replace(poscr, poscr == 9, NA)) %>% # 9 indicates the sample was missing 
  mutate(poseh = replace(poseh, poseh == 9, NA)) 

# Load measles data from final vacciantion records. 
vax = read.csv(file = here("data/bangl/vaccination/final", "washb_bangl_vax_records_midline.csv"))

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
  arrange(block_r) %>%
  ungroup()

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


######################
# join together 
method1_rao = method1 %>%
  dplyr::select(block_r, rao) %>%
  left_join(block_wealth, by = "block_r")  %>%
  mutate(zero = mean_wealthscore - min(mean_wealthscore)) %>%
  mutate(rescale = zero/max(zero))
head(method1_rao)

pc_bangl = round(cor(method1_rao$rao, method1_rao$rescale, method = "spearman"), digits =2 )
cor.test(method1_rao$rao, method1_rao$rescale, method = "spearman")

bangl_wealth = ggplot(data = method1_rao, aes(x = rescale, y = rao )) + 
  geom_point(cex = 3) +
  ylab("Rao's quadratic index") +
  theme_bw() +
  geom_smooth( col = "darkblue") +
  xlab("Mean wealth score") + ggtitle("a    Bangladesh") +
  annotate("text", x = 0.95, y = max(method1_rao$rao, na.rm = TRUE), 
           label = paste("ρ =", pc_bangl), hjust = 1, size = 6.5) +
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

bangl_wealth



# plot effieciency curve 
pt = .75
block_50 <- list()
efficiency_list <- list()
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
  
  order_measles = print(measles_motivated$block_r)
  
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
  
  # wealth
  measles_wealth_motivated = compare_strategy %>% 
    filter(pathogen == "measles") %>%
    left_join(method1_rao, by = "block_r") %>%
    arrange(mean_wealthscore) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  print( head(rao_wealth_motivated))
  
  measles = ggplot(data = measles_motivated, aes(x = block_label/block_n*90, y = target*100)) +
    geom_line(aes(col = "Measles"), cex = 2) +
    geom_line(data = rao_measles_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Rao"), cex = 2) + 
     geom_line(data = measles_wealth_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Wealth"), cex = 1.2) + # can remove 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent vaccinated (%)") +
    # scale_color_manual(values = c("Measles" = "black", "Rao" = "#FF4F00")) +
    scale_color_manual(values = c("Measles" = "black", "Rao" = "#FF4F00", "Wealth" = "darkblue", lty = "dashed")) +
    guides(color =guide_legend(title="Block prioritization strategy")) +
   # geom_abline(slope=100/90, intercept = 0, lwd =1, lty = 2,  col = "gray72") +
    ggtitle("c      Measles")+
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
  
  ascaris_wealth_motivated = compare_strategy %>% 
    filter(pathogen == "al") %>%
    left_join(method1_rao, by = "block_r") %>%
    arrange(mean_wealthscore) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  
  ascaris = ggplot(data = ascaris_motivated, aes(x = block_label/block_n*90, y = target*100)) +
    geom_line(aes(col = "Ascaris"), cex = 2) +
    geom_line(data = rao_ascaris_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Rao"), cex = 2) + 
  #  geom_line(data = ascaris_measles_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Measles"), cex = 2) + 
     geom_line(data = ascaris_wealth_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Wealth"), cex = 1.2, lty = "dashed") + # can remove 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("Ascaris" = "black", "Rao" = "#FF4F00", "Wealth" = "darkblue")) +
    guides(color =guide_legend(title="Block prioritization strategy")) +
  #  geom_abline(slope=100/90, intercept = 0, lwd =1, lty = 2,  col = "gray72") +
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
  
  trichurus_wealth_motivated = compare_strategy %>% 
    filter(pathogen == "tt") %>%
    left_join(method1_rao, by = "block_r") %>%
    arrange(mean_wealthscore) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  trichurus = ggplot(data = trichurus_motivated, aes(x = block_label/block_n*90, y = target*100)) +
    geom_line(aes(col = "T. trichuris"), cex = 2) +
    geom_line(data = rao_trichurus_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Rao"), cex = 2) + 
  #  geom_line(data = trichurus_measles_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Measles"), cex = 2) + 
    geom_line(data = trichurus_wealth_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Wealth"), cex = 1.2, lty = "dashed") + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("T. trichuris" = "black", "Rao" = "#FF4F00", "Wealth" = "darkblue")) +
    guides(color =guide_legend(title="Block prioritization strategy")) +
   # geom_abline(slope=100/90, intercept = 0, lwd =1, lty = 2,  col = "gray72") +
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
  
  # wealth
  hookworm_wealth_motivated = compare_strategy %>% 
    filter(pathogen == "hw") %>%
    left_join(method1_rao, by = "block_r") %>%
    arrange(mean_wealthscore) %>%
    mutate(total_fraction = sum(fraction)) %>%
    mutate(cum = cumsum(fraction)) %>%
    mutate(target = cum/total_fraction) %>%
    mutate(block_label = seq(1, block_n, 1)) %>%
    mutate(above = ifelse(target < pt, 1, 2))
  
  
  hookworm = ggplot(data = hookworm_motivated, aes(x = block_label/block_n*90, y = target*100)) +
    geom_line(aes(col = "Hookworm"), cex = 2) +
    geom_line(data = rao_hookworm_motivated , aes(x = block_label/block_n*90, y = target*100, col = "Rao"), cex = 2) + 
   # geom_line(data = hookworm_measles_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Measles"), cex = 2) + 
   geom_line(data = hookworm_wealth_motivated, aes(x = block_label/block_n*90, y = target*100, col = "Wealth"), cex = 1.2, lty = "dashed") + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    # scale_color_manual(values = c("Hookworm" = "black", "Rao" = "#FF4F00", "Measles" = "#0094C6")) +
    scale_color_manual(values = c("Hookworm" = "black", "Rao" = "#FF4F00", "Wealth" = "darkblue")) +
    guides(color =guide_legend(title="Block prioritization strategy")) +
  #  geom_abline(slope=100/90, intercept = 0, lwd =1, lty = 2, col = "gray72") +
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
method1_efficiency <- plot_grid(efficiency_list = efficiency_list[[1]])
method1_efficiency 

# make legend for strategy 
legend3 = ggplot(data = parasite_prevalence_block) +
  geom_line(aes(x = block_r, y = presence, col = "Wealth-motivated"), cex = 2) + 
  geom_line(aes(x = block_r, y = fraction, col = "Single-pathogen motivated"), cex = 2) + 
  geom_line(aes(x = block_r, y = even_prevalence, col = "Rao-motivated"), cex = 2) + 
  geom_line(aes(x = block_r, y = pop_presence* even_prevalence, col = "Single-pathogen motivated"), cex = 2) + 
  scale_color_manual(values = c( "Rao-motivated" = "#FF4F00", 
                                 "Wealth-motivated" = "darkblue", "Single-pathogen motivated" = "black")) +
  guides(color =guide_legend(title="Strategy"))  +
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, color = "black"))+ 
  geom_blank() 
legend3
legend_only3 <- gtable::gtable_filter(ggplotGrob(legend3), "guide-box")




wealth_targeting_bangl <- cowplot::plot_grid(
  method1_efficiency,   # Main plot
  legend_only3,           # Extracted legend
  ncol = 1,              # Layout with one column
  rel_heights = c(1, 0.1)  # Adjust relative heights to make space for the legend
)

wealth_targeting_bangl


bangl_wealth_plots = plot_grid(bangl_wealth, wealth_targeting_bangl)
bangl_wealth_plots



###########################################

# Cambodia 
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
  ungroup()
head(method1)
summary(method1$rao)  



############## Wealth 
# Load the .dta file
cam_wealth_data <- read_dta(file = here("projects/6-multipathogen-burden/data/cambodia", "KHPR61FL.DTA"))
head(cam_wealth_data)

cam_wealth_dat = cam_wealth_data %>%
  dplyr::select(hv001, hv270, hv271) %>%
  mutate(dhsclust = hv001) %>%
  group_by(dhsclust) %>%
  summarize(across(hv270:hv271, ~ mean(.x)))

gps_cambodia =   readr::read_csv(file = here("projects/6-multipathogen-burden/data/cambodia", "cambodia_ea_dhs.csv")) 
cam_gps_wealth = left_join(gps_cambodia ,  cam_wealth_dat, by = "dhsclust") %>%
  left_join(gps_dat, by = "dhsclust") %>%
  mutate(psuid = psuid + 1) %>%
  left_join(method1, by = "psuid") 
head(cam_gps_wealth)

# only variables we need
rao_wealth_cambodia = cam_gps_wealth %>%
  dplyr::select(dhsclust, psuid, rao, hv271)
head(rao_wealth_cambodia)

#a = ggplot(data = cam_gps_wealth, aes(x = hv270, y = rao)) +
#  geom_point()+
#  ylab("Rao's quadratic entropy") +
#  theme_bw() +
#  geom_smooth( ) +
#  xlab("Mean wealth score (quintile)") +
#  ggtitle("Cambodia")
head(cam_gps_wealth)
pc_cam = round(cor(cam_gps_wealth$rao, cam_gps_wealth$hv271, method = "spearman"), digits =2 )
cor.test(cam_gps_wealth$rao, cam_gps_wealth$hv271, method = "spearman")

cam_wealth = ggplot(data = cam_gps_wealth, aes(x = (hv271 - min(hv271))/max(hv271 - min(hv271)), y = rao)) +
  geom_point(cex = 3)+
  ylab("Rao's quadratic index") +
  theme_bw() +
  geom_smooth( col = "darkblue") +
  xlab("Mean wealth score") +
  ggtitle("b    Cambodia") +
  annotate("text", x = 0.95, y = max(cam_gps_wealth$rao, na.rm = TRUE), 
           label = paste("ρ =", pc_cam), hjust = 1, size = 6.5) +
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
        panel.border = element_rect(colour = "black", fill=NA))  +
   scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
 
cam_wealth

# plot effiency 

# Compare method efficiency using # of blocks 
pt = .75 # targeting threshold for disease, was .5

block_50 <- list()
efficiency_list <- list()
methods = list(method1 = method1)

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
    mutate(above = ifelse(cum < pt, 1, 2))
  # print(rao_malf_motivated$psuid)
  
  # print(head(rao_malf_motivated))
  
  rao_malf_50 <- which(rao_malf_motivated$above == 2)[1]

  
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
  #  geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "P. falciparum"), cex = 2) +
    geom_line(data = rao_malf_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
   geom_line(data =  malf_wealth_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Wealth"), cex = 1.2, lty = "dashed") + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("P. falciparum" = "black", "Rao" = "#FF4F00",
                                  "Wealth" = "darkblue")) + # "wealth" = "darkblue"
    guides(color =guide_legend(title="Strategy")) +
    ggtitle("d      P. falciparum")+
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
 #   geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "Strongyloides"), cex = 2) +
    geom_line(data = rao_strong_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
   geom_line(data = strong_wealth_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Wealth"), cex = 1.2, lty = "dashed") + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("Strongyloides" = "black", "Rao" = "#FF4F00", "Wealth" = "darkblue" )) + #"wealth" = "darkblue"
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
  #  geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "LF"), cex = 2) +
    geom_line(data = rao_lymph_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
   geom_line(data = lf_wealth_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Wealth"), cex = 1.2, lty = "dashed") + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("LF" = "black", "Rao" = "#FF4F00", "Wealth" = "darkblue")) + # 
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
  #  geom_abline(slope=1, intercept = 0, lwd =1, lty = 2, col = "gray72") +
    geom_line(aes(col = "P. vivax"), cex = 2) +
    geom_line(data = rao_malv_motivated , aes(x = block_label/psuid_n*100, y = target*100, col = "Rao"), cex = 2) + 
     geom_line(data = malv_wealth_motivated, aes(x = block_label/psuid_n*100, y = target*100, col = "Wealth"), cex = 1.2, lty = "dashed") + 
    theme_bw() +
    xlab("Clusters targeted") +
    ylab("Percent disease targeted (%)") +
    scale_color_manual(values = c("P. vivax" = "black", "Rao" = "#FF4F00", "Wealth" = "darkblue" )) + # "wealth" = "darkblue"
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

# Plot efficiency of single-pathogen strategies vs rao strategies. 
#method1_efficiency <- plot_grid(efficiency_list = efficiency_list[[1]], ncol = 2)
method1_efficiency <- efficiency_list[[1]]
method1_efficiency 


# make legend for strategy 
legend1 = ggplot(data = parasite_prevalence_block) +
  geom_line(aes(x = psuid, y = presence, col = "Rao-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = presence, col = "Wealth-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = fraction, col = "P. falciparum-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = even_prevalence, col = "P. vivax-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = pop_presence, col = "LF-motivated"), cex = 2) + 
  geom_line(aes(x = psuid, y = pop_presence* even_prevalence, col = "Single-pathogen motivated"), cex = 2) + 
  scale_color_manual(values = c("Wealth-motivated" = "darkblue" , "Rao-motivated" = "#FF4F00",  "Single-pathogen motivated" = "black")) +
  guides(color =guide_legend(title="Strategy"))  +
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10, color = "black"))+ 
  geom_blank() 
legend1
legend_only <- gtable::gtable_filter(ggplotGrob(legend1), "guide-box")

wealth_targeting_cam <- cowplot::plot_grid(
  method1_efficiency,   # Main plot
  legend_only3,           # Extracted legend
  ncol = 1,              # Layout with one column
  rel_heights = c(1, 0.1)  # Adjust relative heights to make space for the legend
)

wealth_targeting_cam


cam_wealth_plots = plot_grid(cam_wealth,wealth_targeting_cam )
cam_wealth_plots


plot_grid(
         bangl_wealth_plots,cam_wealth_plots,  ncol = 1 )
# 1500 vs 1300


