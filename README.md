A framework for guiding integrated disease control measures through multipathogen surveillance 

Samantha J. Bents, sjbents@stanford.edu

Overview
Global health programs have traditionally focused on single diseases. There is potential for synergy through integrated intervention delivery, particularly in areas with overlapping geographic disease burden, but there is limited methodology developed for guiding integrated delivery and assessing efficiency gains through integration. Here, we applied a measure of diversity, Rao’s quadratic index, to quantify multipathogen burden across two large-scale surveys: Bangladesh (90 clusters, 2,396 children) and Cambodia (100 clusters, 2,150 women). In both settings, we observed geographic clustering of multiple pathogens, indicating potential for more efficient, integrated disease control strategies. We assessed the efficiency of a multipathogen-targeted strategy compared to traditional single-pathogen approaches by calculating the percent reduction in the number of spatial clusters needed to reach 75% of the disease burden (infections or unvaccinated individuals) in a hypothetical intervention. In Bangladesh, integrating deworming with measles vaccination guided by Rao’s quadratic index improved efficiency by 15% for Ascaris lumbricoides, 31% for hookworm, and 38% for Trichuris trichiura, compared to a measles-focused approach. In Cambodia, a Rao-guided strategy performed similarly to the best single-pathogen strategy for Strongyloides stercoralis, and reduced the number of spatial clusters that would need to be targeted by 57% (lymphatic filariasis), 83% (Plasmodium falciparum), and 59% (Plasmodium vivax). We also found that higher multipathogen burden was significantly associated with lower household wealth, suggesting that Rao-guided strategies may more effectively reach under-resourced populations. These findings support the use of multipathogen burden metrics to guide integrated program delivery, offering potential for greater efficiency in disease control. 

Software requirements

All analyses were conducted in R Studio Version 2024.12.1+563 (2024.12.1+563).

Installation Guide

R Studio can be installed at this link: https://posit.co/download/rstudio-desktop/. Typical installation is under ten minutes using an updated operating system.

Demo

The analysis is separated into four R scripts. 

bangl_manuscript_analysis: R code required to calculate and map multipathogen burden using Bangladesh multipathogen survelliance data. Expected run time: 25 minutes.

cambodia_manuscript_analysis: R code required to calculate and map multipathogen burden using Cambodia multiplex serological data. Expected run time: 25 minutes.

compare_metrics: R code to compare efficiency of strategies guided by Shannon, Gini-Simpson, and Rao diversity. Expected run time: 5 minutes.

wealth_index: R code to compare multipathogen burden to indepdent measures of household wealth in Cambodia and Bangladesh. Expected run time: 5 minutes.



