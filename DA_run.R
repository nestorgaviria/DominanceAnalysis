# code to run the dominance analysis 
# this calls the necessary functions and then applies them 

# VERSION 29.01.2026


# ==============================================================================
#*******************************************************************************
#************  1. SETUP & LIBRARIES      *******************************
#*******************************************************************************
# ==============================================================================
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(readxl, dplyr, tidyr, glmnet, caret, DALEX, ggplot2, ggrepel, tibble)

library(relaimpo)
library(ggplot2)
library(patchwork)
library(readxl)
library(tibble)
library(dplyr)
library(tidyr)


# LOAD FUNCTIONS 
source("/Users/nestormac/Documents/Coding/Research/DA_OLS/DA_function.R")
source("/Users/nestormac/Documents/Coding/Research/DA_OLS/MdlsGenerator.R")

#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
#*******************************************************************************
#************  2.  Define the basic configuration  *********************
#*******************************************************************************
#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
config <- list(
  response_var       = "DR_tkm2y",
  response_transform = "log",       # Options: "log", "none", "sqrt"
  n_top_models    = 50,     # <---
  
  pool = list(
    # 1. Topography: Terrain steepnes & flatness, potential energy for moving material
    Topography = c("Ksn", 
                   "LocalRelief",
                   "Grad_prc90", "Grad_prc10" 
                   ),
     
    # 2. Seismic: Endogenic forces (Shaking)
     Seismic = c("PGA50pc50yr", "Seis_Energy_Mean"),

    # 3. Land cover: The surface interface (vegetation, Ice, or Barrenland)
    Land_cover  = c("Barren_Vegetated",
                    "glacier_fr"),
    # 
    # 4. Surface materials: Material resistance
    Surface_materials   = c("Erodibility_lithological", "Prop_Granitic",
                             "Prop_UnconsolidatedDeposits", "SoilThickness_m"),
     
    # 5. Climate: Atmospheric forcing, transporting agents
    Climate     = c("MAP", "AnnualMax_daily_pr",
                    "MAT_K", "AridityIndex")
    ),
  
  # Define nature of predictors that will be used in the models for their transformation
  # predictors that are strictly Positive 
  Positive_to_Log =  c(
    "Slope", "Ksn", "LocalRelief", "ElevationRange", 
    "PGA50pc50yr", "Seis_Freq_Mean", "Seis_Energy_Mean",
    "Grad_prc90", "Grad_prc10",
    "SlabDepth", "SlabDip", "AreaSlope20AreaSlope2",
    "NPP_kgC_m2","Barren_Vegetated","Grass_ForestShrubs",
    "Erodibility_lithological", "SoilThickness_m",
    "MAP", "AnnualMax_daily_pr", "Freq_FreezeThawDays",
    "MAT_K", "Precip_Seasonality", "AridityIndex"
  ),
  
  # Predictors that are proportions (0 to 1) -> Logit(x)
  Proportions_to_Logit =  c(
    "fr_Slope_Low_2deg", "fr_Slope_High_20deg",
    "LandSl_Hazard_Very_High",
    "glacier_fr",  "SnowCover_MeanAnnualfr", 
    "Prop_Granitic", "Prop_UnconsolidatedDeposits"
  )
)

#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
#*******************************************************************************
#************  3.  load the data that will be analysed  *********************
#*******************************************************************************
#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
df_meteoric             <- read_excel("/Users/nestormac/Documents/Research/Berillyium/DataforPaper/XXX_Data_XXX/Data_for_Models/df_DR/df_DRMeteoric_Model.xlsx")
df_insitu               <- read_excel("/Users/nestormac/Documents/Research/Berillyium/DataforPaper/XXX_Data_XXX/Data_for_Models/df_DR/df_DRInSitu_Model.xlsx")
df_Gauging              <- read_excel("/Users/nestormac/Documents/Research/Berillyium/DataforPaper/XXX_Data_XXX/Data_for_Models/df_DR/df_DR_gauging_Model.xlsx")


#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
#*******************************************************************************
#************  4.  Generate models  *********************
#*******************************************************************************
#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
Models_ensmbl <- generate_models_ensemble(config$pool, n_random_models = 500, 
                                                        seed = 123)

#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
#*******************************************************************************
#************  5.  RUN.      *********************
#*******************************************************************************
#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
results_meteoric <- Mdls_dominance_Analysis_OLS_AIC(
  df = df_meteoric,
  model_list = Models_ensmbl,
  dataset_name = "Meteoric 10Be",
  cfg = config
)

results_InSitu<- Mdls_dominance_Analysis_OLS_AIC(
  df = df_insitu,
  model_list = Models_ensmbl,
  dataset_name = "In-Situ 10Be",
  cfg = config
)

results_gauging <- Mdls_dominance_Analysis_OLS_AIC(
  df = df_Gauging,
  model_list = Models_ensmbl,
  dataset_name = "Gauging",
  cfg = config
)

#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
#*******************************************************************************
#************  6.  group results  *********************
#*******************************************************************************
#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈

print(results_meteoric$dominance_summary)
print(results_InSitu$dominance_summary)
print(results_gauging$dominance_summary)


# Combine & Plot
all_dominance <- bind_rows(results_meteoric$dominance_summary, results_InSitu$dominance_summary, results_gauging$dominance_summary)
print(all_dominance)


#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
#*******************************************************************************
#************  7. PLOT          *********************
#*******************************************************************************
#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈

# Define Colors (Consistent Palette)
DRMethods_colors <- c("Meteoric 10Be" = "#009590", 
                      "In-Situ 10Be"   = "#DEC0A3", 
                      "Gauging"       = "#f5f5f5")

category_mapping <- c(
  "Climate" = "Climate",
  "Land_cover" = "Land cover",
  "Topography" = "Topography",
  "Seismic" = "Seismic",
  "Surface_materials" = "Surface materials"
)

all_dominance$Class <- category_mapping[all_dominance$Class]
all_dominance$Class <- factor(
  all_dominance$Class, 
  levels = c( 'Seismic','Surface materials','Land cover', 'Climate',  'Topography')
)

all_dominance$Dataset <- factor(
  all_dominance$Dataset, 
  levels = c( 'Meteoric 10Be', 'In-Situ 10Be', 'Gauging')
)


DA_plot<-ggplot(all_dominance, aes(x = Class, 
                                   y = Mean_Importance, fill = Dataset)) +
  #geom_bar(stat = "identity", width = 0.7, alpha = 0.9) +
  # A. The Stem (From 1.0 to Value)
  geom_segment(aes(x = Class, xend = Class, 
                   y = 0, yend = Mean_Importance), 
               linewidth = 1, color = "grey70") +
  geom_errorbar(aes(ymin =CI_lower, 
                    ymax = CI_upper), 
                width = 0.3,
                linewidth = 0.33, color = "grey30") +
  # B. The Head
  
  geom_point(size = 2.5, shape = 21, color = "grey30", alpha=0.8, stroke = 0.5) +
  scale_fill_manual(values = DRMethods_colors) +
  facet_wrap(~Dataset) +
  coord_flip() +
  #scale_fill_brewer(palette = "Set2") +
  #scale_y_continuous(labels = scales::percent) +
  scale_y_continuous(limits = c(0, 59), 
                     expand = c(0, 0)) +
  labs(title = "Drivers of Denudation: Dominance Analysis",
       #subtitle = "LMG Relative Importance",
       y = "General Dominance \n(% of Explained Variance)", 
       x = "") +
  theme_bw() +
  theme( 
    aspect.ratio = 0.9/1,
    #Set the font type for the pot
    text = element_text(family="sans"),
    #Delete grid from background
    strip.text = element_text(face="bold"),
    panel.grid.major.x = element_line(linewidth = 0.6, color = 'grey90'), 
    panel.grid.minor.x = element_line(linewidth = 0.1, color = 'grey87'),
    panel.grid.major.y = element_blank(),
    #panel.grid.major.y =  element_line(linewidth = 0.2, color = 'grey80', linetype = "dotted"), 
    panel.grid.minor.y = element_blank(),
    #panel.grid.major.y = element_blank(),
    #panel.grid.minor.y = element_blank(),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    #Title
    plot.title = element_text(size=12, colour = "grey30"),
    plot.subtitle = element_text(size=10, colour = "grey40"),
    axis.text.y   = element_text(size=10, colour = "grey40"),
    axis.text.x   = element_text(size=10, colour = "grey40", angle = 0, hjust = 0.5),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    
    axis.title.x  = element_text(size=10, colour = "grey40"),
    axis.title.y  = element_blank(),
    #axis.line=element_blank(),
    plot.margin = unit(c(0.1,0.2,0.1,0.1), "cm"),
    
    legend.position = "none",
  )

DA_plot

output_dir<-"/Users/nestormac/Documents/Research/Berillyium/Plots/2025_Be_paper_plots/RawR"

ggsave(DA_plot, 
       filename = file.path(output_dir,"DominanceAnalysis_CI.pdf"), 
       width = 13, height = 7.7, dpi = 300, units = "cm")

ggsave(DA_plot, 
       filename = file.path(output_dir,"DominanceAnalysis_CI.png"), 
       width = 13, height = 7.7, dpi = 300, units = "cm")


