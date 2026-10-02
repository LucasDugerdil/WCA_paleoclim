#### Meta data ####
# Project: Lake Fazilman Holocene climate, vegetation and human impact reconstructions
# Script: main.R
#
# Use:
# Change the Boolean variables (e.g., Import = T) to import, reconstruct past environmental
# variables and drow the figure related to the article (Dugerdil et al., 2026a)
# 
# Author: Lucas Dugerdil
# ORCID: 0000-0003-0266-564X
#
# Affiliations:
# 1) Univ. Lyon, ENS de Lyon, Université Lyon 1, CNRS, UMR 5276 LGL-TPE,
#    F-69364, Lyon, France
# 2) Université de Montpellier, CNRS, IRD, EPHE, UMR 5554 ISEM,
#    Montpellier, France
#
# Description:
# This script loads and cleans raw proxy data (XRF, magnetic susceptibility,
# spectrocolorimetry, pollen, NPPs, lipid biomarkers) from the Lake Fazilman
# sediment core to prepare datasets for subsequent analyses.
#
# Related publication:
# Dugerdil, L. (2026). Mid-Holocene wet optimum and Early-Late Holocene arid
# phases shaped steppe-forest vegetation and human societies of Uzbekistan:
# multi-proxy evidences from Lake Fazilman. CATENA. https://doi.org/10.1016/j.catena.2025.109759
#
# License:
# Creative Commons Attribution 4.0 International (CC BY 4.0)
#
# Citation:
# If you use this code, please cite:
# Lucas Dugerdil. (2025). LucasDugerdil/Paleo_Fazilman: v1.0.1 (v1.0.1). 
# Zenodo. https://doi.org/10.5281/zenodo.17381815
#
# Created: 2025-08-01
# Last modified: 2025-12-27

#### Import data ####
Import = T
if(Import == T){
  #### Scripts ####
  set.seed(123)
  source("Scripts/setup_packages.R")
  source("Scripts/graphical_functions.R")
  source("Scripts/stat_functions.R")
  if(!exists("map_empty_world")){source("Scripts/Load_maps.R")}
  
  #### WCA metadata ####
  Meta.data <- readRDS("Import/test/Site_full_list_clim_selected.Rds")
  Meta.data.no.filter <- readRDS("Import/test/Site_full_no_filter.Rds")
  List.Ensemble <- readRDS("Import/test/List_Ensemble.Rds")
  List.Ensemble.Filtered <- readRDS("Import/test/List_Ensemble_Filtered.Rds")
  List.Ensemble.clean <- readRDS("Import/test/List_Ensemble_clean.Rds")
  
  #### TraCE-21K data ####
  Tra.cauc <- readRDS("Import/TraCE21K/List_TraCE_Caucasus.Rds")
  Tra.Cas <- readRDS("Import/TraCE21K/List_TraCE_Caspian.Rds")
  Tra.Ir <- readRDS("Import/TraCE21K/List_TraCE_Iran.Rds")
  Tra.ETSH <- readRDS("Import/TraCE21K/List_East_TSH_Iran.Rds")
  Tra.WTSH <- readRDS("Import/TraCE21K/List_TraCE_West_TSH.Rds")
  Meta.data.model.cauc <- readRDS("Import/TraCE21K/Metadata_Caucasus_models.Rds")
  Meta.data.model.Cas <- readRDS("Import/TraCE21K/Metadata_Caspian_models.Rds")
  Meta.data.model.Ir <- readRDS("Import/TraCE21K/Metadata_Iran_models.Rds")
  Meta.data.model.WTSH <- readRDS("Import/TraCE21K/Metadata_West_TSH_models.Rds")
  Meta.data.model.ETSH <- readRDS("Import/TraCE21K/Metadata_East_TSH_models.Rds")
  
  #### Path to create ####
  dir.create("Figures/Main_Figures", recursive = T, showWarnings = F)
  dir.create("Figures/SI", recursive = T, showWarnings = F)
  dir.create("Results/WCA", recursive = T, showWarnings = F)
}

#### Parameters and settings ####
Param.clim = c("MAAT", "MAP") # Climate parameters selection
Limits = c(500, 12500)        # Time limits of the study
Windows.length.final = 200    # Lenght of the time windows for binning
Scaling.final = T             # Apply min-max scaling or not
Anomaly.final = F             # Use anomaly or not
Add.smooth.final = F          # Add a loess smooth or not
Smooth.final = 0.4            # Loess smooth parameter            

#### Select figures to plots ####
Figures = T

Tab.1 = F; Fig.3 = F; Fig.4 = F; Fig.5 = F; Fig.6 = F; Fig.7 = F; Fig.8 = F; Fig.9 = F; Fig.10 = F; Fig.11 = F; Fig.12 = F; Fig.14 = T
# Tab.1 = T; Fig.3 = T; Fig.4 = T; Fig.5 = T; Fig.6 = T; Fig.7 = T; Fig.8 = T; Fig.9 = T; Fig.10 = T; Fig.11 = T; Fig.12 = T; Fig.14 = F
# Tab.1 = T; Fig.3 = F; Fig.4 = F; Fig.5 = F; Fig.6 = F; Fig.7 = F; Fig.8 = F; Fig.9 = F; Fig.10 = F; Fig.11 = F; Fig.12 = F; Fig.14 = F

# Fig.S9 = T

#### Figures and Tables ####
if(Figures == T){
  #### Fig. 5 + S4 (Caucasus) ####
  if(Fig.5 == T){
    #### Settings ####
    Limits = c(0, 12000)
    Windows.length = 500
    Add.TraCE = T
    Keep.area = "Caucasus"
    Remove.sites <- c("Amtkel", "Amtkel_1", "Amtkel_3", "Sakhare", "Kumisi", "Sibista")
    
    #### Filter ####
    if(is.null(Keep.area) == F){Meta.data.cauc <- Meta.data.no.filter[Meta.data.no.filter$Area %in% Keep.area,]}
    if(is.null(Remove.sites) == F){Meta.data.cauc <- Meta.data.cauc[!Meta.data.cauc$Site %in% Remove.sites,]}
    LEC.cauc.raw <- List.Ensemble[Meta.data.cauc$ID]
    LEC.cauc.raw <- Filter(function(df) {any(!is.na(df[, -1, drop = F]))}, LEC.cauc.raw)
    
    LEC.cauc <- Stacking.quantif(Imput.list = LEC.cauc.raw, Limits = Limits, Filter.cut.off = 200, Interp.step = 100,
                                 Return.list = T, Interpolation = F, Windows.length = Windows.length, Keep.clim = paste(Param.clim, "FT_Ensemble", sep = "_"), New.param.names = Param.clim)
    
    #### TraCE extraction (synthese WCA) ####
    if(Add.TraCE == T){
      LEC.cauc <- c(LEC.cauc, Tra.cauc)
      Group.by <- "Recons.type"
    }
    else{Group.by <- NULL}
    
    #### Multi proxy plot ####
    LEC.cauc <- Stacking.quantif(Imput.list = LEC.cauc.raw, Limits = Limits, Return.list = T, Windows.length = Windows.length)
    Meta.data.cauc2 <- Meta.data.cauc[Meta.data.cauc$ID %in% names(LEC.cauc),]
    LEC.cauc.MAAT <- Prep.2.multiproxy.plot(LEC.cauc, MD = Meta.data.cauc2, Param.clim = "MAAT", Keep.CI = T)
    LEC.cauc.MAP <- Prep.2.multiproxy.plot(LEC.cauc, MD = Meta.data.cauc2, Param.clim = "MAP", Keep.CI = T)
    
    
    p1 <- Plot.multiprox(Mlake = LEC.cauc.MAAT, Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), 
                         Smooth.param = 0.4, Title = "(A) MAAT-proxies", Scale.y.free = T, Clean.labels = T, Add.parenthesis = T, Add.CI = T,
                         Panel.space.y = -4, Event.clim = c( 4200, 8200))
    p2 <- Plot.multiprox(Mlake = LEC.cauc.MAP, Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), 
                         Smooth.param = 0.4, Title = "(A) MAP-proxies", Scale.y.free = T, Clean.labels = T, Add.parenthesis = T, Add.CI = T,
                         Panel.space.y = -4, Event.clim = c( 4200, 8200))
    p3 <- p1 + p2
    W = 900; H = 1000; Save.plot = "Figures/SI/Figure_S4.pdf"
    ggsave(filename = Save.plot, p3, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Multi proxy plot (selected) ####
    Selection.MAAT <- c(4,3,7,8,9,11)
    Selection.MAP <- c(3,7,13,9,12,8)
    Clim.surf.cauc <- Meta.data.cauc2[!duplicated(Meta.data.cauc2$Site_full),]
    Clim.surf.cauc$Lake <- Clim.surf.cauc$Site_full
    
    pMAAT <- Plot.multiprox(Mlake = LEC.cauc.MAAT[Selection.MAAT], Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                            Smooth.param = 0.4, Title = "(A) Caucasus - temperature", 
                            Scale.y.free = T, Add.CI = T, Facet.tit.color = T,
                            Panel.space.y = -10, Size.panel.tit = 11, Leg.pos = "none", Show.line = F,
                            Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("W", "C", "W"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                            Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    pMAP <- Plot.multiprox(Mlake = LEC.cauc.MAP[Selection.MAP], Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                           Smooth.param = 0.4, Title = "(B) Caucasus - moisture", Scale.y.free = T, Add.CI = T, Panel.space.y = -5,
                           Size.panel.tit = 11, Leg.pos = "none", Show.line = F, Facet.tit.color = T,
                           Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("D", "Wt", "D"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                           Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    W = 900; H = 800; Save.plot = "Figures/Main_Figures/Figure_5.pdf"
    p.Select <- pMAAT | pMAP
    ggsave(filename = Save.plot, p.Select, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Export results ####
    dir.create("Results/Caucasus/Site/", recursive = T)
    saveRDS(Meta.data.cauc, "Results/Caucasus/Site/Cores_Caucasus_synthese_WCA.Rds")
  }
  
  #### Fig. 6 + S5 (Iran) ####
  if(Fig.6 == T){
    #### Settings ####
    Limits = c(0, 12000)
    Windows.length = 500
    Remove.sites <- c("Gol_E_Zard", "Sibaki", "Kandahar", "Maharlou", "Parishan", "Seistan", "Pir_Ghar")
    Keep.area <- c("Iran")
    Add.TraCE = T
    
    #### Filter ####
    if(is.null(Keep.area) == F){Meta.data.Ir <- Meta.data.no.filter[Meta.data.no.filter$Area %in% Keep.area,]}
    if(is.null(Remove.sites) == F){Meta.data.Ir <- Meta.data.Ir[!Meta.data.Ir$Site %in% Remove.sites,]}
    LEC.Ir.raw <- List.Ensemble[Meta.data.Ir$ID]
    LEC.Ir.raw <- Filter(function(df) {any(!is.na(df[, -1, drop = F]))}, LEC.Ir.raw)
    
    LEC.Ir <- Stacking.quantif(Imput.list = LEC.Ir.raw, Limits = Limits, Filter.cut.off = 200, Interp.step = 100, Scaling = F,
                               Return.list = T, Interpolation = F, Windows.length = Windows.length, Keep.clim = paste(Param.clim, "FT_Ensemble", sep = "_"), New.param.names = Param.clim)
    
    #### TraCE extraction (synthese WCA) ####
    if(Add.TraCE == T){
      LEC.Ir <- c(LEC.Ir, Tra.Ir)
      Group.by <- "Recons.type"
    }
    else{Group.by <- NULL}
    
    #### Multi proxy plot ####
    LEC.Ir <- Stacking.quantif(Imput.list = LEC.Ir.raw, Limits = Limits, Return.list = T, Windows.length = Windows.length)
    Meta.data.Ir2 <- Meta.data.Ir[Meta.data.Ir$ID %in% names(LEC.Ir),]
    LEC.Ir.MAAT <- Prep.2.multiproxy.plot(LEC.Ir, MD = Meta.data.Ir2, Param.clim = "MAAT", Keep.CI = T)
    LEC.Ir.MAP <- Prep.2.multiproxy.plot(LEC.Ir, MD = Meta.data.Ir2, Param.clim = "MAP", Keep.CI = T)
    
    p1 <- Plot.multiprox(Mlake = LEC.Ir.MAAT, Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), 
                         Smooth.param = 0.4, Title = "(A) MAAT-proxies", Scale.y.free = T, Clean.labels = T, Add.parenthesis = T, Add.CI = T,
                         Panel.space.y = 10, Event.clim = c( 4200, 8200))
    p2 <- Plot.multiprox(Mlake = LEC.Ir.MAP, Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), 
                         Smooth.param = 0.4, Title = "(A) MAP-proxies", Scale.y.free = T, Clean.labels = T, Add.parenthesis = T, Add.CI = T,
                         Panel.space.y = -1, Event.clim = c( 4200, 8200))
    p3 <- p1 + p2
    W = 900; H = 850; Save.plot = "Figures/SI/Figure_S5.pdf"
    ggsave(filename = Save.plot, p3, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Multi proxy plot (selected) ####
    Selection.MAAT <- c(3,1,4,2,5)
    Selection.MAP <- c(3,7,10,1,6,9,12,11)
    
    Clim.surf.Ir <- Meta.data.Ir2[!duplicated(Meta.data.Ir2$Site_full),]
    Clim.surf.Ir$Lake <- Clim.surf.Ir$Site_full
    
    pMAAT <- Plot.multiprox(Mlake = LEC.Ir.MAAT[Selection.MAAT], Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                            Smooth.param = 0.4, Title = "(A) Iran - temperature", 
                            Scale.y.free = T, Add.CI = T, 
                            Panel.space.y = 0, Size.panel.tit = 9, Leg.pos = "none", Show.line = F, Facet.tit.color = T,
                            Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("W", "C", "W"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                            Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    pMAP <- Plot.multiprox(Mlake = LEC.Ir.MAP[Selection.MAP], Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                           Smooth.param = 0.4, Title = "(B) Iran - moisture", Scale.y.free = T, Add.CI = T, 
                           # Surf.val = Clim.surf.Ir,
                           Panel.space.y = -10,
                           Size.panel.tit = 9, Leg.pos = "none", Show.line = F, Facet.tit.color = T,
                           Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("D", "Wt", "D"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                           Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    p.Select <- pMAAT | pMAP
    W = 900; H = 850; Save.plot = "Figures/Main_Figures/Figure_6.pdf"
    ggsave(filename = Save.plot, p.Select, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Export results ####
    dir.create("Results/Iran/Site/", recursive = T)
    saveRDS(Meta.data.Ir, "Results/Iran/Site/Cores_Iran_synthese_WCA.Rds")
  }
  
  #### Fig. 7 + S6 (Caspienne-Aral) ####
  if(Fig.7 == T){
    #### Settings ####
    Limits = c(0, 12000)
    Windows.length = 500
    Add.TraCE = T
    Re.extract.TraCE = T
    Keep.area = c("Caspian", "Touran")
    Remove.sites <- c("Amirkola", "Anzali", "Astaneh", "Langarud", "Klaus", "Aral_07")
    
    #### Filter ####
    if(is.null(Keep.area) == F){Meta.data.Cas <- Meta.data.no.filter[Meta.data.no.filter$Area %in% Keep.area,]}
    if(is.null(Remove.sites) == F){Meta.data.Cas <- Meta.data.Cas[!Meta.data.Cas$Site %in% Remove.sites,]}
    LEC.Cas.raw <- List.Ensemble[Meta.data.Cas$ID]
    LEC.Cas.raw <- Filter(function(df) {any(!is.na(df[, -1, drop = F]))}, LEC.Cas.raw)
    
    LEC.Cas <- Stacking.quantif(Imput.list = LEC.Cas.raw, Limits = Limits, Filter.cut.off = 200, Interp.step = 100,
                                Return.list = T, Interpolation = F, Windows.length = Windows.length, Keep.clim = paste(Param.clim, "FT_Ensemble", sep = "_"), New.param.names = Param.clim)
    
    Meta.data.Cas$Site_full[Meta.data.Cas$Site_full == "GS18"] <- "Caspian (middle)"
    Meta.data.Cas$Site_full[Meta.data.Cas$Site_full == "South Caspian (CP14, GS_05)"] <- "Caspian (south)"
    Meta.data.Cas$Site_full[Meta.data.Cas$Site_full == "Golestan (Gomishan, Shahkileh)"] <- "Golestan lagoon"
    
    #### TraCE extraction (synthese WCA) ####
    if(Add.TraCE == T){
      LEC.Cas <- c(LEC.Cas, Tra.Cas)
      Group.by <- "Recons.type"
    }
    else{Group.by <- NULL}
    
    #### Multi proxy plot ####
    LEC.Cas <- Stacking.quantif(Imput.list = LEC.Cas.raw, Limits = Limits, Return.list = T, Windows.length = Windows.length)
    Meta.data.Cas2 <- Meta.data.Cas[Meta.data.Cas$ID %in% names(LEC.Cas),]
    LEC.Cas.MAAT <- Prep.2.multiproxy.plot(LEC.Cas, MD = Meta.data.Cas2, Param.clim = "MAAT", Keep.CI = T)
    LEC.Cas.MAP <- Prep.2.multiproxy.plot(LEC.Cas, MD = Meta.data.Cas2, Param.clim = "MAP", Keep.CI = T)
    
    p1 <- Plot.multiprox(Mlake = LEC.Cas.MAAT, Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), 
                         Smooth.param = 0.4, Title = "(A) MAAT-proxies", Add.CI = T,
                         Panel.space.y = 5, Event.clim = c( 4200, 8200))
    p2 <- Plot.multiprox(Mlake = LEC.Cas.MAP, Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), 
                         Smooth.param = 0.4, Title = "(A) MAP-proxies", Bar.size = 1, Add.CI = T,
                         Panel.space.y = -4, Event.clim = c( 4200, 8200))
    
    W = 900; H = 1000; Save.plot = "Figures/SI/Figure_S6.pdf"
    p3 <- p1 + p2
    ggsave(filename = Save.plot, p3, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Multi proxy plot (selected) ####
    Selection.MAAT <- c(1,3,2,4,5)
    Selection.MAP <- c(9,8,5,6,4,11)
    # Selection.MAP <- c(9,8,7,5,6,4)
    Clim.surf.Cas <- Meta.data.Cas2[!duplicated(Meta.data.Cas2$Site_full),]
    Clim.surf.Cas$Lake <- Clim.surf.Cas$Site_full
    
    pMAAT <- Plot.multiprox(Mlake = LEC.Cas.MAAT[Selection.MAAT], Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                            Smooth.param = 0.4, Title = "(A) Caspian/Aral - temperature", 
                            Scale.y.free = T, Add.CI = T, Facet.tit.color = T,
                            Panel.space.y = -5, Size.panel.tit = 10, Leg.pos = "none", Show.line = F,
                            Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("W", "C", "W"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                            Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    pMAP <- Plot.multiprox(Mlake = LEC.Cas.MAP[Selection.MAP], Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                           Smooth.param = 0.4, Title = "(B) Caspian/Aral - moisture", 
                           Scale.y.free = T, Facet.tit.color = T, Add.CI = T, Panel.space.y = -5, Size.panel.tit = 10, Leg.pos = "none", Show.line = F,
                           Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("D", "Wt", "D"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                           Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    W = 900; H = 800; Save.plot = "Figures/Main_Figures/Figure_7.pdf"
    p.Select <- pMAAT | pMAP
    ggsave(filename = Save.plot, p.Select, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Export results ####
    dir.create("Results/Caspian/Site/", recursive = T)
    saveRDS(Meta.data.Cas, "Results/Caspian/Site/Cores_Caspian_synthese_WCA.Rds")
  }
  
  #### Fig. 8 + S7 (West. TSH) ####
  if(Fig.8 == T){
    #### Settings ####
    Limits = c(0, 12000)
    Windows.length = 500
    Add.TraCE = T
    Re.extract.TraCE = F
    Keep.area = c("West. TSH")
    Remove.sites <- c("Bir_Uja", "Ak_Terek", "Nizhnee", "Ortok", "Ashubulak", "Verkhnee", "Issyk_Kul")
    
    #### Filter ####
    if(is.null(Keep.area) == F){Meta.data.WTSH <- Meta.data.no.filter[Meta.data.no.filter$Area %in% Keep.area,]}
    if(is.null(Remove.sites) == F){Meta.data.WTSH <- Meta.data.WTSH[!Meta.data.WTSH$Site %in% Remove.sites,]}
    LEC.WTSH.raw <- List.Ensemble[Meta.data.WTSH$ID]
    LEC.WTSH.raw <- Filter(function(df) {any(!is.na(df[, -1, drop = F]))}, LEC.WTSH.raw)
    
    LEC.WTSH <- Stacking.quantif(Imput.list = LEC.WTSH.raw, Limits = Limits, Return.list = T, Windows.length = Windows.length, Keep.clim = paste(Param.clim, "FT_Ensemble", sep = "_"), New.param.names = Param.clim)
    
    Meta.data.WTSH$Site_full[Meta.data.WTSH$Site_full == "Issyk-Kul (Issyk_Kul_IK0011, Issyk_Kul_ISKUC087)"] <- "Issyk-Kul"
    Meta.data.WTSH$Site_full[Meta.data.WTSH$Site_full == "Son-Kul (Son_Kul, SONK_11_D)"] <- "Son-Kul"
    Meta.data.WTSH$Site_full[Meta.data.WTSH$Site_full == "Tonnel’naya cave"] <- "Tonnelnaya cave"
    
    #### TraCE extraction (synthese WCA) ####
    if(Add.TraCE == T){
      LEC.WTSH <- c(LEC.WTSH, Tra.WTSH)
      Group.by <- "Recons.type"
    }
    else{Group.by <- NULL}
    
    #### Multi proxy plot (all) ####
    LEC.WTSH <- Stacking.quantif(Imput.list = LEC.WTSH.raw, Limits = Limits, Return.list = T, Windows.length = Windows.length)
    Meta.data.WTSH2 <- Meta.data.WTSH[Meta.data.WTSH$ID %in% names(LEC.WTSH),]
    LEC.WTSH.MAAT <- Prep.2.multiproxy.plot(LEC.WTSH, MD = Meta.data.WTSH2, Param.clim = "MAAT", Keep.CI = T)
    LEC.WTSH.MAP <- Prep.2.multiproxy.plot(LEC.WTSH, MD = Meta.data.WTSH2, Param.clim = "MAP", Keep.CI = T)
    
    p1 <- Plot.multiprox(Mlake = LEC.WTSH.MAAT, Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), 
                         Smooth.param = 0.4, Title = "(A) MAAT-proxies", Scale.y.free = T, Add.CI = T, Facet.tit.color = T, Rotate.panel.tit = T, Add.parenthesis = T, Clean.labels = T,
                         Panel.space.y = -4, Event.clim = c( 4200, 8200))
    p2 <- Plot.multiprox(Mlake = LEC.WTSH.MAP, Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), Size.panel.tit = 6, 
                         Smooth.param = 0.4, Title = "(A) MAP-proxies", Bar.size = 1, Scale.y.free = T, Add.CI = T, Facet.tit.color = T, Rotate.panel.tit = T, Add.parenthesis = T, Clean.labels = T,
                         Panel.space.y = -4, Event.clim = c( 4200, 8200))
    p3 <- p1 + p2
    W = 1000; H = 1600; Save.plot = "Figures/SI/Figure_S7.pdf"
    ggsave(filename = Save.plot, p3, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Multi proxy plot (selected) ####
    Selection.MAAT <- c(2,7,6,1,3,4,5)
    Selection.MAP <- c(16,4,7,2,6,5,19,8,10,15)
    Clim.surf.WTSH <- Meta.data.WTSH2[!duplicated(Meta.data.WTSH2$Site_full),]
    Clim.surf.WTSH$Lake <- Clim.surf.WTSH$Site_full
    
    pMAAT <- Plot.multiprox(Mlake = LEC.WTSH.MAAT[Selection.MAAT], Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                            Smooth.param = 0.4, Title = "(A) Western Tian Shan - temperature", 
                            Scale.y.free = T, Facet.tit.color = T, Add.CI = T, 
                            Panel.space.y = -10, Size.panel.tit = 12, Leg.pos = "none", Show.line = F,
                            Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("W", "C", "W"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                            Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    pMAP <- Plot.multiprox(Mlake = LEC.WTSH.MAP[Selection.MAP], Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                           Smooth.param = 0.4, Title = "(B) Western Tian Shan - moisture", 
                           Scale.y.free = T, Add.CI = T, Panel.space.y = -5, Facet.tit.color = T,
                           Size.panel.tit = 10, Leg.pos = "none", Show.line = F,
                           Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("D", "Wt", "D"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                           Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    H = 1100; W = 1000; Save.plot = "Figures/Main_Figures/Figure_8.pdf"
    
    p.Select <- pMAAT | pMAP
    ggsave(filename = Save.plot, p.Select, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Export results ####
    dir.create("Results/West_TSH/Site/", recursive = T)
    saveRDS(Meta.data.WTSH, "Results/West_TSH/Site/Cores_WTSH_synthese_WCA.Rds")
  }
  
  #### Fig. 9 + S8 (East. TSH) ####
  if(Fig.9 == T){
    #### Settings ####
    Limits = c(0, 12000)
    Windows.length = 500
    Add.TraCE = T
    Re.extract.TraCE = F
    Keep.area = c("East. TSH")
    Remove.sites <- c("Xiao_Yourdusi", "Xiaoxigou", "Chaiwopu_CKF", "KS15", "XEBLK", "XY17", "Huashuwozi", "Daxigou")
    
    #### Filter ####
    if(is.null(Keep.area) == F){Meta.data.ETSH <- Meta.data.no.filter[Meta.data.no.filter$Area %in% Keep.area,]}
    if(is.null(Remove.sites) == F){Meta.data.ETSH <- Meta.data.ETSH[!Meta.data.ETSH$Site %in% Remove.sites,]}
    LEC.ETSH.raw <- List.Ensemble[Meta.data.ETSH$ID]
    LEC.ETSH.raw <- Filter(function(df) {any(!is.na(df[, -1, drop = F]))}, LEC.ETSH.raw)
    
    LEC.ETSH <- Stacking.quantif(Imput.list = LEC.ETSH.raw, Limits = Limits, Return.list = T, Windows.length = Windows.length, Keep.clim = paste(Param.clim, "FT_Ensemble", sep = "_"), New.param.names = Param.clim)
    
    Meta.data.ETSH$Site_full[Meta.data.ETSH$Site_full == "Yili River (Yili, Yili_TKP)"] <- "Yili River"
    
    #### TraCE extraction (synthese WCA) ####
    if(Add.TraCE == T){
      LEC.ETSH <- c(LEC.ETSH, Tra.ETSH)
      Group.by <- "Recons.type"
    }
    else{Group.by <- NULL}
    
    #### Multi proxy plot (all) ####
    LEC.ETSH <- Stacking.quantif(Imput.list = LEC.ETSH.raw, Limits = Limits, Return.list = T, Windows.length = Windows.length)
    Meta.data.ETSH2 <- Meta.data.ETSH[Meta.data.ETSH$ID %in% names(LEC.ETSH),]
    LEC.ETSH.MAAT <- Prep.2.multiproxy.plot(LEC.ETSH, MD = Meta.data.ETSH2, Param.clim = "MAAT", Keep.CI = T)
    LEC.ETSH.MAP <- Prep.2.multiproxy.plot(LEC.ETSH, MD = Meta.data.ETSH2, Param.clim = "MAP", Keep.CI = T)
    
    p1 <- Plot.multiprox(Mlake = LEC.ETSH.MAAT, Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), 
                         Smooth.param = 0.4, Title = "(A) MAAT-proxies", Scale.y.free = T, Add.CI = T, Facet.tit.color = T, Rotate.panel.tit = T, Add.parenthesis = T, Clean.labels = T, 
                         Panel.space.y = -4, Event.clim = c( 4200, 8200))
    p2 <- Plot.multiprox(Mlake = LEC.ETSH.MAP, Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Leg.pos = "bottom", Color.by.proxy = T, Age.depth = c("Age"), Size.panel.tit = 10, 
                         Smooth.param = 0.4, Title = "(A) MAP-proxies", Bar.size = 1, Scale.y.free = T, Add.CI = T, Facet.tit.color = T, Rotate.panel.tit = T, Add.parenthesis = T, Clean.labels = T, 
                         Panel.space.y = -4, Event.clim = c( 4200, 8200))
    p3 <- p1 + p2
    W = 800; H = 1200; Save.plot = "Figures/SI/Figure_S8.pdf"
    ggsave(filename = Save.plot, p3, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Multi proxy plot (selected) ####
    Selection.MAAT <- c(4,5,3,2,1)
    Selection.MAP <- c(6,4,8,5,9,7,1,2)
    Clim.surf.ETSH <- Meta.data.ETSH2[!duplicated(Meta.data.ETSH2$Site_full),]
    Clim.surf.ETSH$Lake <- Clim.surf.ETSH$Site_full
    
    pMAAT <- Plot.multiprox(Mlake = LEC.ETSH.MAAT[Selection.MAAT], Param.clim = c("MAAT"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                            Smooth.param = 0.4, Title = "(A) Eastern Tian Shan - temperature", 
                            Scale.y.free = T, Add.CI = T, 
                            Panel.space.y = -5, Size.panel.tit = 12, Leg.pos = "none", Show.line = F, Facet.tit.color = T,
                            Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("W", "C", "C"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                            Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    pMAP <- Plot.multiprox(Mlake = LEC.ETSH.MAP[Selection.MAP], Param.clim = c("MAP"), Manual.sort = T, Limites = Limits, Color.by.proxy = T, Age.depth = c("Age"), 
                           Smooth.param = 0.4, Title = "(B) Eastern Tian Shan - moisture", 
                           Scale.y.free = T, Add.CI = T, Panel.space.y = -10,
                           Size.panel.tit = 12, Leg.pos = "none", Show.line = F, Facet.tit.color = T,
                           Event.clim = c( 4200, 8200), Name.zone = c("Late Hol.", "Mid-Hol.", "Early Hol."), Temp.zone = c("Wt", "D", "D"), Zone.clim  = c(0, 4200, 4200, 8200, 8200, 12000), Zone.clim.box = F,
                           Add.parenthesis = T, Clean.labels = T, Alpha.zone.clim = 0, Rotate.panel.tit = T)
    
    H = 1100; W = 1000; Save.plot = "Figures/Main_Figures/Figure_9.pdf"
    p.Select <- pMAAT | pMAP
    ggsave(filename = Save.plot, p.Select, width = W*0.026458333, height = H*0.026458333, units = "cm")
    
    #### Export ####
    dir.create("Results/East_TSH/Site/", recursive = T)
    saveRDS(Meta.data.ETSH, "Results/East_TSH/Site/Cores_ETSH_synthese_WCA.Rds")
    
  }
  
  #### Import manual section of sites ####
  Meta.data.Ir <- readRDS("Results/Iran/Site/Cores_Iran_synthese_WCA.Rds")
  Meta.data.cauc <- readRDS("Results/Caucasus/Site/Cores_Caucasus_synthese_WCA.Rds")
  Meta.data.ETSH <- readRDS("Results/East_TSH/Site/Cores_ETSH_synthese_WCA.Rds")
  Meta.data.WTSH <- readRDS("Results/West_TSH/Site/Cores_WTSH_synthese_WCA.Rds")
  Meta.data.Cas <- readRDS("Results/Caspian/Site/Cores_Caspian_synthese_WCA.Rds")
  Sub.MD.WCA <- rbind(Meta.data.Ir, Meta.data.cauc, Meta.data.ETSH, Meta.data.WTSH, Meta.data.Cas)
  
  Sub.MD.WCA$Area[Sub.MD.WCA$Area == "Touran"] <- "Caspian"
  Sub.MD.WCA$Area[Sub.MD.WCA$Area == "Caspian"] <- "Caspian-Aral"
  
  Meta.data$Area[Meta.data$Area == "Touran"] <- "Caspian"
  Meta.data$Area[Meta.data$Area == "Caspian"] <- "Caspian-Aral"
  Meta.data.WCA <- Meta.data[Meta.data$Area %in% c("Iran", "Caucasus", "Caspian-Aral", "West. TSH", "East. TSH"),]
  Meta.data.WCA.no.filter <- Meta.data.no.filter[Meta.data.no.filter$Area %in% c("Iran", "Caucasus", "Caspian-Aral", "West. TSH", "East. TSH"),]
  Meta.data.WCA.no.filter <- Meta.data.WCA.no.filter[Meta.data.WCA.no.filter$Acces == "T",]
  
  List.Ensemble.clean.WCA.selected <- List.Ensemble.clean[Meta.data.WCA$ID]
  List.Ensemble.clean.WCA.selected <- Filter(function(df) {any(!is.na(df[, -1, drop = F]))}, List.Ensemble.clean.WCA.selected)
  
  Meta.data.WCA.selected <- Meta.data.WCA[Meta.data.WCA$ID %in% names(List.Ensemble.clean.WCA.selected),]
  Meta.data.WCA.selected <- Meta.data.WCA.selected[Meta.data.WCA.selected$ID %in% Sub.MD.WCA$ID,]
  
  List.Ensemble.clean.WCA.selected$`477|Wenqu|Pol`$MAAT <- NA
  List.Ensemble.clean.WCA.selected <- List.Ensemble.clean[Meta.data.WCA$ID]
  Sub.MD.WCA.in.stack <- Meta.data.WCA[Meta.data.WCA$ID %in% names(List.Ensemble.clean.WCA.selected),]
  
  saveRDS(List.Ensemble.clean.WCA.selected, "Results/WCA/List_Ensemble_clean_WCA_selected.Rds")
  saveRDS(Sub.MD.WCA.in.stack, "Results/WCA/Meta_data_WCA_study_stacked_list.Rds")
  saveRDS(Meta.data.WCA.selected, "Results/WCA/Meta_data_WCA_study_selected.Rds")
  
  List.Ensemble.clean.WCA.selected <- c(List.Ensemble.clean.WCA.selected, Tra.Ir, Tra.ETSH, Tra.Cas, Tra.cauc, Tra.WTSH)
  Meta.data.WCA.selected$Recons.type <-"1. Data"
  Meta.data.WCA.selected <- rbind(Meta.data.WCA.selected, Meta.data.model.Ir, Meta.data.model.ETSH, Meta.data.model.Cas, Meta.data.model.cauc, Meta.data.model.WTSH)
  Meta.data.WCA.selected$Area[Meta.data.WCA.selected$Area == "Touran"] <- "Caspian-Aral"
  Meta.data.WCA.selected$Area[Meta.data.WCA.selected$Area == "Caspian"] <- "Caspian-Aral"
  
  #### Fig. 3 (Hol. maps) ####
  if(Fig.3 == T){
    Data.map.WCA <- List.Ensemble.clean.WCA.selected[!grepl("\\|TRA", names(List.Ensemble.clean.WCA.selected))]
    MD.map.WCA <- Meta.data.WCA.selected
    
    Data.map.WCA <- Data.map.WCA[grepl("\\|Pol", names(Data.map.WCA))]
    Data.map.WCA <- Data.map.WCA[!grepl("SonKu", names(Data.map.WCA))]
    
    A <- Map.panel.climat(
      List.df = Data.map.WCA,
      Timespan = c(200, 4200, 8200, 11700),
      Meta.data = MD.map.WCA,
      Display.plot = T, Altitudinal.gradient = F,
      Period.lab = c("Late Hol.", "Mid-Hol.", "Early Hol."),
      Color.by.rank = T, Add.emptysite = T, Use.Meta.data.4.empty.site = T,
      Map.dot.size = 4.5, Empty.site.shape = 4, DEM.alpha = .6,
      Emphase.proxy = "Pollen", Remove.overlap = T, Show.label = F, Overlap.seuil = 80000,
      Title.size = 12, Dot.alpha = 1, 
      Add.DEM = F, # the DEM has been removed in this repository because of the raster size
      Map.leg.pos = c(0.002, 0.999), Map.limits = c(41,89,26,49),
      Save.plot = "Figures/Main_Figures/Figure_3.pdf", H = 750, W = 850
    )
  }
    
  #### Fig. 4 (alt. gradient) ####
  if(Fig.4 == T){
    Data.map.WCA <- List.Ensemble.clean.WCA.selected[!grepl("\\|TRA", names(List.Ensemble.clean.WCA.selected))]
    MD.map.WCA <- Meta.data.WCA.selected
    
    Data.map.WCA <- Data.map.WCA[grepl("\\|Pol", names(Data.map.WCA))]
    Data.map.WCA <- Data.map.WCA[!grepl("SonKu", names(Data.map.WCA))]
    
    B <- Map.panel.climat(
      List.df = Data.map.WCA,
      Timespan = c(200, 4200, 8200, 11700),
      Meta.data = MD.map.WCA, 
      Display.plot = F, Show.stats = T,
      Altitudinal.gradient = T, Elevation.range = c(-300, 3800),
      Period.lab = c("Late Hol.", "Mid-Hol.", "Early Hol."),
      Color.by.rank = T, Use.Meta.data.4.empty.site = T,
      Map.dot.size = 4.5, Empty.site.shape = 4, DEM.alpha = .6,
      Emphase.proxy = "Pollen", Remove.overlap = T, Show.label = F, Overlap.seuil = 80000,
      Title.size = 12, Dot.alpha = 1, Add.DEM = F,
      Map.leg.pos = c(0.002, 0.999), Map.limits = c(41,89,26,49),
      Save.plot = "Figures/Main_Figures/Figure_4.pdf", H = 500, W = 850
    )
  }
  
  #### Fig. 10 (stacks) ####
  if(Fig.10 == T){
    Gr.area.WCA <- Plot.stack.groups(List.clim = List.Ensemble.clean.WCA.selected, Meta.data = Meta.data.WCA.selected, Group = "Area", 
                                     Subgroup = "Recons.type", Limits = Limits,
                                     Plot.2.merge = F, Display.map = T, Map.limits = c(41,89,26,49),
                                     Add.count.per.bin = F, Windows.length = Windows.length.final, 
                                     Save.data.Rds = "Results/WCA/Holocene_WCA_clusters_zone.Rds",
                                     Show.strip = T, Facet.by.group = T, Display.CI = T, Leg.ncol = NULL, Width.pour.map = 30, CI.alpha = .4,
                                     Facet.order = c("Iran", "Caucasus", "Caspian-Aral", "Touran", "West. TSH", "East. TSH"), 
                                     Remove.extrems = F, Color.by.clim = T, Error.alpha = .2, Add.error.zone = T,
                                     Add.vlines = c(4200, 8200), Scaling = Scaling.final, Anomaly = Anomaly.final, Smooth.param = Smooth.final, Add.smooth = Add.smooth.final,
                                     Name.param.clim = c("Temperature", "Moisture"), Add.emptysite = T, Empty.site.from.Meta.data = T, Map.dot.size = 3.2,
                                     Save.plot = "Figures/Main_Figures/Figure_10.pdf", H = 1000, W = 1115,
                                     Save.plot.Rds = "Results/WCA/Holocene_ACA_clusters_WCA_plot.Rds")
  }
  
  #### Fig. 11 (mat. correl.) ####
  if(Fig.11 == T){
    #### Matrice correlation (prep) ####
    Order.zn <- c("Iran", "Caucasus", "Caspian-Aral", "West. TSH", "East. TSH")
    Keep.ages <- intersect(unique(Gr.area.WCA.data$Age), unique(Gr.area.WCA.mode$Age))
    Insol <- Insol.Uz[Insol.Uz$Age %in% Keep.ages,] 
    Insol <- subset(Insol, select = - c(Age))
    
    Gr.area.WCA.data.MAP <- Gr.area.WCA.data[c("Age", "MAP", "Group")] 
    Gr.area.WCA.data.MAP <- dcast(Gr.area.WCA.data.MAP, Age ~ Group, value.var = "MAP")
    Gr.area.WCA.data.MAP <- subset(Gr.area.WCA.data.MAP, select = - c(Age))
    Gr.area.WCA.data.MAP <- Gr.area.WCA.data.MAP[Order.zn]
    
    Gr.area.WCA.mode.MAP <- Gr.area.WCA.mode[c("Age", "MAP", "Group")] 
    Gr.area.WCA.mode.MAP <- dcast(Gr.area.WCA.mode.MAP, Age ~ Group, value.var = "MAP")
    Gr.area.WCA.mode.MAP <- subset(Gr.area.WCA.mode.MAP, select = - c(Age))
    Gr.area.WCA.mode.MAP <- Gr.area.WCA.mode.MAP[Order.zn]
    
    Gr.area.WCA.data.MAP <- Gr.area.WCA.data.MAP[intersect(row.names(Gr.area.WCA.data.MAP), row.names(Gr.area.WCA.mode.MAP)),]
    Gr.area.WCA.mode.MAP <- Gr.area.WCA.mode.MAP[intersect(row.names(Gr.area.WCA.data.MAP), row.names(Gr.area.WCA.mode.MAP)),]
    row.names(Insol) <- row.names(Gr.area.WCA.data.MAP)
    
    Gr.area.WCA.data.MAAT <- Gr.area.WCA.data[c("Age", "MAAT", "Group")] 
    Gr.area.WCA.data.MAAT <- dcast(Gr.area.WCA.data.MAAT, Age ~ Group, value.var = "MAAT")
    Gr.area.WCA.data.MAAT <- subset(Gr.area.WCA.data.MAAT, select = - c(Age))
    Gr.area.WCA.data.MAAT <- Gr.area.WCA.data.MAAT[Order.zn]
    
    Gr.area.WCA.mode.MAAT <- Gr.area.WCA.mode[c("Age", "MAAT", "Group")] 
    Gr.area.WCA.mode.MAAT <- dcast(Gr.area.WCA.mode.MAAT, Age ~ Group, value.var = "MAAT")
    Gr.area.WCA.mode.MAAT <- subset(Gr.area.WCA.mode.MAAT, select = - c(Age))
    Gr.area.WCA.mode.MAAT <- Gr.area.WCA.mode.MAAT[Order.zn]
    
    Gr.area.WCA.data.MAAT <- Gr.area.WCA.data.MAAT[intersect(row.names(Gr.area.WCA.data.MAAT), row.names(Gr.area.WCA.mode.MAAT)),]
    Gr.area.WCA.mode.MAAT <- Gr.area.WCA.mode.MAAT[intersect(row.names(Gr.area.WCA.data.MAAT), row.names(Gr.area.WCA.mode.MAAT)),]
    row.names(Insol) <- row.names(Gr.area.WCA.data.MAAT)
    
    Insol <- Insol[c(4,2,3,1)]
    
    #### Matrice correlation (plot) ####
    Perm = T; My_perm = 9#999
    H = 500*.8; W = 800*.8; Save.matcor.WCA = "Figures/Main_Figures/Figure_11.pdf"
    pdf(file = Save.matcor.WCA, width = W*0.01041666666667, height = H*0.01041666666667)
    par(mfrow = c(2,3))
    my_color <- colorRampPalette(c("darkblue", "grey95", "grey95", "darkred"))(20)
    
    M1 <- Mat.corel.CWT.clim(Gr.area.WCA.data.MAAT, Gr.area.WCA.data.MAAT, Display = "lower",
                             Display.pval = "pch", my_color = my_color, Xlab.rot = 0, Xlab.pos = "t",
                             Permutation.test = Perm, Nb.permutations = My_perm, I.confiance = .95,
                             Disp.R = "number", Label = F, Average = F,  Bar.pos = "n", Title = "(A) Data vs. Data")
    
    M2 <- Mat.corel.CWT.clim(Gr.area.WCA.data.MAAT, Gr.area.WCA.mode.MAAT, #Display = "lower",
                             Display.pval = "pch", my_color = my_color, Xlab.rot = 0, Xlab.pos = "t", 
                             Permutation.test = Perm, Nb.permutations = My_perm, I.confiance = .95,
                             Disp.R = "number", Label = F, Average = F,  Bar.pos = "n", Title = "(B) Data vs. TraCE-21K")
    
    M3 <- Mat.corel.CWT.clim(Gr.area.WCA.data.MAAT, Insol, #Display = "lower",
                             Display.pval = "pch", my_color = my_color, Xlab.rot = 0, Xlab.pos = "t", 
                             Permutation.test = Perm, Nb.permutations = My_perm, I.confiance = .95, 
                             Disp.R = "number", Label = F, Average = F,  Bar.pos = "n", Title = "(C) Data vs. Insolation")
    
    M4 <- Mat.corel.CWT.clim(Gr.area.WCA.data.MAP, Gr.area.WCA.data.MAP, Display = "lower",
                             Display.pval = "pch", my_color = my_color, Xlab.rot = 0, Xlab.pos = "t",
                             Permutation.test = Perm, Nb.permutations = My_perm, I.confiance = .95,
                             Disp.R = "number", Label = F, Average = F,  Bar.pos = "n")
    
    M5 <- Mat.corel.CWT.clim(Gr.area.WCA.data.MAP, Gr.area.WCA.mode.MAP, #Display = "lower",
                             Display.pval = "pch", my_color = my_color, Xlab.rot = 0, Xlab.pos = "n",
                             Permutation.test = Perm, Nb.permutations = My_perm, I.confiance = .95,
                             Disp.R = "number", Label = F, Average = F,  Bar.pos = "n")
    
    M6 <- Mat.corel.CWT.clim(Gr.area.WCA.data.MAP, Insol, #Display = "lower",
                             Display.pval = "pch", my_color = my_color, Xlab.rot = 0, Xlab.pos = "n",
                             Permutation.test = Perm, Nb.permutations = My_perm, I.confiance = .95,
                             Disp.R = "number", Label = F, Average = F,  Bar.pos = "n")
    
    dev.off()
    
  }
  
  #### Fig. 12 (stack vs. insol.) ####
  if(Fig.12 == T){
    Compar.stack.inso.WCA <- Compar.stack.inso(Mdata = Gr.area.WCA$`1. Data`, Minso = Insol.Uz,
                                               List.keep = list(c("MAP", "East. TSH", "Winter"), 
                                                                c("MAP", "East. TSH", "Summer"),
                                                                c("MAP", "Iran", "Spring"),
                                                                c("MAP", "West. TSH", "Summer"),
                                                                c("MAP", "West. TSH", "Spring")), 
                                               Line.inso.width = 1.5, CI.alpha = .2,
                                               Facet.order = c("Iran", "West. TSH", "East. TSH"), Leg.pos = "bottom", Leg.ncol = 4, Title.size = 10,
                                               Reverse.val = c("Spring"),
                                               Add.vlines = c(4200, 8200), Limits = c(0, 12000),
                                               Remove.ages = c(600, 12000, 12200),
                                               Save.plot = "Figures/Main_Figures/Figure_12.pdf", H = 700, W = 500)
  }
  
  #### Fig. 14 (map bias) ####
  if(Fig.14 == T){
    #### Settings ####
    Meta.data = Meta.data.WCA.no.filter
    
    Proxy.type = "Proxy"
    Map.dot.size = 13
    Map.leg.pos = "bottomright"
    Map.limits = c(49.5, 94.5, 24, 53)
    Show.country.borders = T
    Plot.Density = T
    Plot.Diversity = T
    Diversity.panel.kernel = T
    Diversity.panel.buffer = F
    Add.DEM = T
    DEM.alpha = .5
    Kernel.radius = 250 # km
    Density.resolution = 200 # Nb x and y cells
    Save.plot = "Figures/Main_Figures/Figure_14.pdf"
    H = 600*.9
    W = 720*.9
    
    #### Load SIG data ####
    Path.SIG = "/media/lucas-dugerdil/DDT4/Documents/Recherche/SIG/"
    if(exists("study_area") == F){
      study_area <- st_read(paste0(Path.SIG, "Projets/ACA/Western_Asia/WA_zone_caspian.shp"))
      study_area <- st_set_crs(study_area, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      Caspian <- st_read(paste0(Path.SIG, "Projets/Circum_caspienne/Caspian_extend_political/Caspian_actual_political.shp"))
      Caspian <- st_set_crs(Caspian, "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
      
    }
    
    names(Meta.data)[grep("^lon|^Lon", names(Meta.data))] <- "Longitude"
    names(Meta.data)[grep("^lat|^Lat", names(Meta.data))] <- "Latitude"
    proxy_sf <- st_as_sf(Meta.data, coords = c("Longitude", "Latitude"), crs = 4326)
    
    aea_crs <- "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=85 +datum=WGS84 +units=m +no_defs"
    study_area_proj <- st_transform(study_area, crs = aea_crs)
    proxy_proj      <- st_transform(proxy_sf, crs = aea_crs)
    
    study_vect <- vect(study_area_proj)
    proxy_vect <- vect(proxy_proj)
    r_template <- rast(study_vect, resolution = Density.resolution*1000)
    proxy_raster <- rasterize(proxy_vect, r_template, field = 1)
    
    #### Distance ####
    if(Plot.Density == F){
      dist_raster <- distance(proxy_raster)
      dist_raster_km <- dist_raster / 1000
      dist_masked <- mask(dist_raster_km, study_vect)
      dist_df <- as.data.frame(dist_masked, xy = TRUE, na.rm = TRUE)
      colnames(dist_df)[3] <- "distance_km"
      My_plot <- geom_raster(data = dist_df, aes(x = x, y = y, fill = distance_km))
      Fill_scale <- scale_fill_gradientn(name = "Distance to proxy (km)", colours = c("grey90", "#FF9E3D", "darkorange", "darkred")) 
    }
    
    #### Density (kernel) ####
    if(Plot.Density == T){
      library(spatstat)
      library(scales)
      ppp_points <- as.ppp(proxy_proj["geometry"])
      Window(ppp_points) <- as.owin(study_area_proj)
      density_spatstat <- density(ppp_points, dimyx = Density.resolution)
      density_stars <- stars::st_as_stars(density_spatstat)
      density_sf <- st_as_sf(density_stars)
      density_sf <- st_set_crs(density_sf, "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=85 +datum=WGS84 +units=m +no_defs")
      density_sf$v <- density_sf$v * 1e12  # now ~0–100 or so
      My_plot <- geom_sf(data = density_sf, aes(fill = v), colour = NA)
      
      Fill_scale <-   scale_fill_gradientn(
        name = expression("Proxy density" ~ (10^12 ~ "sites.km"^-2)),
        # colours = c("white", "grey90", "#FF9E3D", "darkorange", "darkred"), transform = "sqrt")
        colours = c("white", "#ff80e5ff", "#d400aaff", "#800066ff", "#550044ff"), transform = "sqrt")
    }
    
    #### Limits ####
    if(is.null(Map.limits) == F){
      bbox_ll <- st_bbox(c(
        xmin = Map.limits[1],
        xmax = Map.limits[2],
        ymin = Map.limits[3],
        ymax = Map.limits[4]), crs = 4326)
      
      bbox_ll_sf <- st_as_sfc(bbox_ll)
      bbox_proj <- st_transform(bbox_ll_sf, st_crs(study_area_proj))
      bbox_proj_vals <- st_bbox(bbox_proj)
      
      My_coords <- coord_sf(
        xlim = c(bbox_proj_vals["xmin"], bbox_proj_vals["xmax"]),
        ylim = c(bbox_proj_vals["ymin"], bbox_proj_vals["ymax"]), expand = F)
    }
    else{My_coords <- NULL}
    
    #### Diversity (local buffer approach) ####
    if(Diversity.panel.buffer == T){
      proxy_types <- unique(proxy_proj[[Proxy.type]])
      r_list <- list()
      
      names(proxy_proj)[names(proxy_proj) == Proxy.type] <- "Proxy"
      for(ptype in proxy_types){
        p_sf <- filter(proxy_proj, Proxy == ptype)
        p_vect <- vect(p_sf[,0])  # keep only geometry
        r <- rasterize(p_vect, r_template, field=1)
        r[is.na(r)] <- 0
        r_list[[ptype]] <- r
      }
      
      proxy_stack <- rast(r_list)  # multilayer raster
      names(proxy_stack) <- proxy_types
      
      radius_m <- Kernel.radius*1000  # 200 km
      w <- focalMat(proxy_stack[[1]], radius_m, type="circle")
      
      proxy_smooth_list <- lapply(1:nlyr(proxy_stack), function(i){
        focal(proxy_stack[[i]], w, fun=sum, na.rm=TRUE)
      })
      
      proxy_stack_smooth <- rast(proxy_smooth_list)
      names(proxy_stack_smooth) <- proxy_types
      
      shannon_raster <- app(proxy_stack_smooth, fun=function(x){
        p <- x / sum(x)
        p <- p[p>0]
        if(length(p)==0) return(0)
        -sum(p*log(p))
      })
      shannon_raster <- mask(shannon_raster, study_vect)
      dens_df <- as.data.frame(shannon_raster, xy=TRUE, na.rm=TRUE)
      colnames(dens_df)[3] <- "shannon"
      
      # Scale_fill_2 <- scale_fill_gradientn(name = "Proxy diversity (Shannon)", colours = c("grey90", "#26a186", "royalblue", "#4a51b7", "#461d64"))
      Scale_fill_2 <- scale_fill_gradientn(name = "Proxy diversity (Shannon)", colours = c("grey90", "#26a186","#4a51b7", "#461d64"))
      Plot_div <- geom_raster(data = dens_df, aes(x = x, y = y, fill = shannon), alpha = Alpha.layer)
    }
    else{Scale_fill_2 <- NULL}
    
    #### Diversity (kernel probability approach) ####
    if(Diversity.panel.kernel == T){
      library(spatstat)
      proxy_types <- unique(proxy_proj[[Proxy.type]])
      
      ppp_points <- as.ppp(proxy_proj["geometry"])
      Window(ppp_points) <- as.owin(study_area_proj)
      dens_list <- list()
      
      for(ptype in proxy_types){
        p_sub <- proxy_proj[proxy_proj[[Proxy.type]] == ptype, ]
        
        if(nrow(p_sub) > 0){
          ppp_sub <- as.ppp(p_sub["geometry"])
          Window(ppp_sub) <- as.owin(study_area_proj)
          dens <- density(ppp_sub, sigma = Kernel.radius*1000, dimyx = Density.resolution)
          dens_list[[ptype]] <- dens
        }
      }
      
      stars_list <- lapply(dens_list, stars::st_as_stars)
      dens_stack <- do.call(c, stars_list)
      dens_df <- as.data.frame(dens_stack, xy = TRUE)
      dens_mat <- as.matrix(dens_df[ , -(1:2)])
      row_sum <- rowSums(dens_mat)
      
      prob_mat <- dens_mat / row_sum
      prob_mat[is.na(prob_mat)] <- 0
      
      shannon <- -rowSums(prob_mat * log(prob_mat + 1e-12))
      shannon[is.nan(shannon)] <- 0
      shannon[shannon == 0] <- NA
      
      # Threshold
      total_density <- rowSums(dens_mat)
      density_norm <- total_density / max(total_density, na.rm = T)
      # threshold <- quantile(total_density, 0.1, na.rm = T)  # or fixed value
      # # threshold <- 0.5
      # shannon[total_density < threshold] <- NA
      shannon_weighted <- shannon * density_norm
      shannon_weighted <- Scaling.fun(shannon_weighted, 0, 1)
      dens_df$shannon2 <- shannon_weighted
      # effective_div <- exp(shannon)
      # dens_df$div_weighted <- effective_div * density_norm
      # total_density <- rowSums(dens_mat)
      # density_norm <- total_density / max(total_density, na.rm = T)
      # shannon_weighted <- shannon * density_norm
      # dens_df$shannon <- shannon_weighted
      # dens_df <- dens_df[!is.na(dens_df$shannon), ]
      
      dens_df$shannon <- shannon
      # stop()
      # dens_df <- dens_df[c(1,2,ncol(dens_df))]
      
      dens_df$alpha_val <- as.character(ifelse(dens_df$shannon2 > 0.4, 5,
                                               ifelse(dens_df$shannon2 > 0.3, 4, 
                                                      ifelse(dens_df$shannon2 > 0.1, 3, 
                                                             ifelse(dens_df$shannon2 > 0.05, 2, 
                                                                    ifelse(dens_df$shannon2 > 0.03, 1, 0))))))
      dens_df$alpha_val[is.na(dens_df$alpha_val)] <- 0
      Scale_fill_2 <- scale_fill_gradientn(name = "Proxy diversity (Shannon)", colours = c("grey90", "#26a186","#4a51b7", "#461d64"), na.value = NA)
      Alpha.layer = 0.8 
      Plot_div <- geom_raster(data = dens_df, aes(x = x, y = y, fill = shannon), alpha = Alpha.layer)
    }
    
    #### Display DEM ####
    if(Add.DEM == T){
      library(raster)
      Path.world.DEM = paste(Path.SIG, "MNT/WorldClim2_MNT/wc2.1_10m_elev.tif", sep = "")
      DEM = raster::raster(Path.world.DEM)
      DEM <- projectRaster(DEM, crs = st_crs(study_area_proj)$proj4string, method="bilinear")
      DEM <- raster::crop(DEM, raster::extent(bbox_proj_vals["xmin"], bbox_proj_vals["xmax"], bbox_proj_vals["ymin"], bbox_proj_vals["ymax"]))
      DEM_lowres <- raster::aggregate(DEM, fact = 1, fun = mean, na.rm = TRUE)
      
      breaks <- c(-Inf, 0, 250, 500, 1000, 2000, 3000, 4000, 5000, Inf)
      DEM_disc <- raster::cut(DEM_lowres, breaks = breaks, include.lowest = T)
      
      DEM_df <- raster::as.data.frame(DEM_disc, xy = TRUE, na.rm = TRUE)
      colnames(DEM_df) <- c("x", "y", "elev_class")
      
      density_sf$alpha_val <- as.character(ifelse(density_sf$v > 25, 5, 
                                                  ifelse(density_sf$v > 20, 4, 
                                                         ifelse(density_sf$v > 15, 3, 
                                                                ifelse(density_sf$v > 10, 2, 
                                                                       ifelse(density_sf$v > 5, 1, 0))))))
      My_plot <- geom_sf(data = density_sf, aes(fill = v, alpha = alpha_val), colour = NA)
      if(Diversity.panel.kernel == T){Plot_div <- geom_raster(data = dens_df, aes(x = x, y = y, fill = shannon, alpha = alpha_val))}
      alpha.scale <- scale_alpha_discrete(c("0" = 0.9, "1" = 0.7, "2" = 0.8, "3" = 0.9, "4" = 0.95, "5" = 1), guide = "none")
      My_DEM <-  geom_raster(data = DEM_df, aes(x = x, y = y, fill = elev_class), alpha = DEM.alpha)
      DEM_scale <- scale_fill_gradientn(colors = c(NA, "grey90", "grey80", "grey60", "tan3", "sienna4", "white"), guide = "none")
      Alpha.layer = 0.8
    }
    else{My_DEM <- NULL; DEM_scale <- NULL; Alpha.layer = 1; alpha.scale <- NULL}
    
    #### Add neighbors ####
    if(Show.country.borders == T){
      library("rnaturalearth")
      library("rnaturalearthdata")
      world <- ne_countries(scale = "medium", returnclass = "sf")
      world_proj <- st_transform(world, st_crs(study_area_proj))
      Country_out <- geom_sf(data = world_proj, fill = "grey", color = "grey40", linewidth = 0.4)
      Country_in <- geom_sf(data = world_proj, fill = NA, color = "grey40", linewidth = 0.4)
      Country_in2 <- geom_sf(data = world_proj, fill = NA, color = "grey40", linewidth = 0.4)
    }
    else{Country_out <- NULL; Country_in <- NULL}
    
    #### Graphical params ####
    if(Map.leg.pos[1] == "bottomright"){
      Map.leg.pos = c(0.999, 0.001);
      My_just = c("right", "bottom")
    }
    else{My_just = c(0, 1)}
    
    My_theme <- theme(plot.background = element_blank(), panel.grid = element_blank(), panel.border = element_rect(colour = "black", linewidth = 1),
                      legend.position = Map.leg.pos, legend.justification = My_just, legend.box = "horizontal",
                      panel.background = element_rect(fill = "#7BACCF"), axis.title = element_blank(),
                      legend.background = element_rect(fill = alpha("white", 1), color = "black"), legend.box.background = element_blank())
    
    My_guide <- guides(
      fill = guide_colorbar(direction = "horizontal", title.position = "top",
                            title.hjust = 0, 
                            barwidth = unit(5.5, "cm"), barheight = unit(0.5, "cm")))
    
    #### Map plot (distance ou density) ####
    p.map <- ggplot() +
      Country_out + 
      alpha.scale +
      My_DEM + DEM_scale +
      new_scale_fill() +
      geom_sf(data = Caspian, shape = 4, fill = "#7BACCF", color = "black", linewidth = 0.6) +
      My_plot +
      Country_in +
      geom_sf(data = Caspian, shape = 4, fill = NA, color = "black", linewidth = 0.6) +
      Fill_scale +
      geom_sf(data = proxy_proj, shape = 4, fill = "black", color = "black", size = 2, stroke = 0.4) +
      geom_sf(data = study_area_proj, fill = NA, color = "black", linewidth = 0.6) +
      My_guide + My_coords +
      theme_minimal() + My_theme
    
    #### Map (diversity) ####
    if(Plot.Diversity == T){
      p.map.2 <- ggplot() +
        Country_out +
        alpha.scale +
        geom_sf(data = study_area_proj, fill= "white", color="black", linewidth=0.8) +
        # Country_in2 +
        My_DEM + DEM_scale +
        geom_sf(data = Caspian, shape = 4, fill = "#7BACCF", color = "black", linewidth = 0.6) +
        new_scale_fill()+
        Plot_div +
        Scale_fill_2 +
        Country_in2 +
        geom_sf(data = Caspian, shape = 4, fill = NA, color = "black", linewidth = 0.6) +
        
        geom_sf(data = proxy_proj, shape = 4, fill = "black", color = "black", size = 2, stroke = 0.4) +
        geom_sf(data = study_area_proj, fill=NA, color="black", linewidth=0.8) +
        My_coords + My_guide +
        theme_minimal() + My_theme
      
      p.map <- p.map + p.map.2
      W = 2*W
    }
    #### Export ####
    if(is.null(Save.plot) == F){
      if(is.null(W) == F & is.null(H) == F){ggsave(p.map, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", limitsize = F)}
      else{ggsave(Save.plot)}}
  }
  #### Fig. S9 (insolation) ####
  if(Fig.S9 == T){
    Insol.Uz <- Seasonal.insol(Import.path = "Import/Insolation/Laskar_2004", Imported.raw.data = T,
                               Save.Rds = "Results/Insolation/Orbit_seasonnal.Rds",
                               Limits = c(500, 12000), Remove.annual = T, CI.lines = F,
                               Show.CI = T, Show.each.month = T, Anomaly = T, Date.Anomaly = -100,
                               Manual.y.line = c(500, 4200, 8200), Show.Facet = F, Keep.only = c("Summer", "Winter", "Spring", "Autumn"),
                               H = 500, W = 700, Save.plot = "Figures/SI/Figure_S9.pdf")
  }
  else{Insol.Uz <- readRDS("Results/Insolation/Orbit_seasonnal.Rds")}
  
  Gr.area.WCA <- readRDS("Results/WCA/Holocene_WCA_clusters_zone.Rds")
  Sub.MD.WCA.in.stack <- readRDS("Results/WCA/Meta_data_WCA_study_stacked_list.Rds")
  Meta.data.WCA.selected <- readRDS("Results/WCA/Meta_data_WCA_study_selected.Rds")
  List.Ensemble.clean.WCA.selected <- readRDS("Results/WCA/List_Ensemble_clean_WCA_selected.Rds")
  
  Gr.area.WCA <- split(Gr.area.WCA, Gr.area.WCA$Subgroup)
  Gr.area.WCA.data <- Gr.area.WCA$`1. Data`
  Gr.area.WCA.mode <- Gr.area.WCA$`2. Model`
  
  #### Export Latex table for sites ####
  if(Tab.1 == T){
    Sub.MD.WCA <- Sub.MD.WCA[Sub.MD.WCA$Proxy_classes != "TraCE-21Ka",]
    print(paste("**** In the WCA study we have", nrow(Sub.MD.WCA), "sites. ****"))
    print(paste("**** In the WCA study we have", nrow(Sub.MD.WCA[Sub.MD.WCA$Proxy_classes == "Pollen",]), "sites. ****"))
    print(paste("**** In the WCA STACK we have", nrow(Sub.MD.WCA.in.stack), "sites. ****"))
    print(paste("**** In the WCA STACK we have", nrow(Sub.MD.WCA.in.stack[Sub.MD.WCA.in.stack$Proxy_classes == "Pollen",]), "sites. ****"))
    
    dir.create("Tables/", showWarnings = F)
    A <- Table.Latex.settings(Table = Sub.MD.WCA, Round.coord = 3, Save.csv = T, Save.Latex = T, Names.2.rows = T, Abbrev.OK = T, Print.tab.head = T,
                              Group.order = c("Caucasus", "Iran", "Caspian-Aral", "West. TSH", "East. TSH"),
                              Keep.columns = c("Area", "ID", "Site_full", "Latitude", "Longitude", "Altitude", "Archive", "Proxy", "Begin", "End", "Res", "Param.clim", "DB", "Reference"),
                              Rename.columns = c("Region", "Site ID", "Site name", "Lat. (°N)", "Long. (°E)", "Elev. (m~a.s.l.)", "Archive type", "Proxy type", "Time begin", "Time end", "Res. (yrs)", "Climate proxy", "Database", "References"),
                              LateX.caption = "Settings of the palaeoclimate archives from the Western part of the Arid Central Asia, selected in this study.",
                              Save.path = "Tables/Table_1",
                              Latex.lab = "Table_Setting")}
  
  
}

########################################### A FAIRE ############################################################
A.faire = F
if(A.faire == T){

#### Load functions ####
# if (!exists("Clim.param.extraction", mode = "function")){source("Scripts/Load_func.R"); Load.my.functions.2(Script = "Scripts/Climat_extract.R", Function.name = "Clim.param.extraction")}
# if (!exists("Extract.TraCE.21ka", mode = "function")){source("Scripts/Load_func.R"); Load.my.functions.2(Script = "Scripts/TraCE21ka_extract.R", Function.name = "Extract.TraCE.21ka")}
# if (!exists("Fossil.MAT.prep", mode = "function")){source("Scripts/Load_func.R"); Load.my.functions.2(Script = "Scripts/Pollen_fun_trans.R", Function.name = "Fossil.MAT.prep")}
# if (!exists("Fossil.corresp.surface", mode = "function")){source("Scripts/Load_func.R"); Load.my.functions.2(Script = "Scripts/Pollen_fun_trans.R", Function.name = "Fossil.corresp.surface")}
# if (!exists("FT.core", mode = "function")){source("Scripts/Load_func.R"); Load.my.functions.2(Script = "Scripts/Pollen_fun_trans.R", Function.name = "FT.core")}

#### Set NULL values ####
Run.process = F; Run.plots = F; Display.errors = F; Recalculate.models = F; Recalculate.Ensembles = F; Fuzzy.Cm = F; Extract.Clim.WC = F

Min.end = NULL; Min.begin = NULL; Min.Resolution = NULL; Limits = NULL
Type_conversion = "Corresp_FT_Odile"
Keep.proxy = NULL; Keep.archive = NULL; Keep.DB = NULL; Remove.proxy = NULL; Remove.archive = NULL; Test.Subset = NULL
Remove.site = NULL; Keep.area = NULL; Remove.area = NULL; Keep.site = NULL; Only.sites = NULL; Keep.country = NULL; Remove.country = NULL

#### Import Training databases ####
if(exists("MAT.ACADB") == F){
  print("**** Import FT DB surface. ****")
  Version.DB <- "DB5498"
  MAT.ACADB    <- readRDS("Results/ACA/Pollen/Func_trans/Surface/MAT_ACADB.Rds")
  WAPLS.ACADB  <- readRDS("Results/ACA/Pollen/Func_trans/Surface/WAPLS_ACADB.Rds")
  BRT.ACADB    <- readRDS("Results/ACA/Pollen/Func_trans/Surface/BRT_ACADB.Rds")
  MAT.TUSDB    <- readRDS("Results/Uzbekistan/Pollen/Func_trans/Surface/MAT_TUSDB.Rds")
  WAPLS.TUSDB  <- readRDS("Results/Uzbekistan/Pollen/Func_trans/Surface/WAPLS_TUSDB.Rds")
  RF.TUSDB     <- readRDS("Results/Uzbekistan/Pollen/Func_trans/Surface/RF_TUSDB.Rds")
  BRT.TUSDB    <- readRDS("Results/Uzbekistan/Pollen/Func_trans/Surface/BRT_TUSDB.Rds")
  MAT.COSTDB   <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "MAT_COSTDB.Rds", sep = "/"))
  WAPLS.COSTDB <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "WAPLS_COSTDB.Rds", sep = "/"))
  BRT.COSTDB   <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "BRT_COSTDB.Rds", sep = "/"))
  MAT.TAIGDB   <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "MAT_TAIGDB.Rds", sep = "/"))
  WAPLS.TAIGDB <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "WAPLS_TAIGDB.Rds", sep = "/"))
  RF.TAIGDB    <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "RF_TAIGDB.Rds", sep = "/"))
  BRT.TAIGDB   <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "BRT_TAIGDB.Rds", sep = "/"))
  MAT.WASTDB   <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "MAT_WASTDB.Rds", sep = "/"))
  WAPLS.WASTDB <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "WAPLS_WASTDB.Rds", sep = "/"))
  RF.WASTDB    <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "RF_WASTDB.Rds", sep = "/"))
  BRT.WASTDB   <- readRDS(paste("Results/World_DB/Pollen/Func_trans/Surface", Version.DB, "BRT_WASTDB.Rds", sep = "/"))
  MAT.CAUCDB    <- readRDS("Results/Caucasus/Pollen/Func_trans/Surface/MAT_CAUCDB.Rds")
  WAPLS.CAUCDB  <- readRDS("Results/Caucasus/Pollen/Func_trans/Surface/WAPLS_CAUCDB.Rds")
  BRT.CAUCDB    <- readRDS("Results/Caucasus/Pollen/Func_trans/Surface/BRT_CAUCDB.Rds")
  
  # List.FT.calib <- list("BRT.TUSDB"   = BRT.TUSDB, "WAPLS.TUSDB" = WAPLS.TUSDB, "MAT.TUSDB"   = MAT.TUSDB)
  # List.FT.calib <- list("BRT.ACADB"   = BRT.ACADB, "WAPLS.ACADB" = WAPLS.ACADB, "MAT.ACADB"   = MAT.ACADB)
  List.FT.calib <- list("BRT.TAIGDB"   = BRT.TAIGDB, "WAPLS.TAIGDB" = WAPLS.TAIGDB, "MAT.TAIGDB"   = MAT.TAIGDB)
  # List.FT.calib <- list("BRT.COSTDB"   = BRT.COSTDB, "WAPLS.COSTDB" = WAPLS.COSTDB, "MAT.COSTDB"   = MAT.COSTDB)
  # List.FT.calib <- list("BRT.CAUCDB"   = BRT.CAUCDB, "WAPLS.CAUCDB" = WAPLS.CAUCDB, "MAT.CAUCDB"   = MAT.CAUCDB)
  
  List.FT.calib <- list("BRT.TUSDB"   = BRT.TUSDB, "BRT.ACADB"  = BRT.ACADB, "BRT.COSTDB" = BRT.COSTDB, "BRT.TAIGDB" = BRT.TAIGDB, "BRT.CAUCDB" = BRT.CAUCDB,
                        "WAPLS.TUSDB" = WAPLS.TUSDB, "WAPLS.ACADB"= WAPLS.ACADB, "WAPLS.COSTDB" = WAPLS.COSTDB, "WAPLS.TAIGDB" = WAPLS.TAIGDB, "WAPLS.CAUCDB" = WAPLS.CAUCDB,
                        "MAT.TUSDB"   = MAT.TUSDB, "MAT.ACADB"  = MAT.ACADB, "MAT.COSTDB" = MAT.COSTDB, "MAT.TAIGDB" = MAT.TAIGDB, "MAT.CAUCDB" = MAT.CAUCDB)
}

if(exists("Calib.list") == F){
  List.pollen.aqua <- data.frame(read.csv(file="Import/ACA/Site/Paleo_sequences/ACA_remove_pollen.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
  List.loess.PCA <- data.frame(read.csv(file="Import/ACA/Geoch/Loess/Loess_PCA_direction.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
  DB.odile.Po.COST <- data.frame(read.csv(file="Results/World_DB/Pollen/Func_trans/Surface/SIG/Odile_COST_Po.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
  Uz.MP.FT <- readRDS("Results/Uzbekistan/Pollen/Func_trans/Surface/Uz_MP_FT_clean.Rds")
  DB.odile.Po.TAIGDB <- data.frame(read.csv(file="Results/World_DB/Pollen/Func_trans/Surface/SIG/Odile_TAIGDB_Po.csv", sep=",", dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) # DB Odile, world
  ACADB.Po <- readRDS("Results/ACA/Pollen/Func_trans/Surface/ACA_MP_FT_clean.Rds")
  CAUCDB.Po <- readRDS("Results/Caucasus/Pollen/Func_trans/Surface/CAUC_MP_FT_clean.Rds")
  
  Calib.list <- list(COSTDB = DB.odile.Po.COST, TAIGDB = DB.odile.Po.TAIGDB, TUSDB = Uz.MP.FT, ACADB = ACADB.Po, CAUCDB = CAUCDB.Po)
}

#### Table conversion pollen ####
Tab.Cor.Name.Arm <-data.frame(read.csv(file="Import/Armenia/Pollen/Zarishat/Corresp_pollen_FT_Armenia.csv",sep=",",dec=".",header=T,row.names=1, stringsAsFactors = FALSE))
Tab.Cor.Name.Aral  <-   data.frame(read.csv(file="Import/Kazakhstan/Pollen/Cores/Corresp_pollen_FT_Aral.csv",sep=",",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) 
Tab.Cor.Name.Legacy <-  data.frame(read.csv(file="Import/ACA/Pollen/LegacyPollen_ACA_corresp_tab_auto_manuel.csv",sep=",",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) [c(2:5)]
Tab.Cor.Name.Neotoma <- data.frame(read.csv(file="Import/World_DB/Pollen/Neotoma/Indexes/Corresp_pollen_FT_Neotoma.csv",sep=",",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) 
Tab.Cor.Name.ERAPDB <- data.frame(read.csv(file="Import/World_DB/Pollen/ERAPDB/Indexes/Corresp_pollen_FT_ERAPDB.csv",sep=",",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE))
Tab.Cor.Name.Gottingen <- setNames(data.frame(read.csv(file="Import/Mongolia/Pollen/Func_trans/Corresp_pollen_MAT_Mongolia_Gottingen.csv",sep=",",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE))[c(2:4)], names(Tab.Cor.Name.Aral))
Tab.Cor.Name.EMBSeCBIO  <- data.frame(read.csv(file="Import/World_DB/Pollen/EMBSeCBIO/EMBSeCBIO_ACA_corresp_tab_clean.csv",sep=",",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE)) 

names(Tab.Cor.Name.Legacy)[names(Tab.Cor.Name.Legacy) == "Short.names"] <- "Nom"
Tab.Cor.Name.list <- list(Tab.Cor.Name.Aral, Tab.Cor.Name.Neotoma, Tab.Cor.Name.Legacy, Tab.Cor.Name.ERAPDB[c(1:3)], Tab.Cor.Name.Gottingen, Tab.Cor.Name.EMBSeCBIO, Tab.Cor.Name.Arm)
Tab.Cor.Name <- suppressMessages(suppressWarnings(Reduce(function(x, y) full_join(x, y), Tab.Cor.Name.list)))
Tab.Cor.Name <- Tab.Cor.Name[!duplicated(Tab.Cor.Name$Nom),]
# Tab.Cor.Name <- Tab.Cor.Name[c(1,4,2,3,5)]
A <- Tab.Cor.Name[duplicated(Tab.Cor.Name$Nom) | duplicated(Tab.Cor.Name$Nom, fromLast = T),]
A <- A[order(A$Nom),]

A <- setdiff(names(ACADB.Po), unique(Tab.Cor.Name$Corresp_FT_Odile))
# print(A)

Check.error <- Tab.Cor.Name[c(1,2)]
Check.error$Simple_name <- gsub("-|/", " ", Check.error$Nom)
Check.error$Simple_name <- gsub("^cf |^cf\\. |type| sp$| sp\\.|\\.| undiff|\\(|\\)|\\?", "", Check.error$Simple_name)
Check.error$Simple_name <- gsub(" $", "", Check.error$Simple_name)
Check.error$Simple_name <- gsub("^ ", "", Check.error$Simple_name)
Check.error$Simple_name <- gsub("  |   ", " ", Check.error$Simple_name)
Check.error$Simple_name <- gsub(" sp$", "", Check.error$Simple_name)
Check.error$Simple_name <- gsub(" subsp.*", "", Check.error$Simple_name)
Check.error$Simple_name <- gsub(" sect.*", "", Check.error$Simple_name)
Check.error$Simple_name <- gsub(" subg.*", "", Check.error$Simple_name)
Check.error$Family <- gsub(" .*", "", Check.error$Simple_name)
Check.error$Genus <- Check.error$Family
Check.error$Genus[grepl("eae$", Check.error$Genus)] <- ""
Check.error$Family[!grepl("eae$", Check.error$Family)] <- ""

# First error
A <- Check.error[c(2,3)]
B <- A[(duplicated(A$Corresp_FT_Odile) | duplicated(A$Corresp_FT_Odile, fromLast = T)) & !duplicated(A),]
B <- B[order(B$Corresp_FT_Odile),]

lapply(Tab.Cor.Name.list, function(x) x[x$Corresp_FT_Odile == "ZIZYPHUS",])

# Check.missing = "Results/ACA/Pollen_missing"
# Missing.path <- list.files(path = Check.missing, pattern = ".csv", full.names = T, all.files = T, recursive = F)
# Check.missing <- bind_rows(lapply(Missing.path, read.csv))
# Check.missing <- Check.missing[!duplicated(Check.missing$x),]
# write.table(Check.missing, "Results/ACA/Pollen_missing/merged.csv", sep = ",")

Table.SI.pollentypes = F
if(Table.SI.pollentypes == T){
  #### Setting ####
  
  Save.csv = F; Save.Latex = T; Latex.lab = "SI_pollen_table_FT"; Add.lines = F
  LateX.caption = "Matching table between row pollen type and aggregated pollen types reliable for pollen-based climate reconstructions."
  Save.path = "Results/ACA/Export_LaTeX/Papier_synthese_WCA/SI_table_pollen_type_FT"
  
  #### Cleaning ####
  Clean.tab.PT <- Tab.Cor.Name[c("Nom", Type_conversion)]
  Clean.tab.PT <- Clean.tab.PT[!is.na(Clean.tab.PT[[Type_conversion]]),]
  Clean.tab.PT <- Clean.tab.PT[Clean.tab.PT[[Type_conversion]] != "ABS",]
  names(Clean.tab.PT)[2] <- "TYPE"
  
  TLegacy <-  data.frame(read.csv(file="Import/ACA/Pollen/LegacyPollen_ACA_corresp_tab_auto_manuel.csv",sep=",",dec=".",header=T, row.names = 1, stringsAsFactors = FALSE))
  Clean.tab.PT$Legacy <- TLegacy$Nom[match(Clean.tab.PT$Nom, TLegacy$Short.names)]
  Clean.tab.PT$Nom[!is.na(Clean.tab.PT$Legacy)] <- Clean.tab.PT$Legacy[!is.na(Clean.tab.PT$Legacy)]
  Clean.tab.PT <- subset(Clean.tab.PT, select = - c(Legacy))
  Clean.tab.PT$Nom <- gsub("  |   |   |^cf |^cf. |\\d| reworked|\\?|\\(|\\)|µ|/| I$| II$| III$", " ", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- gsub(" $|  |^ ", " ", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- gsub("\\+", "-", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- gsub(" type| t$", "-type", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- gsub(" subg ", " subg. ", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- gsub(" sect ", " sect. ", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- gsub(" var ", " var. ", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- gsub("ë", "e", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- gsub(" (subsp|sub|ssp)\\.? ", " subsp. ", Clean.tab.PT$Nom)
  
  Clean.tab.PT$Nom[Clean.tab.PT$Nom == "Monocotyledones"] <- "Monocotyledons"
  Clean.tab.PT$Nom[Clean.tab.PT$Nom == "Pinus haplostellate"] <- "Pinus Haploxylon"
  Clean.tab.PT$Nom[Clean.tab.PT$Nom == "Pinus haplostele"] <- "Pinus Haploxylon"
  
  Clean.tab.PT <- Clean.tab.PT[!duplicated(Clean.tab.PT$Nom),]
  Clean.tab.PT$Sim_name <- gsub("-type| spp.| sp.| undiff.| group| undiff", "", Clean.tab.PT$Nom)
  Clean.tab.PT <- Clean.tab.PT[!duplicated(Clean.tab.PT$Sim_name),]
  Clean.tab.PT$Nom <- gsub(" $|  |^ ", "", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- gsub(" sp$| undiff$| group$| sp\\.$", "", Clean.tab.PT$Nom)
  Clean.tab.PT$Nom <- paste0(toupper(substr(Clean.tab.PT$Nom, 1, 1)), substr(Clean.tab.PT$Nom, 2, nchar(Clean.tab.PT$Nom)))
  Clean.tab.PT <- subset(Clean.tab.PT, select = - c(Sim_name))
  Clean.tab.PT <- Clean.tab.PT[!duplicated(Clean.tab.PT$Nom),]
  
  Clean.tab.PT <- Clean.tab.PT[order(Clean.tab.PT$Nom),]
  Clean.tab.PT <- Clean.tab.PT[!grepl("DELETE", Clean.tab.PT$Nom),]
  # Clean.tab.PT$TYPE <- paste0(toupper(substr(tolower(Clean.tab.PT$TYPE), 1, 1)), substr(tolower(Clean.tab.PT$TYPE), 2, nchar(tolower(Clean.tab.PT$TYPE))))
  Clean.tab.PT$TYPE <- paste("\\textbf{", Clean.tab.PT$TYPE, "}", sep = "")
  Clean.tab.PT <- Clean.tab.PT[nchar(Clean.tab.PT$Nom) >= 4,]
  Clean.tab.PT <- aggregate(Nom ~ TYPE, data = Clean.tab.PT, FUN = function(x) paste(x, collapse = ", "))
  
  #### Save and export ####
  if(is.null(Save.path) == F){
    dir.create(file.path(dirname(Save.path), ""), showWarnings = F, recursive = T)
    if(Save.csv == T){
      Save.path.csv <- paste(Save.path, "csv", sep = ".")
      write.table(Clean.tab.PT, file = Save.path.csv, row.names=T, col.names=NA, sep=",", dec = ".")
    }
    
    if(Save.Latex == T){
      library(xtable)
      if(Add.lines == T){Position.line <- c(-1, seq(0, nrow(Clean.tab.PT)))}
      else{Position.line <- c(-1, 0, nrow(Clean.tab.PT))}
      
      
      Tlatex <- xtable(Clean.tab.PT, caption = LateX.caption, type = "latex", label = Latex.lab)
      print.xtable(Tlatex, file = Save.path.Latex, booktabs = T, include.rownames = F, comment = F,
                   caption.placement = "top", sanitize.text.function = function(x){x}, 
                   tabular.environment = "longtable", floating = F,
                   scale.box = '\textwidth',
                   hline.after = Position.line)
    }}
  
}





#### Application ACA full ####
ACA.2026 = F
if(ACA.2026 == T){
  #### Import metadata ####
  if(exists("Meta.data.raw") == F){
    # All the cores for the whole process
    Meta.data.raw <- data.frame(readODS::read_ods("Import/ACA/Site/Paleo_sequences/Site_full_list.ods", col_types = "cnniccccccccccnccccc", na = "NA"))
  }
  # Only cores OK in chrono
  # Meta.data.raw <- data.frame(read.table("Import/ACA/Site/Paleo_sequences/Site_full_list_clim_selected.csv", sep = "\t", header = T))
  # Meta.data.raw <- Meta.data.raw[c(1:17)]
  # Meta.data.raw$Acces[Meta.data.raw$Acces == "TRUE"] <- "T"
  # Meta.data.raw$Acces[Meta.data.raw$Acces == "FALSE"] <- "F"
  
  #### Settings ####
  Run.process = F
  Run.plots = F
  Recalculate.models = F
  Recalculate.models.only.news = F; Apply.local.DB = T; Remove.pollen.aqua = T; Display.errors = T
  Recalculate.Ensembles = F; DES.Ensemble = T; DES.CO = F
  Extract.Clim.WC = F; Combine.by.sites = T; Scaling = T
  Fast.plot.check = F
  Fuzzy.Cm = F; DTW.test = F
  Plot.stacked.by.groups = T; Stats.cluster.power = F; Plot.codebar = F; Map.clim.Hol = F; Interractive.map = F
  Param.clim = c("MAAT", "MAP")
  Filter.chrono.before = F
  Min.end = 7800; Min.begin = 3000; Min.Resolution = 1000; Limits = c(100, 12500)
  # Min.end = 9200; Min.begin = 500; Min.Resolution = 1000; Limits = c(100, 12000)
  # Min.end = 3800; Min.begin = 2500; Min.Resolution = 1000; Limits = c(200, 11700)
  Interpolation = T; Filtering = 600
  # Keep.DB = c("ALDS")
  # Keep.archive = c("Loess")
  # Keep.archive = c("Speleothem")
  # Remove.archive = c("Loess", "Ice core", "Speleothem")
  # Test.Subset = c(110,232,245,421)
  # Keep.proxy = c("Pollen")
  # Remove.proxy = c("Pollen")
  # Remove.proxy = c("GDGT", "Pollen", "Lake_level")
  # Keep.proxy = c("Lake_level")
  # Keep.proxy = c("GDGT")
  # Remove.area = c("Iran", "Touran")
  # Keep.area= c("Caucasus", "Iran", "West. TSH", "Caspian")
  # Keep.area= c("Caucasus", "Iran", "Caspian")
  # Keep.area= c("Caucasus", "Iran", "Caspian", "West. TSH")
  # Keep.area= c("Caucasus", "Iran", "Caspian", "West. TSH", "Touran", "East. TSH")
  # Keep.area= c("Caucasus", "Iran", "Caspian", "West. TSH", "Touran", "East. TSH", "Altai", "Kazakh Ste.")
  # Keep.area = c("Caspian")
  # Keep.area = c("West. TSH")
  # Keep.site = c("Kesang")
  # Keep.country = c("Pakistan")
  
  Save.working.figures = "Figures/ACA/SI/"
  Save.update.metadata = "Results/ACA/Metadata_ACA_sites.Rds"
  Save.missing.tax = "Results/ACA/Pollen_missing/"
  # SI.fig.Rds = c("Didachara")
  
  #### Whole process ####
  if(Run.process == T){
    #### Metadata process ####
    # Meta.data.ACA <- Meta.data.raw
    Meta.data.ACA <- Meta.data.raw[Meta.data.raw$Zone != "Out",]
    Meta.data <- Meta.data.ACA[which(Meta.data.ACA$Acces == "T" | Meta.data.ACA$Processed == "T"),]
    Meta.data <- Meta.data[which(Meta.data$Proxy != ""),]
    
    Meta.data <- tidyr::separate_rows(Meta.data, Proxy, sep = ",")
    Meta.data$Proxy <- gsub("\\s", "", Meta.data$Proxy)
    Meta.data$Proxy[grep("n-alkane", Meta.data$Proxy)] <- "n-alkane"
    Meta.data$Site <- gsub("\\s|-|’", "_", Meta.data$Site)
    Site.simple <- gsub("_", "", Meta.data$Site)
    Meta.data$ID <- paste(seq(1:nrow(Meta.data)), substr(Site.simple, 1, 3+min(nchar(Site.simple), na.rm = T)), substr(Meta.data$Proxy, 1, 4), sep = "|")
    Meta.data$ID <- gsub("Poll", "Pol", Meta.data$ID)
    Meta.data$ID <- gsub("\\|Lake", "\\|Lev", Meta.data$ID)
    Meta.data$ID <- gsub("\\|Pale", "\\|Loess", Meta.data$ID)
    
    Meta.data$Path.raw <- paste("Import", Meta.data$Country, Meta.data$Proxy, Meta.data$Site, sep = "/")
    Meta.data$Path.clean <- paste("Results", Meta.data$Country, Meta.data$Proxy, Meta.data$Site, sep = "/")
    Meta.data$Path.clean[Meta.data$Archive == "Speleothem"] <- paste("Results", Meta.data$Country[Meta.data$Archive == "Speleothem"], "Geoch/d18O", Meta.data$Site[Meta.data$Archive == "Speleothem"], sep = "/")
    Meta.data$Path.clean[Meta.data$Archive == "Ice core"] <- paste("Results", Meta.data$Country[Meta.data$Archive == "Ice core"], "Geoch/Ice_core", Meta.data$Site[Meta.data$Archive == "Ice core"], sep = "/")
    Meta.data$Path.clean <- gsub("Pollen", "Pollen/Cores", Meta.data$Path.clean)
    Meta.data$Path.clean <- gsub("Dust", "Geoch/Dust", Meta.data$Path.clean)
    Meta.data$Path.clean[grepl("LegacyPollen2.0", Meta.data$DB)] <- gsub("Cores", "Cores/LegacyPollen2.0", Meta.data$Path.clean[grepl("LegacyPollen2.0", Meta.data$DB)])
    Meta.data$Path.clean[grepl("EMBSeCBIO", Meta.data$DB)] <- gsub("Cores", "Cores/EMBSeCBIO", Meta.data$Path.clean[grepl("EMBSeCBIO", Meta.data$DB)])
    Meta.data$Path.clean <- gsub("Paleodust", "Geoch/Loess", Meta.data$Path.clean)
    Meta.data$Path.clean[Meta.data$Archive == "Loess"] <- gsub("d13C", "Geoch/Loess", Meta.data$Path.clean[Meta.data$Archive == "Loess"])
    
    Meta.data$Path.clean <- gsub("GDGT", "GDGT/Cores", Meta.data$Path.clean)
    Meta.data$Path.clean[Meta.data$Archive %in% c("Lake", 'Alluvial/Fluvial')] <- gsub("d18O", "Geoch/d18O", Meta.data$Path.clean[Meta.data$Archive %in% c("Lake", "Alluvial/Fluvial")])
    # Meta.data$Path.clean[Meta.data$Archive == "Lake"] <- gsub("Isotopes", "Geoch/d18O", Meta.data$Path.clean[Meta.data$Archive == "Lake"])
    
    Meta.data$Path.clean[Meta.data$DB == "ALDS"] <- paste("Results", Meta.data$Country[Meta.data$DB == "ALDS"], "Geoch/Loess", Meta.data$Site[Meta.data$DB == "ALDS"], sep = "/")
    Meta.data$Path.clean[Meta.data$Proxy == "Lake_level"] <- gsub("Lake_level", "Geomorpho/Lake_level", Meta.data$Path.clean[Meta.data$Proxy == "Lake_level"])
    Meta.data$Path.clean[Meta.data$Proxy %in% c("n-alkane", "n-alkane+d13C+d2H")] <- gsub("(n-alkane)|(+d13C+d2H)", "Geoch/Alkane", Meta.data$Path.clean[Meta.data$Proxy %in% c("n-alkane", "n-alkane+d13C+d2H")])
    
    Meta.data$Path.model <- paste("Results", Meta.data$Country, Meta.data$Proxy, Meta.data$Site, sep = "/")
    Meta.data$Path.model <- gsub("Pollen", "Pollen/Func_trans", Meta.data$Path.model)
    Meta.data$Path.model <- gsub("Dust", "Geoch/Dust", Meta.data$Path.model)
    Meta.data$Path.model[Meta.data$Archive == "Speleothem"] <- paste("Results", Meta.data$Country[Meta.data$Archive == "Speleothem"], "Geoch/d18O", Meta.data$Site[Meta.data$Archive == "Speleothem"], sep = "/")
    Meta.data$Path.model[Meta.data$Archive %in% c("Lake", "Alluvial/Fluvial")] <- gsub("d18O", "Geoch/d18O", Meta.data$Path.model[Meta.data$Archive %in% c("Lake", "Alluvial/Fluvial")])
    # Meta.data$Path.model[Meta.data$Archive == "Lake"] <- paste("Results", Meta.data$Country[Meta.data$Archive == "Lake"], "Geoch/d18O", Meta.data$Site[Meta.data$Archive == "Speleothem"], sep = "/")
    Meta.data$Path.model[Meta.data$Archive == "Ice core"] <- paste("Results", Meta.data$Country[Meta.data$Archive == "Ice core"], "Geoch/Ice_core", Meta.data$Site[Meta.data$Archive == "Ice core"], sep = "/")
    Meta.data$Path.model <- gsub("Paleodust", "Geoch/Loess", Meta.data$Path.model)
    Meta.data$Path.model[Meta.data$DB == "ALDS"] <- paste("Results", Meta.data$Country[Meta.data$DB == "ALDS"], "Geoch/Loess", Meta.data$Site[Meta.data$DB == "ALDS"], sep = "/")
    Meta.data$Path.model <- gsub("GDGT", "GDGT/Cores", Meta.data$Path.model)
    Meta.data$Path.model[Meta.data$Proxy == "Lake_level"] <- gsub("Lake_level", "Geomorpho/Lake_level", Meta.data$Path.model[Meta.data$Proxy == "Lake_level"])
    Meta.data$Path.model[Meta.data$Proxy %in% c("n-alkane", "n-alkane+d13C+d2H")] <- gsub("(n-alkane)|(+d13C+d2H)", "Geoch/Alkane", Meta.data$Path.model[Meta.data$Proxy %in% c("n-alkane", "n-alkane+d13C+d2H")])
    Meta.data$Path.model[Meta.data$Archive == "Loess"] <- gsub("d13C", "Geoch/Loess", Meta.data$Path.model[Meta.data$Archive == "Loess"])
    # Meta.data$Path.model[Meta.data$Archive == "Lake"] <- gsub("Isotopes", "Geoch/d18O", Meta.data$Path.model[Meta.data$Archive == "Lake"])
    # stop()
    
    #### Actual clim ####
    if(any(is.na(Meta.data$Latitude)) & Display.errors == T){
      print("**** The following sites don't have coordinates and have been removed from the analysis. ****")
      print(Meta.data$Site[is.na(Meta.data$Latitude)])}
    Meta.data <- data.frame(Meta.data[!is.na(Meta.data$Latitude),])
    
    if(Extract.Clim.WC == T){
      print("**** Extraction of elevation + actual param. clim from DEM and worldclim2. ****")
      Mcoord <- Meta.data
      row.names(Mcoord) <- Mcoord$ID
      Mclim <- Clim.param.extraction(M = Mcoord, 
                                     All.param = F, Aridity = T, Altitude = T, Season = F, MAF = F,
                                     Seasonality = F, Biome = F, Map.display = T, Clim.display = F, Land.cover = F,
                                     Save.path = "Import/ACA/Site/Paleo_sequences/Site_full_list_clim.csv",
                                     Save.plot = "Figures/ACA/Maps/Paleo_samples_sites.pdf", H = 1000, W = 1500)}
    if(Extract.Clim.WC == F){Mclim <- data.frame(read.csv(file="Import/ACA/Site/Paleo_sequences/Site_full_list_clim.csv",sep=",",dec=".",header=T,row.names=1))}
    
    if(identical(nrow(Meta.data), nrow(Mclim)) == F){
      print("**** You try to add DEM altitude from former extraction. It is not possible. You need to try 'Extract.Clim.WC = T' to solve the problem. ****")
    }
    else{
      Meta.data$Altitude_DEM <- Mclim$Altitude
      Meta.data$Altitude[is.na(Meta.data$Altitude)] <- Meta.data$Altitude_DEM[is.na(Meta.data$Altitude)]
      Meta.data$Elevation_classes <- ifelse(Meta.data$Altitude <= 1000, "Low-elevation",
                                            ifelse(Meta.data$Altitude <= 2000, "Mid-elevation", "High-elevation")) 
      
      Meta.data$Elevation_classes <- ifelse(Meta.data$Altitude <= 2000, "Low-elevation", "High-elevation")
      
      Meta.data <- cbind(Meta.data, Mclim[Param.clim])
      row.names(Meta.data) <- NULL
    }
    
    Meta.data$Proxy_classes <- Meta.data$Proxy
    Meta.data$Proxy_classes[Meta.data$Archive == "Speleothem"] <- Meta.data$Archive[Meta.data$Archive == "Speleothem"]
    Meta.data$Proxy_classes[Meta.data$Archive == "Ice core"] <- Meta.data$Archive[Meta.data$Archive == "Ice core"]
    
    #### Filter proxy / archive to keep or remove ####
    if(is.null(Test.Subset) == F){Meta.data <- Meta.data[Test.Subset,]}
    if(is.null(Keep.archive) == F){Meta.data <- Meta.data[Meta.data$Archive %in% Keep.archive,]}
    if(is.null(Keep.proxy) == F){Meta.data <- Meta.data[Meta.data$Proxy %in% Keep.proxy,]}
    if(is.null(Remove.proxy) == F){Meta.data <- Meta.data[!Meta.data$Proxy %in% Remove.proxy,]}
    if(is.null(Remove.archive) == F){Meta.data <- Meta.data[!Meta.data$Archive %in% Remove.archive,]}
    if(is.null(Keep.DB) == F){Meta.data <- Meta.data[Meta.data$DB %in% Keep.DB,]}
    if(is.null(Keep.area) == F){Meta.data <- Meta.data[Meta.data$Area %in% Keep.area,]}
    if(is.null(Remove.area) == F){Meta.data <- Meta.data[!Meta.data$Area %in% Remove.area,]}
    if(is.null(Remove.site) == F){Meta.data <- Meta.data[!Meta.data$Site %in% Remove.site,]}
    if(is.null(Keep.site) == F){Meta.data <- Meta.data[Meta.data$Site %in% Keep.site,]}
    if(is.null(Remove.country) == F){Meta.data <- Meta.data[!Meta.data$Country %in% Remove.country,]}
    if(is.null(Keep.country) == F){Meta.data <- Meta.data[Meta.data$Country %in% Keep.country,]}
    # Meta.data <- Meta.data[Meta.data$Archive == "Speleothem" | Meta.data$Proxy == "GDGT",]   # A VIRER
    # Meta.data <- Meta.data[Meta.data$Archive == "Speleothem",]   # A VIRER
    # Meta.data <- Meta.data[Meta.data$Proxy == "",]   # A VIRER
    # Meta.data <- Meta.data[Meta.data$Proxy != "Lake_level",]   # A VIRER
    # Meta.data <- Meta.data[Meta.data$Proxy != "Pollen",]   # A VIRER
    # Meta.data <- Meta.data[Meta.data$Proxy == "Pollen",]   # A VIRER
    # Meta.data <- Meta.data[grepl("LegacyPollen2.0", Meta.data$DB),]   # A VIRER
    # Meta.data <- Meta.data[grepl("GDGT", Meta.data$Proxy),]   # A VIRER
    
    
    #### Check path to data ####
    PB.path.clean <- Meta.data$Path.clean[dir.exists(Meta.data$Path.clean) == F]
    if(Display.errors == T & length(PB.path.clean) > 0){
      OK.path.clean <- Meta.data$Path.clean[dir.exists(Meta.data$Path.clean) == T]
      print("**** The following pathes do not exist. ****")
      print(PB.path.clean[order(PB.path.clean)])
      print(paste("****", round(length(PB.path.clean)/(length(PB.path.clean) + length(OK.path.clean)), digits = 2)*100, "% of the sites not found. ****"))
    }
    Meta.data <- Meta.data[dir.exists(Meta.data$Path.clean) == T,]
    
    #### Filter cores too short (before FT) ####
    if((is.null(Min.end) == F | is.null(Min.begin) == F) & Filter.chrono.before == T){
      Meta.data <- Meta.data[Meta.data$Begin < Min.begin,]
      Meta.data <- Meta.data[Meta.data$End > Min.end,]
    }
    # stop()
    #### Process FT ####
    for(i in 1:nrow(Meta.data)){
      #### Prep. data for each i loop ####
      ID <- Meta.data$ID[i]
      Path.clean <- Meta.data$Path.clean[i]
      Save.path <- Meta.data$Path.model[i]
      Proxy <- Meta.data$Proxy[i]
      Site <- Meta.data$Site[i]
      Area <- Meta.data$Area[i]
      Archive <- Meta.data$Archive[i]
      Country <- Meta.data$Country[i]
      Interpretation <- Meta.data$Param.clim[i]
      DB <- Meta.data$DB[i]
      
      #### Speleothems ####
      if(Archive == "Speleothem" & Recalculate.models == T){
        #### Import data ####
        if(Display.errors == T){print(paste(Site, " (", ID, ") ", "- FT calculation.", sep = ""))}
        Speleo.path <- list.files(path = Path.clean, pattern = "speleothem.Rds", full.names = T, all.files = T, recursive = F)
        Data.spe <- readRDS(Speleo.path)
        Data.spe[grep("d18O", names(Data.spe))] <- - Data.spe[grep("d18O", names(Data.spe))]
        Data.spe[grep("d13C", names(Data.spe))] <- - Data.spe[grep("d13C", names(Data.spe))]
        Data.spe[grep("Sr_Ca", names(Data.spe))] <- - Data.spe[grep("Sr_Ca", names(Data.spe))]
        
        #### Proxy to climate ####
        vec.I <- trimws(strsplit(Interpretation, "\\+")[[1]])
        vec.P <- trimws(strsplit(Proxy, "\\+")[[1]])
        
        if(length(vec.I) == 0){vec.I <- "MAP"}
        print(vec.I)
        if(length(vec.I) == 1){vec.I <- rep(vec.I, length(vec.P))}
        if(length(vec.I) == length(vec.P)){
          if(any(vec.I %in% c("MPCOQ"))){vec.I[vec.I %in% c("MPCOQ")] <- "MAP"}
          if(any(vec.I %in% c("MTCOQ"))){vec.I[vec.I %in% c("MTCOQ")] <- "MAAT"}
          # Interp <- setNames(data.frame(t(data.frame(vec.I))), vec.P)
          
          MAP.j <- which(vec.I == "MAP")
          Keep.for.MAP <- vec.P[MAP.j]
          MAAT.j <- which(vec.I == "MAAT")
          Keep.for.MAAT <- vec.P[MAP.j]
          
          # Keep.for.MAP <- colnames(Interp[names(Interp) == "MAP"])
          # print(Keep.for.MAP)
          # print(MAP.j)
          # print(vec.P)
        }
        
        
        #### PCA ####
        if(is.null(Save.working.figures) == F){
          Save.PCA <- paste(Save.working.figures, "PCA/Speleothem/PCA_", Site, ".pdf", sep = "")
          Save.Proxy <- paste(Save.working.figures, "Proxy/Speleothem/Proxy_", Site, ".pdf", sep = "")
        }
        else{Save.PCA <- NULL; Save.Proxy <- NULL}
        
        Data.spe <- Data.spe[, colSums(is.na(Data.spe)) < nrow(Data.spe)]
        Data.spe <- Data.spe[!grepl("MaxI|MinI", names(Data.spe))]
        Data.spe <- Data.spe[!grepl("Depth", names(Data.spe))]
        Data.spe <- na.omit(Data.spe)
        
        if(nrow(Data.spe) > 2){
          Data.spe <- PCA.components(Data.spe, Save.PCA = Save.PCA)
        }
        
        if(any("PC1" %in% names(Data.spe))){
          if(unique(vec.I) == "MAP"){names(Data.spe) <- gsub("PC1", "MAP", names(Data.spe))}
          if(unique(vec.I) == "MAAT"){names(Data.spe) <- gsub("PC1", "MAAT", names(Data.spe))}
          
        }
        else{
          if(unique(vec.I) == "MAP"){names(Data.spe) <- gsub("d18O", "MAP", names(Data.spe))}
          if(unique(vec.I) == "MAAT"){names(Data.spe) <- gsub("d18O", "MAAT", names(Data.spe))}
        }
        
        # print(head(Data.spe))
        saveRDS(Data.spe, gsub(".Rds", "_clim_param.Rds", Speleo.path))
        
      }
      
      #### Loess ####
      if(Archive == "Loess" & DB == "ALDS" & Recalculate.models == T){
        if(Display.errors == T){print(paste(Site, " (", ID, ") ", "- FT calculation.", sep = ""))}
        
        Data.loess <- NULL; Loess.path <- NULL
        Loess.path <- list.files(path = Path.clean, pattern = "_loess.Rds", full.names = T, all.files = T, recursive = F)
        if(length(Loess.path) > 0){Data.loess <- readRDS(Loess.path)}
        if(is.null(Data.loess) == F){
          
          if(is.null(Save.working.figures) == F){
            Save.PCA <- paste(Save.working.figures, "PCA/ALDS/PCA_", Site, ".pdf", sep = "")
            Save.Proxy <- paste(Save.working.figures, "Proxy/ALDS/Proxy_", Site, ".pdf", sep = "")
          }
          else{Save.PCA <- NULL; Save.Proxy <- NULL}
          
          if(all(is.na(colSums(Data.loess, na.rm = T)))){next}
          Data.loess <- Data.loess[!grepl("Age_max|Age_min|Depth", names(Data.loess))]
          Data.loess <- Data.loess[, colSums(is.na(Data.loess), na.rm = T) < nrow(Data.loess)]
          Data.loess <- na.omit(Data.loess)
          if(!"Age" %in% names(Data.loess)){next}
          if(nrow(Data.loess) == 0){next}
          Data.loess <- PCA.components(Data.loess, Save.PCA = Save.PCA)
          
          List.loess.PCA.i <- List.loess.PCA[row.names(List.loess.PCA) == Site,]
          # print(List.loess.PCA.i)
          if(length(List.loess.PCA.i) > 0){
            if(List.loess.PCA.i == "x"){
              if(Display.errors == T){print("**** PCA reversed ****")}
              Data.loess$PC1 <- - Data.loess$PC1
            }
          }
          
          if(is.null(Save.Proxy) == F){
            M <- melt(Data.loess[!names(Data.loess) %in% c("Depth", "Age_min", "Age_max")], id = "Age")
            p <- ggplot(M, aes(x = Age, y = value))+facet_grid(variable ~., scale = "free_y")+geom_point()+geom_line()
            ggsave(p, file = Save.Proxy, width = 600*0.026458333, height = 800*0.026458333, units = "cm", limitsize = F)
          }
          
          if(any("PC1" %in% names(Data.loess))){names(Data.loess) <- gsub("PC1", "MAP", names(Data.loess))}
          else{names(Data.loess) <- gsub("MS_low", "MAP", names(Data.loess))}
          
          saveRDS(Data.loess, gsub(".Rds", "_clim_param.Rds", Loess.path))
          
        }
      }
      
      #### Pollen ####
      if(Proxy == "Pollen"){
        #### No recalculation ####
        FT.path <- list.files(path = Path.clean, pattern = ".csv", full.names = T, all.files = T, recursive = F)
        
        if(length(FT.path) == 0){
          if(Display.errors == T){print(paste("**** No row pollen data found to calculate FT in ", Site, ". ****", sep = ""))}
          next}
        MA <- FT.path[grep("MA_", FT.path)]
        MP <- FT.path[grep("MP_", FT.path)]
        MP <- data.frame(read.csv(file = MP, sep=",",dec=".",header=T, stringsAsFactors = F, row.names = 1))
        MA <- data.frame(read.csv(file = MA, sep=",",dec=".",header=T, stringsAsFactors = F, row.names = 1))
        
        if(is.null(Save.working.figures) == F){Save.DP <- paste(Save.working.figures, "DP/DP_", Site, ".pdf", sep = "")}
        else{Save.DP <- NULL}
        
        if(!dir.exists(dirname(Save.path))){dir.create(dirname(Save.path), recursive = T)}
        
        #### Calculation missing FT ####
        if(Recalculate.models == T){
          #### Check correspondance paleo/surf ####
          if(Display.errors == T){print(paste(Site, " (", ID, ") ", "- FT calculation.", sep = ""))}
          if(any(names(MP) %in% unique(Tab.Cor.Name$Nom)) == F){
            if(Display.errors == T){print(paste("**** The pollen data of", Site, "are not in the good shape to calculate FTs. ****"))}
            next}
          MP.conv <- Fossil.MAT.prep(MP, MA, Type_MAT = Type_conversion, Displot = T,  Corresp_name = Tab.Cor.Name, 
                                     Show.message = F, Save.plot = Save.DP, H = 500, W = 2000, Seuil.p = 0.005, 
                                     Save.missing = paste(Save.missing.tax, Site, ".csv", sep = ""))
          
          if(Remove.pollen.aqua == T){
            List.pollen.aqua.i <- List.pollen.aqua[row.names(List.pollen.aqua) == Site,]
            if(nrow(List.pollen.aqua.i) > 0){
              if(Display.errors == T){print("**** Aquatic pollen removed. ****")}
              List.pollen.aqua.i <- names(List.pollen.aqua.i)[which(List.pollen.aqua.i == "x")]
              MP.conv <- MP.conv[!names(MP.conv) %in% List.pollen.aqua.i]
              MP.conv <- MP.conv/rowSums(MP.conv)
            }
          }
          saveRDS(MP.conv, paste(Path.clean, "/", Site, "_MP_conv.Rds", sep = ''))
          
          DB.cal <- gsub(".*\\.", "", names(List.FT.calib))
          Method <- gsub("\\..*", "", names(List.FT.calib))
          
          #### Loop on each DB calib ####
          for(j in unique(DB.cal)){
            WAPLS.j = List.FT.calib[intersect(which(Method == "WAPLS"), grep(j, names(List.FT.calib)))]
            MAT.j = List.FT.calib[intersect(which(Method == "MAT"), grep(j, names(List.FT.calib)))]
            BRT.j = List.FT.calib[intersect(which(Method == "BRT"), grep(j, names(List.FT.calib)))]
            RF.j = List.FT.calib[intersect(which(Method == "RF"), grep(j, names(List.FT.calib)))]
            
            if(is.null(WAPLS.j) == F){MSURF <- data.frame(WAPLS.j[[1]][[1]]$y)}
            else{
              if(is.null(MAT.j) == F){MSURF <- data.frame(MAT.j[[1]][[1]]$y)}
              else{
                print("**** La j'ai un PB: il faut modifier l'entrainement de BRT et RF pour mettre la matrice de surface à l'intérieur ****")
                next}}
            
            if(Recalculate.models.only.news == T){
              FT.ENS.path <- list.files(path = Save.path, pattern = ".Rds", full.names = T, all.files = T, recursive = F)
              if(any(grepl(j, FT.ENS.path))){next}}
            
            if(Apply.local.DB == T){
              if(unique(DB.cal) == "CAUCDB"){if(! Area %in% c("Caucasus", "Caspian")){next}}
              if(unique(DB.cal) == "COSTDB"){if(! Area %in% c("Caucasus", "Baikal", "Altai", "Caspian", "Iran")){next}}
              if(unique(DB.cal) == "TUSDB"){if(! Country %in% c("Uzbekistan", "Kirghizistan", "Kazakhstan", "Iran")){next}}
              if(unique(DB.cal) == "TAIGDB"){if(! Area %in% c("Altai", "Kazakh Ste.", "Baikal", "MP")){next}}
            }
            
            MP.conv.j <- Fossil.corresp.surface(MP.conv, MSURF)
            if(is.null(WAPLS.j) == F | is.null(MAT.j) == F | is.null(RF.j) == F | is.null(BRT.j) == F){print("zozo")
              FT.core(Model.WAPLS = WAPLS.j, Model.MAT = MAT.j, Model.RF = RF.j, Model.BRT = BRT.j, MCore = MP.conv.j, MAge = MA$Age,
                      Save.RDS = T, Displot = F, Nodeparse = T, #Model.name = names(List.FT.calib[j]), #GDGT = T, GDGT.model = Method.j,
                      # Select.clim = Param.clim,
                      H = 500, W = 500, LakeName = Site,
                      Save.path = paste(Save.path, "/", Site, ".csv", sep = ''))
            }
          }}
      }
    }
    #### Check pathes to models ####
    PB.path.model <- Meta.data$Path.model[dir.exists(Meta.data$Path.model) == F]
    if(Display.errors == T & length(PB.path.model) > 0){
      OK.path.model <- Meta.data$Path.model[dir.exists(Meta.data$Path.model) == T]
      print("**** The following pathes do not exist for ensemble. ****")
      print(PB.path.model)
      print(paste("****", round(length(PB.path.model)/(length(PB.path.model) + length(OK.path.model)), digits = 2)*100, "% of the sites not found for ensemble. ****"))
    }
    Meta.data <- Meta.data[dir.exists(Meta.data$Path.model) == T,]
    
    #### Process Ensemble ####
    List.Ensemble <- list()
    for(i in 1:nrow(Meta.data)){
      #### Settings ####
      ID <- Meta.data$ID[i]
      Path.model <- Meta.data$Path.model[i]
      Path.clean <- Meta.data$Path.clean[i]
      Proxy <- Meta.data$Proxy[i]
      Site <- Meta.data$Site[i]
      Area <- Meta.data$Area[i]
      Archive <- Meta.data$Archive[i]
      DB <- Meta.data$DB[i]
      Surf.clim.i <- Meta.data[i, names(Meta.data) %in% Param.clim]
      row.names(Surf.clim.i) <- Site
      
      #### GDGT ####
      if(Proxy == "GDGT"){
        if(Display.errors == T){print(paste(ID, "-", Site, "- Ensemble"))}
        
        GDGT.path <- list.files(path = Path.model, pattern = ".Rds", full.names = T, all.files = T, recursive = F)
        if(any(grepl("full_models_CB_LS_ens.Rds", GDGT.path))){
          GDGT.path.mn <- GDGT.path[grepl("full_models_CB_LS_ens.Rds", GDGT.path)]
          GDGT.i <- readRDS(GDGT.path.mn)
          GDGT.i <- GDGT.i[c("Age", "MAAT_Ensemble", "MPCOQ_Ensemble")]
          
          if(any(grepl("full_models_CB_LS_ens_MaxI.Rds", GDGT.path))){
            GDGT.path.max <- GDGT.path[grepl("full_models_CB_LS_ens_MaxI.Rds", GDGT.path)]
            GDGT.path.min <- GDGT.path[grepl("full_models_CB_LS_ens_MinI.Rds", GDGT.path)]
            GDGT.i.max <- readRDS(GDGT.path.max)
            GDGT.i.min <- readRDS(GDGT.path.min)
            GDGT.i.max <- GDGT.i.max[c("MAAT_Ensemble", "MPCOQ_Ensemble")]
            GDGT.i.min <- GDGT.i.min[c("MAAT_Ensemble", "MPCOQ_Ensemble")]
            names(GDGT.i.max) <- c("MAAT_Ensemble_MaxI", "MPCOQ_Ensemble_MaxI")
            names(GDGT.i.min) <- c("MAAT_Ensemble_MinI", "MPCOQ_Ensemble_MinI")
            GDGT.i <- cbind(GDGT.i, GDGT.i.max, GDGT.i.min)
            GDGT.i <- GDGT.i[c(1,2,4,6,3,5,7)]
          }
        }
        else{
          if(any(grepl("full_models_CB_LS.Rds", GDGT.path))){
            GDGT.path <- GDGT.path[grepl("full_models_CB_LS.Rds", GDGT.path)]
            GDGT.i <- readRDS(GDGT.path)
            
            GDGT.i <- GDGT.i[c("Age", "MAAT_BRT_ACA", "MPCOQ_BRT_ACA")]
            
            if(any(grepl("full_models_CB_LS_ens_MaxI.Rds", GDGT.path))){
              GDGT.path.max <- GDGT.path[grepl("full_models_CB_LS_ens_MaxI.Rds", GDGT.path)]
              GDGT.path.min <- GDGT.path[grepl("full_models_CB_LS_ens_MinI.Rds", GDGT.path)]
              GDGT.i.max <- readRDS(GDGT.path.max)
              GDGT.i.min <- readRDS(GDGT.path.min)
              GDGT.i.max <- GDGT.i.max[c("MAAT_BRT_ACA", "MPCOQ_BRT_ACA")]
              GDGT.i.min <- GDGT.i.min[c("MAAT_BRT_ACA", "MPCOQ_BRT_ACA")]
              names(GDGT.i.max) <- c("MAAT_Ensemble_MaxI", "MPCOQ_Ensemble_MaxI")
              names(GDGT.i.min) <- c("MAAT_Ensemble_MinI", "MPCOQ_Ensemble_MinI")
              GDGT.i <- cbind(GDGT.i, GDGT.i.max, GDGT.i.min)
              GDGT.i <- GDGT.i[c(1,2,4,6,3,5,7)]
            }
          }
          else{
            GDGT.path <- GDGT.path[grepl("orig_cal.Rds", GDGT.path)]
            GDGT.i <- readRDS(GDGT.path)
            GDGT.i <- GDGT.i[c("Age", "MAAT")]
          }
          
          
        }
        if(ncol(GDGT.i) == 2){names(GDGT.i) <- c("Age", "MAAT_FT_Ensemble")}
        if(ncol(GDGT.i) == 3){names(GDGT.i) <- c("Age", "MAAT_FT_Ensemble", "MAP_FT_Ensemble")}
        if(ncol(GDGT.i) == 7){names(GDGT.i) <- c("Age", "MAAT_FT_Ensemble", "MAAT_FT_Ensemble_MaxI", "MAAT_FT_Ensemble_MinI", "MAP_FT_Ensemble", "MAP_FT_Ensemble_MaxI", "MAP_FT_Ensemble_MinI")}
        
        List.Ensemble[[i]] <- GDGT.i
        names(List.Ensemble)[i] <- ID
      }
      
      #### Speleothems ####
      if(Archive == "Speleothem"){
        # if(Display.errors == T){print(paste(ID, "-", Site, "- Ensemble"))}
        
        Speleo.path <- list.files(path = Path.clean, pattern = "speleothem_clim_param.Rds", full.names = T, all.files = T, recursive = F)
        Data.spe <- readRDS(Speleo.path)
        
        names(Data.spe)[names(Data.spe) == "MAP"] <- "MAP_FT_Ensemble"
        names(Data.spe)[names(Data.spe) == "MAP_MaxI"] <- "MAP_FT_Ensemble_MaxI"
        names(Data.spe)[names(Data.spe) == "MAP_MinI"] <- "MAP_FT_Ensemble_MinI"
        names(Data.spe)[names(Data.spe) == "MAAT"] <- "MAAT_FT_Ensemble"
        names(Data.spe)[names(Data.spe) == "MAAT_MaxI"] <- "MAAT_FT_Ensemble_MaxI"
        names(Data.spe)[names(Data.spe) == "MAAT_MinI"] <- "MAAT_FT_Ensemble_MinI"
        
        if(length(names(Data.spe)[names(Data.spe) == "MAAT_FT_Ensemble"]) == 0){Data.spe$MAAT_FT_Ensemble <- NA}
        Data.spe <- Data.spe[names(Data.spe) %in% c("Age", "MAP_FT_Ensemble", "MAP_FT_Ensemble_MaxI", "MAP_FT_Ensemble_MinI",
                                                    "MAAT_FT_Ensemble", "MAAT_FT_Ensemble_MaxI", "MAAT_FT_Ensemble_MinI")]
        
        Data.spe <- Data.spe[!is.na(Data.spe$Age),]
        Data.spe <- Data.spe[!duplicated(Data.spe$Age),]
        
        List.Ensemble[[i]] <- Data.spe
        names(List.Ensemble)[i] <- ID
      }
      
      #### Ice cores ####
      if(Archive == "Ice core"){
        # if(Display.errors == T){print(paste(ID, "-", Site, "- Ensemble"))}
        Data.ic <- NULL ; ic.path <- NULL
        ic.path <- list.files(path = Path.clean, pattern = "_ice_core.Rds", full.names = T, all.files = T, recursive = F)
        # print(ic.path)
        if(length(ic.path) > 0){Data.ic <- readRDS(ic.path)}
        if(is.null(Data.ic) == F){
          # print(Data.ic)
          Data.ic <- Data.ic[!is.na(Data.ic$Age),]
          Data.ic <- Data.ic[!duplicated(Data.ic$Age),]
          
          List.Ensemble[[i]] <- Data.ic
          names(List.Ensemble)[i] <- ID
        }
      }
      
      #### n-alkanes + isotopie lacustre ####
      if(Proxy == "n-alkane"){
        # if(Display.errors == T){print(paste(ID, "-", Site, "- Ensemble"))}
        Data.nal <- NULL ; ic.path <- NULL
        ic.path <- list.files(path = Path.clean, pattern = "_alkane.Rds", full.names = T, all.files = T, recursive = F)
        # print(Path.clean)
        # print(ic.path)
        if(length(ic.path) > 0){Data.nal <- readRDS(ic.path)}
        if(is.null(Data.nal) == F){
          # print(head(Data.nal))
          Data.nal <- Data.nal[!is.na(Data.nal$Age),]
          Data.nal <- Data.nal[!duplicated(Data.nal$Age),]
          
          List.Ensemble[[i]] <- Data.nal
          names(List.Ensemble)[i] <- ID
        }
      }
      
      #### Bulk/ostracode isotopie lacustre ####
      if(Proxy == "d18O" & Archive %in% c("Lake", "Peatland", "Alluvial/Fluvial")){
        if(Display.errors == T){print(paste(ID, "-", Site, "- Ensemble"))}
        Data.lkiso <- NULL ; lkiso.path <- NULL
        lkiso.path <- list.files(path = Path.clean, pattern = "_dO18.Rds", full.names = T, all.files = T, recursive = F)
        # print(Path.clean)
        # print(lkiso.path)
        if(length(lkiso.path) > 0){Data.lkiso <- readRDS(lkiso.path)}
        if(is.null(Data.lkiso) == F){
          # print(head(Data.lkiso))
          Data.lkiso <- Data.lkiso[!is.na(Data.lkiso$Age),]
          Data.lkiso <- Data.lkiso[!duplicated(Data.lkiso$Age),]
          
          List.Ensemble[[i]] <- Data.lkiso
          names(List.Ensemble)[i] <- ID
        }
      }
      
      #### Bulk/ostracode isotopie lacustre ####
      if(Proxy == "Dust" & Archive %in% c("Lake", "Peatland", "Alluvial/Fluvial")){
        if(Display.errors == T){print(paste(ID, "-", Site, "- Ensemble"))}
        Data.lkiso <- NULL ; lkiso.path <- NULL
        lkiso.path <- list.files(path = Path.clean, pattern = "(_XRF)|(_MS)\\.Rds", full.names = T, all.files = T, recursive = F)
        # print(Path.clean)
        # print(lkiso.path)
        if(length(lkiso.path) > 0){Data.lkiso <- readRDS(lkiso.path)}
        if(is.null(Data.lkiso) == F){
          # print(head(Data.lkiso))
          Data.lkiso <- Data.lkiso[!is.na(Data.lkiso$Age),]
          Data.lkiso <- Data.lkiso[!duplicated(Data.lkiso$Age),]
          
          List.Ensemble[[i]] <- Data.lkiso
          names(List.Ensemble)[i] <- ID
        }
      }
      
      #### Loess ####
      if(Archive == "Loess"){
        if(Display.errors == T){print(paste(ID, "-", Site, "- Ensemble"))}
        Data.loess <- NULL ; Loess.path <- NULL
        if(Proxy == "Paleodust"){Loess.path <- list.files(path = Path.model, pattern = "_loess_clim_param.Rds", full.names = T, all.files = T, recursive = F)}
        if(Proxy == "d13C"){Loess.path <- list.files(path = Path.model, pattern = "_d13C.Rds", full.names = T, all.files = T, recursive = F)}
        
        if(length(Loess.path) > 0){Data.loess <- readRDS(Loess.path)}
        if(is.null(Data.loess) == F){
          
          names(Data.loess)[names(Data.loess) == "MAP"] <- "MAP_FT_Ensemble"
          names(Data.loess)[names(Data.loess) == "MAP_MaxI"] <- "MAP_FT_Ensemble_MaxI"
          names(Data.loess)[names(Data.loess) == "MAP_MinI"] <- "MAP_FT_Ensemble_MinI"
          names(Data.loess)[names(Data.loess) == "MAAT"] <- "MAAT_FT_Ensemble"
          names(Data.loess)[names(Data.loess) == "MAAT_MaxI"] <- "MAAT_FT_Ensemble_MaxI"
          names(Data.loess)[names(Data.loess) == "MAAT_MinI"] <- "MAAT_FT_Ensemble_MinI"
          
          if(length(names(Data.loess)[names(Data.loess) == "MAAT_FT_Ensemble"]) == 0){Data.loess$MAAT_FT_Ensemble <- NA}
          Data.loess <- Data.loess[names(Data.loess) %in% c("Age", "MAP_FT_Ensemble", "MAP_FT_Ensemble_MaxI", "MAP_FT_Ensemble_MinI",
                                                            "MAAT_FT_Ensemble", "MAAT_FT_Ensemble_MaxI", "MAAT_FT_Ensemble_MinI")]
          Data.loess <- Data.loess[!is.na(Data.loess$Age),]
          Data.loess <- Data.loess[!duplicated(Data.loess$Age),]
          
          # print(head(Data.loess))
          List.Ensemble[[i]] <- Data.loess
          names(List.Ensemble)[i] <- ID
        }
      }
      
      #### Lake levels ####
      if(Proxy == "Lake_level"){
        # if(Display.errors == T){print(paste(ID, "-", Site, "- Ensemble"))}
        
        LL.path <- list.files(path = Path.model, pattern = "_lake_level.Rds", full.names = T, all.files = T, recursive = F)
        Data.ll <- data.frame(readRDS(LL.path))
        Data.ll <- Data.ll[!is.na(Data.ll$MAP),]
        if(length(unique(Data.ll$MAP)) == 1){next}
        names(Data.ll)[names(Data.ll) == "MAP"] <- "MAP_FT_Ensemble"
        Data.ll$MAAT_FT_Ensemble <- NA
        Data.ll <- Data.ll[c(1,3,2)]
        List.Ensemble[[i]] <- Data.ll
        names(List.Ensemble)[i] <- ID
      }
      
      #### Pollen ####
      if(Proxy == "Pollen"){
        #### Surface values is missing ? ####
        if(length(Surf.clim.i) == 0){
          print("**** Ensemble calculation impossible, since no climate values are in 'Meta.data', please rerun 'Extract.Clim.WC = F'. ****")
          stop()
        }
        
        #### Try do not have to recalculate ####
        FT.ENS.path <- list.files(path = Path.model, pattern = "Ensemble_FT.Rds", full.names = T, all.files = T, recursive = F)
        if(Recalculate.Ensembles == F & length(FT.ENS.path) == 1){
          FT.Ens.res  <- readRDS(FT.ENS.path)
          if(all(sapply(Param.clim, function(a) any(grepl(a, names(FT.Ens.res)))))){
            Missing.param = F
          }
          else{
            Missing.param = T
          }
        }
        else{Missing.param = F}
        
        #### Recalculate and erase elder Ensemble ####
        if(Recalculate.Ensembles == T | length(FT.ENS.path) == 0 | Missing.param == T){
          files.FT  <- list.files(path = Path.model, pattern = ".Rds", full.names = T, all.files = T, recursive = F)
          if((length(FT.ENS.path) == 0 | Missing.param == T) & Display.errors == T){
            if(length(files.FT) > 0){
              print(paste("**** Since not ensemble was found, we recalculated the ensemble model for", Site, "with", Proxy, "****"))
            }
            
            else{
              print(paste("**** Since not ensemble was found, but we cannot recalculated ensemble for", Site, "with", Proxy, " since no model is found. ****"))
              List.Ensemble[[i]] <- NULL
              next
            }
          }
          if(Display.errors == T){print(paste(ID, "-", Site, "- Ensemble"))}
          
          if(is.null(Save.working.figures) == F){Save.DP <- paste(Save.working.figures, "FT_Ensemble/", Area, "/Ens_", Site, ".pdf", sep = "")}
          else{Save.DP <- NULL}
          
          files.FT  <- files.FT[!grepl("_full|Ensemble|randomFT", files.FT)]
          List.FT.paleo <- do.call(c, lapply(files.FT, readRDS))
          Keep.calib <- intersect(names(List.FT.paleo), names(List.FT.calib))
          Keep.calib <- intersect(names(List.FT.paleo), names(List.FT.calib))
          # print(Keep.calib)
          
          if(length(List.FT.paleo[Keep.calib]) == 0){
            print(paste("**** The calibration", unique(gsub(".*\\.", "", Keep.calib)), " was not found for", Site, "with", Proxy, " then, no Ensemble is calculated. ****"))
            next
          }
          
          #### If Dynamic Ensemble selection is TRUE ####
          if(DES.Ensemble == T){
            DES.Ensemble.i = T
            MP.conv.path <- list.files(path = Path.clean, pattern = "conv.Rds", full.names = T, all.files = T, recursive = F)
            
            if(length(MP.conv.path) == 0){
              FT.path <- list.files(path = Path.clean, pattern = ".csv", full.names = T, all.files = T, recursive = F)
              
              if(length(FT.path) == 0){
                print(paste("**** DES impossible, since no row pollen data found for ", Site, ". ****", sep = ""))
                DES.Ensemble.i = F
                next
              }
              MA <- FT.path[grep("MA_", FT.path)]
              MP <- FT.path[grep("MP_", FT.path)]
              MP <- data.frame(read.csv(file = MP, sep=",",dec=".",header=T, stringsAsFactors = F, row.names = 1))
              MA <- data.frame(read.csv(file = MA, sep=",",dec=".",header=T, stringsAsFactors = F, row.names = 1))
              
              if(any(names(MP) %in% unique(Tab.Cor.Name$Nom)) == F){
                if(Display.errors == T){print(paste("**** The pollen data of", Site, "are not in the good shape to calculate FTs. ****"))}
                next}
              MP.conv <- Fossil.MAT.prep(MP, MA, Type_MAT = Type_conversion, Displot = T,  Corresp_name = Tab.Cor.Name, Show.message = Display.errors)
              saveRDS(MP.conv, paste(Path.clean, "/", Site, "_MP_conv.Rds", sep = ''))
            }
            else{
              MP.conv <- readRDS(MP.conv.path)
            }
            
            if(DES.CO == T){
              CO.i <- 1/length(unique(gsub(".*\\.", "", Keep.calib)))
            }
            else{CO.i <- NULL}
          }
          else{MP.conv <- NULL; DES.Ensemble.i = F}
          
          #### Calculate ensemble ####
          FT.Ens.res <- Pollen.Ensemble.RMSE(MPal = List.FT.paleo[Keep.calib],
                                             Msurf = List.FT.calib[Keep.calib],
                                             Surf.val = Surf.clim.i, Surf.val.max.age = 1000, Mcalib = Calib.list, DES.Ensemble = DES.Ensemble.i, Mpal.raw = MP.conv,
                                             Show.DES.contrib = T, Display.sim = T, DES.Soft.max = 9, DES.distance.metric = "aitchison", DES.KDE = T, DES.cut.off = CO.i,
                                             Save.plot = Save.DP, H = 500, W = 1200, Clim.lab = Param.clim[which(Param.clim %in% names(List.FT.paleo[Keep.calib][[1]]))],
                                             Save.path = paste(Path.model, "/", Site, "_Ensemble_FT.Rds", sep = ""), Name.core = Site, Dot.size = 2,
                                             Plot.x = "Age", Keep.clim = Param.clim[which(Param.clim %in% names(List.FT.paleo[Keep.calib][[1]]))])
        }
        List.Ensemble[[i]] <- FT.Ens.res
        names(List.Ensemble)[i] <- ID
      }
    }
    
    List.Ensemble <- Filter(Negate(is.null), List.Ensemble)
    List.Ensemble <- lapply(List.Ensemble, function(df) {df[, names(df) != "x", drop = F]})
    all_columns <- unique(unlist(lapply(List.Ensemble, names)))
    List.Ensemble <- lapply(List.Ensemble, function(df) {
      if(nrow(df) > 0){
        missing_cols <- setdiff(all_columns, names(df))
        for (col in missing_cols) {df[[col]] <- NA}
        df <- df[all_columns]
        return(df)
      }
    })
    
    Meta.data <- Meta.data[Meta.data$ID %in% names(List.Ensemble),]
    if(Display.errors == T){
      Missing.sites <- Meta.data.ACA[Meta.data.ACA$Proxy == "Pollen",]
      print(paste("**** Total number of pollen sites: ", nrow(Missing.sites)))
      
      print(paste("**** The Ensemble models are ready for", nrow(Meta.data), "sites. ****"))
      print(paste("**** For pollen", nrow(Meta.data[Meta.data$Proxy == "Pollen",]), "sites. ****"))
      print(paste("**** For GDGT", nrow(Meta.data[Meta.data$Proxy == "GDGT",]), "sites. ****"))
      print(paste("**** For speleothems", nrow(Meta.data[Meta.data$Archive == "Speleothem",]), "sites. ****"))
      print(paste("**** For Lake level", nrow(Meta.data[Meta.data$Proxy == "Lake_level",]), "sites. ****"))
      print(paste("**** For Loess dust", nrow(Meta.data[Meta.data$DB == "ALDS",]), "sites. ****"))
    }
    
    #### Combine by sites (e.g. Issyk-Kul...) ####
    if(Combine.by.sites == T){
      Keep.df <- List.Ensemble[Meta.data$ID[Meta.data$Main_site == ""]]
      Comb.df <- List.Ensemble[Meta.data$ID[Meta.data$Main_site != ""]]
      Meta.data.cb <- Meta.data[Meta.data$Main_site != "",]
      
      Comb.df <- lapply(split(Meta.data.cb$ID, Meta.data.cb$Main_site), function(sites_in_area) {
        dfs <- Comb.df[sites_in_area]
        dfs <- dfs[!sapply(dfs, is.null)]
        if (length(dfs) == 0) return(NULL)
        out <- bind_rows(dfs, .id = "SiteName")
        if ("Age" %in% colnames(out)) {out <- arrange(out, Age)}
        out$SiteName <- NULL
        out
      })
      
      grouped_sites <- split(Meta.data.cb$ID, Meta.data.cb$Main_site)
      names(Comb.df) <- sapply(grouped_sites, function(x) x[1])
      List.Ensemble <- c(Keep.df, Comb.df)
      num_order <- as.numeric(sub("\\|.*", "", names(List.Ensemble)))
      List.Ensemble <- List.Ensemble[order(num_order)]
      List.Ensemble <- Filter(Negate(is.null), List.Ensemble)
      
      Meta.data <- Meta.data[!Meta.data$ID %in% Meta.data.cb$ID,]
      Meta.data.cb <- Aggregate.metadata(Meta.data.cb, group_col  = "Main_site", site_col   = "Site",
                                         spatial_cols = c("Longitude", "Latitude"),
                                         mean_cols    = c("Altitude", "Altitude_DEM", "Res", "MAAT", "MAP"),
                                         range_cols   = c("Begin" = "min", "End" = "max"),
                                         first_cols   = c("ID", "Country", "Zone", "Area", "Archive", "Acces", "Processed", "Path.raw", "Path.clean", "Path.model", "Coerence", "Elevation_classes"),
                                         concat_cols  = c("Reference", "Proxy_classes", "Proxy", "DB", "Site_full", "Param.clim"), collapse_sep = (", "))
      
      Meta.data <- rbind(Meta.data, Meta.data.cb)
      Meta.data <- Meta.data[order(as.double(row.names(Meta.data))),]
      Meta.data <- Meta.data[Meta.data$ID %in% names(List.Ensemble),]
    }
    
    #### Filter cores too short (after FT) ####
    List.Ensemble.Filtered <- List.Ensemble
    Time.boun <- sapply(List.Ensemble.Filtered, function(x) c(round(min(x$Age, na.rm = T), digits = -2), round(max(x$Age, na.rm = T), digits = -2)))
    Nsamp <- sapply(List.Ensemble.Filtered, nrow)
    Meta.data$Begin <- Time.boun[1,]
    Meta.data$End <- Time.boun[2,]
    Extend <- Meta.data$End - Meta.data$Begin
    Meta.data$Res <- round(Extend/Nsamp, digits = -1)
    Meta.data.no.filter <- Meta.data
    
    if((is.null(Min.end) == F | is.null(Min.begin) == F) & Filter.chrono.before == F){
      if(Display.errors == T){print("**** Filtering ****")}
      if(is.null(Min.Resolution) == F){
        Keep.highres.cores <- Meta.data$ID[Meta.data$Res <= Min.Resolution]
        Meta.data <- Meta.data[Meta.data$ID %in% Keep.highres.cores,]
        List.Ensemble.Filtered <- List.Ensemble.Filtered[Keep.highres.cores]
      }
      
      if(is.null(Min.begin) == T){Min.begin <- min(Time.boun[1,], na.rm = T)}
      if(is.null(Min.end) == T){Min.end <- min(Time.boun[2,], na.rm = T)}
      Keep.long.cores <- Meta.data$ID[Meta.data$End >= Min.end & Meta.data$Begin <= Min.begin]
      List.Ensemble.Filtered <- List.Ensemble.Filtered[Keep.long.cores]
      List.Ensemble.Filtered <- Filter(Negate(is.null), List.Ensemble.Filtered)
      Meta.data <- Meta.data[Meta.data$ID %in% Keep.long.cores,]
    }
    else{List.Ensemble.Filtered <- List.Ensemble}
    
    if(Display.errors == T){
      print(paste("**** After filtering, n-pollen =", nrow(Meta.data[Meta.data$Proxy == "Pollen",]), "sites. ****"))
      print(paste("**** After filtering, n-GDGT =", nrow(Meta.data[Meta.data$Proxy == "GDGT",]), "sites. ****"))
      print(paste("**** After filtering, n-speleothems =", nrow(Meta.data[Meta.data$Archive == "Speleothem",]), "sites. ****"))
      print(paste("**** After filtering, n-Lake-levels =", nrow(Meta.data[Meta.data$Proxy == "Lake_level",]), "sites. ****"))
      print(paste("**** After filtering, n-Loess-dust =", nrow(Meta.data[Meta.data$DB == "ALDS",]), "sites. ****"))
    }
    
    #### Process Cleaning, interpolation, binning ####
    List.Ensemble.clean <- Stacking.quantif(
      Imput.list = List.Ensemble.Filtered,#[c(1,2)],
      Scaling = Scaling, Binning = F, Anomaly = F, Add.bin.count = F, Bin.sd = F, Detrending = F,
      Limits = Limits, Plot.x = "Age",
      Filter.cut.off = Filtering,
      Interp.step = 50,
      Return.list = T, Interpolation = Interpolation,
      Keep.clim = paste(Param.clim, "FT_Ensemble", sep = "_"), New.param.names = Param.clim
      # New.param.names = gsub("_FT_Ensemble", "", names(List.Ensemble.Filtered[[1]])[-1])
    )
    
    #### Save plot all filtered and clean data ####
    if(is.null(Save.working.figures) == F){
      Plot.test <- lapply(names(List.Ensemble.clean), function(name) {
        df <- List.Ensemble.clean[[name]]
        df$Group <- name
        return(df)
      })
      Plot.test <- bind_rows(Plot.test)
      p <- ggplot(Plot.test, aes(x = Age, y = MAAT, color = Group)) + geom_point() + geom_line() + facet_wrap(Group ~.) + theme(legend.position = "none")
      H = 2000; W = 2000
      ggsave(p, file = paste(Save.working.figures, "Multi-site/Ensemble_ACA.pdf", sep = ""), width = W*0.026458333, height = H*0.026458333, units = "cm", limitsize = F)
      
    }
    
    List.Ensemble.clean <- Filter(Negate(is.null), List.Ensemble.clean)
    LEC.Red <- lapply(names(List.Ensemble.clean), function(name) {
      df <- List.Ensemble.clean[[name]]
      # new_names <- paste(Meta.data$Site[Meta.data$ID == name][[1]], names(df), Meta.data$Proxy[Meta.data$ID == name][[1]], sep = ".") # pour code barre
      new_names <- paste(Meta.data$ID[Meta.data$ID == name][[1]], names(df), sep = ".") # pour Fuzzy C-means
      names(df) <- new_names
      names(df)[grep("Age", names(df))] <- "Age"
      return(df)
    })
    
    write.table(Meta.data, "Import/ACA/Site/Paleo_sequences/Site_full_list_clim_selected.csv", sep = "\t", row.names = F)
    write.table(Meta.data.no.filter, "Results/ACA/Site_full_no_filter.csv", sep = "\t", row.names = F)
    saveRDS(Meta.data, "Import/ACA/Site/Paleo_sequences/Site_full_list_clim_selected.Rds")
    saveRDS(Meta.data.no.filter, "Results/ACA/Site_full_no_filter.Rds")
    
    saveRDS(List.Ensemble, "Results/ACA/Bulk/List.Ensemble.Rds")
    saveRDS(List.Ensemble.Filtered, "Results/ACA/Bulk/List.Ensemble.Filtered.Rds")
    saveRDS(List.Ensemble.clean, "Results/ACA/Bulk/List.Ensemble.clean.Rds")
    saveRDS(LEC.Red, "Results/ACA/Bulk/LEC.Red.Rds")
    
    To.complete <- Meta.data[c("Site", "Longitude", "Latitude", "Begin", "End", "Res")]
    Meta.data.ACA <- suppressMessages(left_join(Meta.data.ACA, To.complete, by = c("Site",	"Longitude",	"Latitude")))
    Meta.data.tofill <- suppressMessages(left_join(Meta.data.raw, Meta.data, by = c("Site",	"Longitude",	"Latitude")))
    Meta.data.tofill <- Meta.data.tofill[order(Meta.data.tofill$Site),]
    # write.table(Meta.data.tofill, "Import/ACA/Site/Paleo_sequences/test.csv", sep = "\t", row.names = F)
    
    #### Fast plot checking ####
    if(Fast.plot.check == T){
      MP <- List.Ensemble.clean
      # MP <- List.Ensemble
      # Select <- c("502|Tonne|d18O", "106|Balik|Pol", "342|Kesan|d18O", "363|YE|Loess", "10|GS18|Pol", "276|LJW|GDGT")
      # MP <- MP[Select]
      MP <- lapply(names(MP), function(name) {
        df <- MP[[name]]
        df$Group <- name
        return(df)})
      MP <- bind_rows(MP)
      MP <- melt(MP, id = c("Age", "Group"))
      MP <- MP[!is.na(MP$value),]
      
      # MP <- MP[MP$Age < 12000,]
      p <- ggplot(MP, aes(x = Age, y = value))+facet_wrap(Group ~ variable, scale = "free_y")+geom_point()+geom_line()
      # devices <- dev.list()
      # if (is.null(devices) || !"X11cairo" %in% names(devices)) {X11()}
      print(p)
    }
    
    #### Export GitHub Test 1 ####
    saveRDS(Meta.data, "Project_GitHub/WCA_paleoclim/Import/Metadata/Site_full_list_clim_selected.Rds")
    # saveRDS(Meta.data.no.filter, "Project_GitHub/WCA_paleoclim/Import/Metadata/Site_full_no_filter.Rds")
    saveRDS(List.Ensemble, "Project_GitHub/WCA_paleoclim/Import/Data/List_Ensemble.Rds")
    # saveRDS(List.Ensemble.Filtered, "Project_GitHub/WCA_paleoclim/Import/Data/List_Ensemble_Filtered.Rds")
    # saveRDS(List.Ensemble.clean, "Project_GitHub/WCA_paleoclim/Import/Bulk/List_Ensemble_clean.Rds")
    # saveRDS(LEC.Red, "Project_GitHub/WCA_paleoclim/Import/Bulk/LEC_Red.Rds")
    
  }
  else{
    Meta.data <- readRDS("Import/ACA/Site/Paleo_sequences/Site_full_list_clim_selected.Rds")
    Meta.data.no.filter <- readRDS("Results/ACA/Site_full_no_filter.Rds")
    List.Ensemble <- readRDS("Results/ACA/Bulk/List.Ensemble.Rds")
    List.Ensemble.Filtered <- readRDS("Results/ACA/Bulk/List.Ensemble.Filtered.Rds")
    List.Ensemble.clean <- readRDS("Results/ACA/Bulk/List.Ensemble.clean.Rds")
    LEC.Red <- readRDS("Results/ACA/Bulk/LEC.Red.Rds")
  }
  
  
  
}
}
