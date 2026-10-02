#### Required packages ####
required_packages <- c(
  # Core
  "dplyr", "reshape", "reshape2", "zoo", "jmuOutlier", 
  # SIG packages
  "sf", "raster", "tidyr", "ggh4x", "spatstat", "terra",
  # Visualization
  "ggplot2", "ggnewscale", "corrplot", "ggcorrplot", "ggthemes",
  "patchwork", "RColorBrewer", "ggimage", "multcompView", "xtable", "maps",
  # Reproducibility
  "renv"
)

required_pak_github <- c("EarthSystemDiagnostics/corit")


#### First installation ####
First = F
if(First == T){
  # renv::init()
  # renv::status()
  
  renv::install(required_packages)
  renv::install(required_pak_github)
  
  missing <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1))]
  
  if (length(missing) > 0L) {stop("The following required packages are not installed: ", paste(missing, collapse = ", "))}
  
  renv::status()
  renv::snapshot()
  renv::snapshot(packages = required_packages, prompt = F)
  
  lock <- renv::lockfile_read()
  
  missing_lock <- setdiff(required_packages, names(lock$Packages))
  if(length(missing_lock) > 0){message("Some package were not properly recorded in the lock file.")}
}

required_packages <- c(required_packages, required_pak_github)

#### Install PKG ####
Install.PKG = F
if(Install.PKG == T){
  install_and_load <- function(pkgs) {
    for (pkg in pkgs) {
      if (!requireNamespace(pkg, quietly = TRUE)){
        message(paste0("Installing missing package: ", pkg))
        install.packages(pkg, dependencies = TRUE)}
      library(pkg, character.only = TRUE)
    }
  }
  
  # Run it
  install_and_load(required_packages)
  
}

#### Install by Renv ####
Install.Renv = F
if(Install.Renv == T){
  if (!requireNamespace("renv", quietly = TRUE)) {install.packages("renv")}
  renv::load()
  status <- renv::status()
  
  if (!isTRUE(status$synchronized)) {
    message("The renv environment is not synchronized with renv.lock.")
    message("Restoring the project environment...")
    renv::restore(prompt = FALSE)
    message("renv environment successfully restored.")
  }
  
  
  load_pak <- function(pkgs) {for (pkg in pkgs) {if(requireNamespace(pkg, quietly = TRUE)){library(pkg, character.only = TRUE)}}}
  
  load_pak(required_packages)
}

message("All required packages are installed and loaded.")
