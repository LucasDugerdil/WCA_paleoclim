# library("ggplot2"); library("reshape2"); library("ggnewscale"); library("patchwork")
# library("e1071"); library("factoextra"); library("dplyr"); library("zoo")
# library("dtwclust"); library("ggrepel")
# library("randomForest"); library("gbm")

required_packages <- c(
  # Core
  "dplyr", "reshape", "reshape2",
  # "tibble", "rioja", "tidypaleo", "randomForest", "gbm",
  # SIG packages
  "rgdal", "sf",
  # "stringr", "reshape2",
  # Visualization
  "ggplot2", "ggnewscale", #"RColorBrewer", "patchwork", "ggthemes", "ggpmisc",
  "patchwork",
  # Reproducibility
  "renv"
)

install_and_load <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste0("Installing missing package: ", pkg))
      install.packages(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

# Run it
install_and_load(required_packages)

message("All required packages are installed and loaded.")

