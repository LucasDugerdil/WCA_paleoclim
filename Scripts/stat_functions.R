Scaling.fun <- function(x = NULL, a = -1, b = 1){x <- a + (b-a)*(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)); return(x)}

Stacking.quantif <- function(Imput.list, Keep.clim = NULL, Scaling = F, Plot.x = "Age", Type.name = NULL, Return.list = F, Add.bin.count = F, Keep.empty.bin = F, Bin.sd = F,
                             Limits = NULL, Binning = F, Anomaly = F, Detrending = F, Windows.length = NULL, Stacking = T, Scale.min = -1, Scale.max = 1,
                             Interpolation = F, Interp.step = 10, Filter.cut.off = 500, New.param.names = NULL, Calib.name = NULL) {
  
  #### Inner functions ####
  # Scaling.fun <- function(x = NULL, a = -1, b = 1){x <- a + (b-a)*(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)); return(x)}
  
  scale_ensemble_df <- function(df, a = -1, b = 1) {
    detect_ensemble_params <- function(names_vec) {
      means <- names_vec[!grepl("_(MaxI|MinI)$", names_vec) & names_vec != Plot.x]
      params <- lapply(means, function(m) {c(mean = m, max  = paste0(m, "_MaxI"), min  = paste0(m, "_MinI"))})
      params[sapply(params, function(x) all(x %in% names_vec))]
    }
    
    Scaling.fun <- function(x = NULL, a = -1, b = 1){x <- a + (b-a)*(x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)); return(x)}
    
    scale_ensemble_triplet <- function(df, triplet, a = -1, b = 1) {
      ref_min <- min(df[, triplet], na.rm = TRUE)
      ref_max <- max(df[, triplet], na.rm = TRUE)
      
      if (!is.finite(ref_min) || ref_min == ref_max) return(df)
      
      scale_fun <- function(x) {a + (b - a) * (x - ref_min) / (ref_max - ref_min)}
      
      df[, triplet] <- lapply(df[, triplet], scale_fun)
      return(df)
    }
    
    triplets <- detect_ensemble_params(names(df))
    
    if(length(triplets) == 0){
      df <- cbind(Age = df[Plot.x], data.frame(apply(df[!names(df) %in% Plot.x], 2, Scaling.fun, a, b)))
    }
    else{
      for(tr in triplets){df <- scale_ensemble_triplet(df, tr, a, b)}}
    
    return(df)
  }
  
  Bin.by.TW <- function(x, Anomaly = F, Add.bin.count = F){
    Nb.param <- ncol(x)
    x <- dplyr::mutate(x, TimeWind = cut(.data[[Plot.x]], breaks = TW.x, include.lowest = F))
    
    if(Add.bin.count == T){
      counts <- as.data.frame(table(x$TimeWind))
      colnames(counts) <- c("Group.1", "N")
    }
    
    if(any(colSums(is.na(x[-c(1, ncol(x))])) == nrow(x))){
      Keep.order.col <- names(x) 
      Keep.col.na <- names(x[which(colSums(is.na(x)) == nrow(x))])
      x <- x[which(colSums(is.na(x)) != nrow(x))]
    }
    else{Keep.col.na <- NULL}
    
    if(Bin.sd == T){
      y <- aggregate(. ~ TimeWind, data = x[-c(1)], FUN = sd, na.rm = T, na.action = na.pass)
      names(y) <- c("Group.1", paste("SD", names(y)[-c(1)], sep = "_"))
    }
    
    x <- aggregate(x, list(x$TimeWind), FUN = mean, na.action = na.pass, na.rm = T)
    
    if(Keep.empty.bin == T){
      Add.missing.TW <- levels(x$Group.1)[which(!levels(x$Group.1) %in% unique(x$Group.1))]
      MNA <- data.frame(matrix(nrow = length(Add.missing.TW), ncol = ncol(x)))
      MNA[c(1)] <- Add.missing.TW
      names(MNA) <- names(x)
      x <- rbind(x, MNA)
    }
    
    if(Add.bin.count == T & Keep.empty.bin == F){counts <- counts[which(counts$Group.1 %in% unique(as.character(x$Group.1))),]}
    
    x[[Plot.x]] <- (as.numeric(gsub(",.*", "", gsub("\\]", "", gsub("\\(", "", x$Group.1))))+as.numeric(gsub(".*,", "", gsub("\\]", "", gsub("\\(", "", x$Group.1)))))/2
    x <- x[order(x[[Plot.x]]),]
    
    if(Bin.sd == T){
      x <- dplyr::full_join(x, y, by = "Group.1")
      
    }
    
    x <- x[!names(x) %in% c("TimeWind", "Group.1")]
    
    if(Anomaly == T){
      if(Bin.sd == F){x <- cbind(Age = x[[Plot.x]], data.frame(apply(x[2:Nb.param], 2, function(x) x <- x - x[1])))}
      else{
        Off.set <- as.vector(unlist(x[1,2:Nb.param]))
        x <- cbind(Age = x[[Plot.x]],
                   data.frame(apply(x[which(names(x) != "Age" & ! grepl("^SD_", names(x)))], 2, function(x) x <- x - x[1])),
                   x[grepl("^SD_", names(x))])
      }}
    
    if(Anomaly == F){x <- cbind(Age = x[[Plot.x]], x[-c(1)])}
    
    if(Add.bin.count == T){x <- cbind(x, N = counts$N)}
    
    if(is.null(Keep.col.na) == F){
      x[Keep.col.na] <- NA
      if(Bin.sd == T){
        x[paste("SD", Keep.col.na, sep = "_")] <- NA
      }
    }
    return(x)
  }
  
  
  Corit.simple <- function (timser, detr = Detrending, method = c("InterpolationMethod", "DirectFiltering", "IntegrandInterpolationMethod", "NoFilter"), 
                            appliedFilter = c("gauss", "runmean", "lowpass"), fc, tn = seq(from = min(index(timser)), to = max(index(timser)), by = 10), dt,
                            int.method = c("linear", "nearest"), k = 5){
    n <- max(index(timser))
    if(method == "InterpolationMethod"){res <- InterpolationMethod(detrTimser(timser, detr), fc, dt, n, int.method, appliedFilter, k)}
    if(method == "DirectFiltering"){res <- DirectFiltering(detrTimser(timser, detr), fc, tn, appliedFilter, k)}
    if(method == "IntegrandInterpolationMethod"){res <- IntegrandInterpolationMethod(detrTimser(timser, detr), fc, tn, appliedFilter, k)}
    if(method == "NoFilter"){res <- detrTimser(timser, detr)}
    return(res)
  }
  
  #### Preparation ####
  Imput.list <- Filter(Negate(is.null), Imput.list)
  if(is.null(New.param.names) == F){New.param.names <- c(Plot.x, New.param.names)}
  if(is.null(names(Imput.list)) == F){Keep.order <- names(Imput.list)}
  
  Imput.list <- lapply(Imput.list, function(x) x[!duplicated(x),])
  Imput.list <- lapply(Imput.list, function(x) x[intersect(Reduce(intersect, lapply(Imput.list, names)), names(x))])
  if(is.null(Calib.name) == F){
    for(i in 1:length(Imput.list)){
      Imput.list[[i]] <- Imput.list[[i]][c(Plot.x, names(Imput.list[[i]])[grep(paste(Calib.name, "$", sep = ""), names(Imput.list[[i]]))])]
      names(Imput.list[[i]]) <- gsub(paste(Calib.name, "$", sep = ""), "", names(Imput.list[[i]]))
    }
  }
  if(is.null(Keep.clim) == F){Imput.list <- lapply(Imput.list, function(x) x[intersect(c(Plot.x, Keep.clim), names(x))])}
  if(is.null(Limits) == F){
    Imput.list <- lapply(Imput.list, function(x){
      x <- x[x[Plot.x] >= Limits[1] & x[Plot.x] <= Limits[length(Limits)],]
      if(nrow(x) == 0){x <- NULL}
      return(x)
    })
    Imput.list <- Filter(Negate(is.null), Imput.list)
  }
  
  if(Scaling == T){Imput.list <- lapply(Imput.list, scale_ensemble_df, a = Scale.min, b = Scale.max)} # New scaling (Par triplet de mean, max, min)
  
  if(Anomaly == T & Binning == F){Imput.list <- lapply(Imput.list, function(y){cbind(y[1], apply(y[2:ncol(y)], 2, function(x) x <- x - x[1]))})}
  
  if(length(Imput.list) == 0){
    print("**** Not of the input data.frame is within the settled limits ! ****")
    return(NULL)}
  
  #### Corit cleaning ####
  if(Interpolation == T){Interpolation.method <- "InterpolationMethod"}
  else{Interpolation.method <- "NoFilter"}
  
  Imput.list <- lapply(Imput.list, function(x){
    if(Interpolation == T){
      x2 <- seq(min(x$Age, na.rm = T), max(x$Age, na.rm = T), by = Interp.step)
      x2 <- data.frame(Age = x2)
    }
    else{x2 <- x}
    
    for(i in 2:ncol(x)){
      if(all(is.na(x[i])) == F){
        P1 <- zoo::zoo(as.vector(unlist(x[i])), order.by = as.vector(unlist(x[1])))
        Cor <- Corit.simple(
          timser = P1,
          detr = Detrending,                       #remove linear trend time series
          method = Interpolation.method,
          appliedFilter = "gauss",
          fc = 1/Filter.cut.off, #cut-off frequency
          dt = Interp.step,    #time step for the interpolation
          int.method = "linear",  #kind of interpolation
        )
        
        x2[1] <- zoo::fortify.zoo(Cor)[1]
        x2[i] <- zoo::fortify.zoo(Cor)[2]
      }
      else{
        x2[1] <- x2[1]
        x2[i] <- NA
        
      }
    }
    names(x2) <- names(x)
    return(x2)
  }
  )
  
  #### Stacking and binning ####
  if(Stacking == T){
    Imput.list <- do.call(rbind, Imput.list)
    if(nrow(Imput.list)>1){
      if(Binning == T){
        if(is.null(Limits) == T){
          print(paste("**** No limits provided. We binned the data from the max(", Plot.x, ") to the min(", Plot.x, "). Use 'Limits' argument for manual limits. ****", sep = ""))
          Limits <- c(min(Imput.list, na.rm = T), max(Imput.list, na.rm = T))
        }
        if(is.null(Windows.length) == T){
          Windows.length <- round((Limits[2] - Limits[1])/15, digits = -1)
          print(paste("**** No 'Windows.length' provide, we binned the data by ", Windows.length, ". Use 'Windows.length' argument else. ****", sep = ""))
        }
        TW.x <- seq(Limits[1], Limits[2], Windows.length)
        Imput.list <- Bin.by.TW(Imput.list, Anomaly = Anomaly, Add.bin.count = Add.bin.count)
        
        if(Bin.sd == T & is.null(New.param.names) == F){New.param.names <- c(New.param.names, paste("SD", New.param.names[-c(1)], sep = "_"))}
        if(Add.bin.count == T & is.null(New.param.names) == F){New.param.names <- c(New.param.names, "N")}
      }
      
      #### Add type (or zone) name and change param.clim names ####
      if(is.null(Type.name) == F){
        Imput.list$Zone <- Type.name
        New.param.names <- c(New.param.names, "Zone")
      }
      
      if(is.null(New.param.names) == F){names(Imput.list) <- New.param.names}
      
      
      #### If not binned, return as list ####
      if(Return.list == T){
        if(Binning == F){
          Site.names <- gsub("\\..*", "", row.names(Imput.list))
          row.names(Imput.list) <- paste0(Site.names, ".Ech_", ave(seq_along(Site.names), Site.names, FUN = seq_along))
          Imput.list <- split(Imput.list, Site.names)
          Imput.list <- Imput.list[Keep.order]
          Imput.list <- Filter(Negate(is.null), Imput.list)
          
        }
        else{print("**** Be carefull! 'Return.list = T' is only for 'Binning = F'. ****")}
      }
      return(Imput.list)
    }
    else{
      print(paste("**** Be carefull! The limits do not contain any data in", Type.name, "****", sep = " "))
      return(NULL)
    }
  }
  
  #### No stacking and binning ####
  if(Stacking == F){
    if(Binning == T){
      for(i in 1:length(Imput.list)){
        Imput.list.i <- Imput.list[[i]]
        Limits.i <- Limits
        Windows.length.i <- Windows.length
        
        if(is.null(Limits.i) == T){
          print(paste("**** No limits provided. We binned the data from the max(", Plot.x, ") to the min(", Plot.x, "). Use 'Limits' argument for manual limits. ****", sep = ""))
          Limits.i <- c(min(Imput.list.i, na.rm = T), max(Imput.list.i, na.rm = T))
        }
        if(is.null(Windows.length.i) == T){
          Windows.length.i <- round((Limits.i[2] - Limits.i[1])/15, digits = -1)
          print(paste("**** No 'Windows.length' provide, we binned the data by ", Windows.length, ". Use 'Windows.length' argument else. ****", sep = ""))
        }
        TW.x <- seq(Limits.i[1], Limits.i[2], Windows.length.i)
        
        Imput.list.i <- Bin.by.TW(Imput.list.i, Anomaly = Anomaly, Add.bin.count = Add.bin.count)
        if(Bin.sd == T & is.null(New.param.names) == F){New.param.names <- c(New.param.names, paste("SD", New.param.names[-c(1)], sep = "_"))}
        if(Add.bin.count == T & is.null(New.param.names) == F){New.param.names <- c(New.param.names, "N")}
        if(is.null(New.param.names) == F){names(Imput.list.i) <- New.param.names}
        
        Imput.list[[i]] <- Imput.list.i
      }
      return(Imput.list)
    }
    else{
      if(Return.list == T){
        Imput.list <- lapply(Imput.list, function(df){
          Site.names <- gsub("\\..*", "", row.names(df))
          row.names(df) <- paste0(Site.names, ".Ech_", ave(seq_along(Site.names), Site.names, FUN = seq_along))
          return(df)
        })
        Imput.list <- Imput.list[Keep.order]
        Imput.list <- Filter(Negate(is.null), Imput.list)
        
        if(is.null(Type.name) == F){
          Imput.list <- lapply(Imput.list, function(df){
            df$Zone <- Type.name
            return(df)})
          New.param.names <- c(New.param.names, "Zone")
        }
        
        if(is.null(New.param.names) == F){
          Imput.list <- lapply(Imput.list, function(df){
            names(df) <- New.param.names
            return(df)})
        }
        
        return(Imput.list)
      }
      else{
        Imput.list <- do.call(rbind, Imput.list)
        
        if(is.null(Type.name) == F){
          Imput.list$Zone <- Type.name
          New.param.names <- c(New.param.names, "Zone")}
        
        if(is.null(New.param.names) == F){names(Imput.list) <- New.param.names}
        
        return(Imput.list)
      }
    }
  }
  
}

Aggregate.metadata <- function(df, group_col, site_col, spatial_cols = NULL, mean_cols = NULL,
                               first_cols = NULL, range_cols = NULL, concat_cols = NULL, collapse_sep = "|") {
  
  #### Warnings ####
  if(!group_col %in% names(df)){print("**** Be carefull ! Error in 'group_col'. ****"); stop()}
  if(!site_col  %in% names(df)){print("**** Be carefull ! Error in 'site_col'. ****"); stop()}
  
  mean_cols <- mean_cols[which(mean_cols %in% names(df))]
  spatial_cols <- spatial_cols[which(spatial_cols %in% names(df))]
  concat_cols <- concat_cols[which(concat_cols %in% names(df))]
  first_cols <- first_cols[which(first_cols %in% names(df))]
  
  Input.cols <- c(group_col, site_col, spatial_cols, mean_cols, first_cols, names(range_cols), concat_cols)
  Keep.col.order <- names(df)
  
  if(length(setdiff(names(df), Input.cols)) > 0){
    print("**** You forgot these columns in the aggregation, please provide them. ****")
    print(setdiff(names(df), Input.cols))}
  
  #### Aggregation ####
  df <- df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      
      #### Spatial barycenter ####
      dplyr::across(all_of(spatial_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
      
      #### Mean aggregation ####
      dplyr::across(all_of(mean_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
      
      #### First-row metadata ####
      dplyr::across(all_of(first_cols), dplyr::first, .names = "{.col}"),
      
      #### Range aggregation (min/max) ####
      dplyr::across(all_of(names(range_cols)), ~ do.call(range_cols[[cur_column()]], list(.x, na.rm = TRUE)),.names = "{.col}"),
      dplyr::across(all_of(concat_cols), ~ paste(unique(na.omit(.x)), collapse = collapse_sep), .names = "{.col}"),
      
      Site_full = paste0(dplyr::first(.data[[group_col]]), " (", paste(unique(.data[[site_col]]), collapse = ", "), ")"),
      Site = first(.data[[site_col]]), .groups = "drop")
  
  #### Clean export ####
  df <- df[match(Keep.col.order, names(df))[!is.na(match(Keep.col.order, names(df)))]]
  df <- data.frame(df)
  row.names(df) <- gsub("\\|.*", "", df$ID)
  return(df)
}

DES.weights <- function(models_cal, Xpred,
                        K_cal = 100, K_pred = 10,
                        causal = F, n_cores = 2,
                        pseudocount = 1e-6,
                        metric = "aitchison",
                        Cut.off = NULL, Soft.max = 1,
                        Local.subsampling = T,
                        parallelize = TRUE) {
  #### Setting values ####
  Xpred <- as.data.frame(Xpred)
  n_models <- length(models_cal)
  n_pred   <- nrow(Xpred)
  
  #### Distance method ####
  compute_distance <- function(Xcal_loc, Xpred_loc, metric, pseudocount=1e-6) {
    #### Aitchison distance ####
    if(metric == "aitchison"){
      clr_Xcal  <- compositions::clr(as.matrix(Xcal_loc) + pseudocount)
      clr_Xpred <- matrix(colMeans(compositions::clr(as.matrix(Xpred_loc) + pseudocount)),
                          nrow = nrow(clr_Xcal), ncol = ncol(clr_Xcal), byrow = TRUE)
      d <- sqrt(rowSums((clr_Xcal - clr_Xpred)^2))}
    
    #### Hellinger distance ####
    if(metric == "hellinger") {
      Xcal_sqrt  <- sqrt(as.matrix(Xcal_loc))
      Xpred_sqrt <- matrix(sqrt(colMeans(as.matrix(Xpred_loc))), nrow=nrow(Xcal_sqrt), ncol=ncol(Xcal_sqrt), byrow=TRUE)
      d <- sqrt(rowSums((Xcal_sqrt - Xpred_sqrt)^2))} 
    
    #### Bray-Curtis distance ####
    if(metric == "bray") {
      bray_curtis_distance <- function(Xcal_loc, Xpred_loc) {
        xpred <- colMeans(Xpred_loc)
        xcal <- as.matrix(Xcal_loc)
        d <- rowSums(abs(xcal - matrix(xpred, nrow=nrow(xcal), ncol=ncol(xcal), byrow=TRUE))) /
          rowSums(xcal + matrix(xpred, nrow=nrow(xcal), ncol=ncol(xcal), byrow=TRUE))
        return(d)
      }
      d <- bray_curtis_distance(Xcal_loc, Xpred_loc)}
    
    #### Jensen-Shannon distance ####
    if(metric == "js") {
      Xcal_mat <- as.matrix(Xcal_loc)
      Xpred_vec <- colMeans(as.matrix(Xpred_loc))
      js_dist <- function(p,q) {
        m <- 0.5*(p+q)
        sqrt(0.5*(sum(p*log(p/m)) + sum(q*log(q/m))))
      }
      d <- apply(Xcal_mat, 1, function(row) js_dist(row + pseudocount, Xpred_vec + pseudocount))
    } 
    
    if(!metric %in% c("js", "aitchison", "bray", "hellinger")) {
      print("**** Unknown metric selected, please choise between 'js', 'aitchison', 'bray' or 'hellinger'. ****")
      stop()}
    return(d)
  }
  
  #### Add missing species ####
  all_species <- unique(c(colnames(Xpred), unlist(lapply(models_cal, colnames))))
  
  # Rebuild Xpred with all species
  for (sp in setdiff(all_species, colnames(Xpred))) {
    Xpred[[sp]] <- 0
  }
  Xpred <- Xpred[, all_species, drop = FALSE]  # reorder consistently
  
  # Rebuild each calibration dataset
  models_cal <- lapply(models_cal, function(Xcal) {
    Xcal <- as.data.frame(Xcal)
    for (sp in setdiff(all_species, colnames(Xcal))) {
      Xcal[[sp]] <- 0
    }
    Xcal <- Xcal[, all_species, drop = FALSE]
    return(Xcal)
  })
  
  #### Function for each past sample ####
  compute_single <- function(i) {
    half_K <- floor(K_pred / 2)
    if(causal){
      start_idx <- max(1, i - K_pred)
      end_idx   <- i
    } else {
      start_idx <- max(1, i - half_K)
      end_idx   <- min(n_pred, i + half_K)
    }
    
    Xpred_loc <- Xpred[start_idx:end_idx, , drop = FALSE]
    sim <- numeric(n_models)
    
    for (m in seq_len(n_models)) {
      Xcal <- models_cal[[m]]
      Kc <- min(K_cal, nrow(Xcal))
      if (Kc == 0) {sim[m] <- 0; next}
      
      #### KDE (kernel density weighting) ####
      if(Local.subsampling == T){
        nnc <- FNN::get.knnx(data = Xcal, query = Xpred[i, , drop = FALSE], k = Kc)$nn.index
        Xcal_loc <- Xcal[nnc, , drop = FALSE]
      }
      else{Xcal_loc <- Xcal}
      
      # Handle single-row Xpred_loc
      if (nrow(Xpred_loc) == 1) {Xpred_loc <- rbind(Xpred_loc, Xpred_loc)}
      
      #### Distance calculation ####
      d <- compute_distance(Xcal_loc, Xpred_loc, metric=metric, pseudocount=pseudocount)
      sim[m] <- 1 / (1 + mean(d))
    }
    
    #### Weighting and Cut-off ####
    sim <- sim^Soft.max
    sim <- sim / sum(sim)
    
    if(is.null(Cut.off) == F){
      if (!is.null(Cut.off)) sim[sim < Cut.off] <- 0
      if (sum(sim) == 0) {sim <- rep(0, n_models)}
      else {sim <- sim / sum(sim)}
    }
    
    return(sim)
  }
  
  #### Apply to each past sample (parallelize or not) ####
  if (parallelize && n_cores > 1) {
    W_list <- parallel::mclapply(1:n_pred, compute_single, mc.cores = n_cores)
  } else {
    W_list <- lapply(1:n_pred, compute_single)
  }
  
  W <- data.frame(do.call(rbind, W_list))
  if (!is.null(names(models_cal))) {
    names(W) <- names(models_cal)
  } else {
    names(W) <- paste0("Model", seq_len(n_models))
  }
  return(W)
}

Pollen.Ensemble.RMSE <- function(MPal, Plot.x, Keep.clim, Msurf, Mpal.raw = NULL, Mcalib = NULL, Ensemble.WM = T, DES.Ensemble = F, GDGT.plot.merge = F, Variance.method = "quadratic", Temp.col.alpha = 0.1, Cores.lab = NULL,
                                 Limites = NULL, Select.interv = 1000, Lab.age.size = 6, Zone.clim = NULL, Name.zone = NULL, Temp.zone = NULL, Clim.lab = NULL, Dot.size = 6,
                                 Manual.vlines = NULL, Smooth.sd = T, Smooth.T = F, Smooth.param = 0.6, Surf.val = NULL, Name.core = NULL, Surf.val.max.age = NULL, Title.display = T,
                                 Manual.y.val = NULL, Axis.right = F, Anomaly = F, Remove.calib = NULL, FT.plot.merge = T,
                                 DES.cut.off = NULL, Show.DES.contrib = F, Display.sim = T, DES.Soft.max = 1, DES.distance.metric = "aitchison", DES.KDE = T, DES.parallelize = F,
                                 Add.lim.space = F, Save.plot = NULL, W = NULL, H = NULL, Save.path = NULL, Show.box = F, Plot.Rds = NULL, return.plot = F){
  #### Check if same number of samples in each paleo reconstructions ####
  if(length(unique(sapply(MPal, nrow))) > 1){
    print(paste("Be carefull! The ages are not similar from the differentes reconstruction for", Name.core))
    all_ages <- sort(unique(unlist(lapply(MPal, function(df) df$Age))))
    MPal <- lapply(MPal, function(df) {
      template <- data.frame(Age = all_ages)
      out <- full_join(template, df, by = "Age")
      return(out)
    })
  }
  
  #### Clean data if pb ####
  if(typeof(MPal[[1]]$Age) != "double"){MPal <- lapply(MPal, function(x){x$Age  <- as.double(x$Age); return(x)})}
  
  #### Weighted variance function ####
  combine_similarity_rmse_models <- function(S, rmses, epsilon = 1e-12) {
    n_pred <- nrow(S)
    datasets <- colnames(S)
    models <- names(rmses)
    
    # Extract dataset name from model names (assumes format "MODEL.DATASET")
    model_datasets <- sub(".*\\.", "", models)
    
    W <- matrix(0, n_pred, length(models))
    
    for (i in 1:n_pred) {
      w_raw <- numeric(length(models))
      
      for (m in seq_along(models)) {
        ds <- model_datasets[m]
        if (!(ds %in% datasets)) next  # skip if dataset not in similarity matrix
        sim_val <- S[i, ds]
        w_raw[m] <- sim_val / (rmses[m]^2 + epsilon)
      }
      
      if (sum(w_raw) == 0) {
        W[i, ] <- rep(1 / length(models), length(models))  # fallback uniform
      } else {
        W[i, ] <- w_raw / sum(w_raw)
      }
    }
    
    colnames(W) <- models
    return(W)
  }
  
  row_weighted_variance <- function(mat, weights, method = c("standard", "unbiased", "quadratic")) {
    method <- match.arg(method)
    
    if (ncol(mat) != length(weights)) {
      stop("Number of columns in matrix must match length of weights.")
    }
    
    apply(mat, 1, function(x) {
      w_mean <- sum(weights * x) / sum(weights)
      
      if (method == "standard") {
        return(sum(weights * (x - w_mean)^2) / sum(weights))
        
      } else if (method == "unbiased") {
        w_sum <- sum(weights)
        w_squared_sum <- sum(weights^2)
        correction = w_sum - (w_squared_sum / w_sum)
        return(sum(weights * (x - w_mean)^2) / correction)
        
      } else if (method == "quadratic") {
        return(sum(weights^2 * (x - w_mean)^2))
      }
    })
  }
  
  row_weighted_variance.DES <- function(mat, weights, method = c("standard", "unbiased", "quadratic")) {
    method <- match.arg(method)
    
    if (!all(dim(mat) == dim(weights))) {
      stop("Matrix of predictions and weights must have the same dimensions.")
    }
    
    sapply(seq_len(nrow(mat)), function(i) {
      x <- mat[i, ]
      w <- weights[i, ]
      w_mean <- sum(w * x) / sum(w)
      
      if (method == "standard") {
        return(sum(w * (x - w_mean)^2) / sum(w))
        
      } else if (method == "unbiased") {
        w_sum <- sum(w)
        w_squared_sum <- sum(w^2)
        correction <- w_sum - (w_squared_sum / w_sum)
        return(sum(w * (x - w_mean)^2) / correction)
        
      } else if (method == "quadratic") {
        return(sum(w^2 * (x - w_mean)^2))
      }
    })
  }
  
  Keep.plot.x <- MPal[[1]][Plot.x]
  
  #### Graphical option ####
  Col.pol <- "#BEA33A"
  if(is.null(Limites) == T){Limites = c(min(MPal[[1]][Plot.x], na.rm = T), max(MPal[[1]][Plot.x], na.rm = T))}
  
  if(Show.box == T){
    My_border <- element_rect(fill = NA, colour = "grey30", linewidth = .3)
    Axis.x <- element_blank(); Axis.y <- element_blank()
  }
  else{
    My_border <- element_blank()
    Axis.x <- element_line(); Axis.y <- element_line()
  }
  
  if(GDGT.plot.merge == T){
    Title.age = NULL
    Axis.x <- element_blank(); Ticks.x <- element_blank(); Axis.lab.x <- element_blank()
    if(Add.lim.space == T){
      Limites[2] <- (Limites[2] + Limites[2]*0.23)}
    Legende.Core.name.chart <- c(17, rep(0, (length(Keep.clim)-1)))
  }
  else{
    Title.age = "Time (cal. year BP)"
    Ticks.x <- element_line(); Axis.lab.x <- element_text(angle = 45, hjust = 1, size = Lab.age.size, colour = "grey55")
    Legende.Core.name.chart <- c(16, rep(0, (length(Keep.clim)-1)))
  }
  
  if(Title.display == F){My_title <- element_blank()}
  else{My_title <- element_text(size = 17, hjust = 0.5)}
  
  Plot.list.col <- list()
  
  My_color_models <- c("WASTDB" = "#e2a064ff", "WAST" = "#e2a064ff", 
                       "NMSDB" = "#F3A481", "MDB" = "#963326", "TUDB" = "#0094AF", "TUSDB" = "#0094AF",
                       "ST" = "#91C4DD", "STDB" = "#91C4DD", 
                       "ACADB" = "#963326", "CAUCDB" = "#0F3361",
                       "COSTDB" = "#f3c768ff", "COST" = "#f3c768ff",
                       "MEDTEMP" = "#F3A481", "TEMPSCAND" = "#91C4DD", 
                       "EAPDB" = "#0F3361", "TAIGDB" = "#32156eff")
  
  #### Climat zone settings ####
  Rect.color.scale <- c("W" = "#E76D51",
                        "C" = "#75AADB",
                        "G" = "grey40",
                        "D" = "#c67f05",
                        "Wt" = "#004266")
  Rect.color.scale <- Rect.color.scale[unique(Temp.zone)]
  
  if(is.null(Zone.clim) ==F){
    yo = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                    xmax = Zone.clim[seq(2,length(Zone.clim), by=2)], 
                    Temp.col = Temp.zone)
    Rect.clim <- geom_rect(data = yo, inherit.aes = F,
                           mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf, fill = Temp.col),
                           alpha = Temp.col.alpha, color = "grey", size = 0.3, linetype = 2, na.rm = T)} 
  else{Rect.clim <- NULL}
  if(is.null(Name.zone) == F){yo2 = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                                               xmax = Zone.clim[seq(2,length(Zone.clim), by=2)],
                                               Temp.col = Temp.zone,
                                               Title.zone = Name.zone,
                                               Model = Model[length(Model)])}
  else{yo2 = data.frame(xmin = 0, xmax = 0, Temp.col = "", Title.zone = "")}
  
  #### Case only some Surface Values ####
  if(is.null(Surf.val) == F){
    A <- setNames(data.frame(matrix(ncol = length(setdiff(Keep.clim, names(Surf.val))), nrow = nrow(Surf.val))), setdiff(Keep.clim, names(Surf.val))) 
    row.names(A) <- row.names(Surf.val)
    A[is.na(A)] <- NA
    A <- cbind(Surf.val, A)
    Surf.val <- A[,sort(names(A))]}
  
  #### Smooth curve ####
  if(Smooth.T == T){
    My_smooth <- geom_smooth(method = "loess", se = Smooth.sd, fullrange = F, level = 0.95, linetype="solid", formula = "y ~ x",
                             size = .8, show.legend = F, fill = Col.pol, color = Col.pol, span = Smooth.param, alpha = 0.2)}
  else{My_smooth <- NULL}
  #### Dynamic Ensemble Selection based on Aitchison distance ####
  if(DES.Ensemble == T){
    if(is.null(Mcalib) == T | is.null(Mpal.raw) == T){
      print("**** Dynamic Ensemble Selection impossible without 'Mcalib' and 'Mpal.raw'. ****")
      Msim <- NULL; DES.Ensemble = F}
    else{
      Keep.calib <- unique(gsub(".*\\.", "", names(MPal)))
      Keep.calib2 <- unique(gsub(".*\\.", "", names(Msurf)))
      Keep.calib <- intersect(Keep.calib, Keep.calib2)
      
      if(length(names(Mcalib)[!names(Mcalib) %in% Keep.calib]) > 0){
        
        print(paste("**** These calib. were removed from the DES, since they were not applied on ", Name.core, ". ****", sep = ""))
        print(names(Mcalib)[!names(Mcalib) %in% Keep.calib])
      }
      
      
      if(all(Keep.calib %in% names(Mcalib)) == F){
        print(paste("**** Dynamic Ensemble Selection impossible without ", Keep.calib[!Keep.calib %in% names(Mcalib)], " given in 'Mcalib'. ****", sep = ""))
        Msim <- NULL
        DES.Ensemble = F
      }
      else{
        Mcalib <- Mcalib[names(Mcalib) %in% Keep.calib]
        
        if(unique(lapply(MPal, nrow)) != nrow(Mpal.raw)){
          print("**** The number of sample in 'MPal' and 'Mpal.raw are different ! ****")
          Msim <- NULL; Msim.m <- NULL; Show.DES.contrib = F; Display.sim = F; DES.Ensemble = F
        }
        else{
          Msim <- DES.weights(
            models_cal = Mcalib, Xpred = Mpal.raw, 
            K_cal = min(sapply(Mcalib, nrow)), K_pred = signif(nrow(Mpal.raw)/10, digits = 0),
            Cut.off = DES.cut.off, Soft.max = DES.Soft.max, causal = F, n_cores = 2, 
            metric = DES.distance.metric, Local.subsampling = DES.KDE, parallelize = DES.parallelize, )
          
          if(Show.DES.contrib == T & Display.sim == T){
            Msim.m <- cbind(Keep.plot.x, Msim)
            Msim.m <- melt(Msim.m, id = Plot.x)
            names(Msim.m)[1] <- "Plot.x"
          }
        }
      }
    }
  }
  else{Msim <- NULL; Msim.m <- NULL}
  
  #### Loop on climate param ####
  Res <- list()
  Res.contrib <- list()
  for(j in Keep.clim){
    List.predict <- list(); RMSE <- c(NA)
    for(i in 1:length(Msurf)){
      Name.model.full <- paste(j, names(Msurf)[i], sep = ".")
      
      Method.i <- names(Msurf)[i]
      if(is.null(Remove.calib) == F){
        print("**** These calibration datasets were removed. ****")
        print(any(Remove.calib %in% Name.model.full))
        if(any(Remove.calib %in% Name.model.full)){next}
      }
      
      if(j %in% names(MPal[[Method.i]]) == F){next}
      
      List.predict[[i]] <- setNames(MPal[[Method.i]][j], Method.i)
      RMSE.i <- Msurf[[Method.i]]$Best.Param
      RMSE.i <- RMSE.i[row.names(RMSE.i) == j, "RMSE.choice"]
      RMSE[i] <- RMSE.i
    }
    List.predict <- Filter(Negate(is.null), List.predict)
    RMSE <- Filter(Negate(is.na), RMSE)
    List.predict <- do.call(cbind, List.predict)
    
    #### Weighted Mean RMSE ####
    if(Ensemble.WM == T){
      #### DES ####
      if(DES.Ensemble == T){
        names(RMSE) <- names(List.predict)
        wg <- combine_similarity_rmse_models(Msim, RMSE)
        
        Ensemble_mean <- rowSums(List.predict * wg)
        prediction_matrix <- as.matrix(List.predict)
        weighted_variance <- row_weighted_variance.DES(prediction_matrix, wg, method = Variance.method)
        
        Res[[paste(j, "FT_Ensemble", sep = "_")]]    <- Ensemble_mean
        Res[[paste(j, "FT_Ensemble_MaxI", sep = "_")]] <- Ensemble_mean + 1.96 * sqrt(weighted_variance)
        Res[[paste(j, "FT_Ensemble_MinI", sep = "_")]] <- Ensemble_mean - 1.96 * sqrt(weighted_variance)
        
        if(Show.DES.contrib == T & Display.sim == F){
          colnames(wg) <- paste(j, colnames(wg), sep = "_")
          Res.contrib[[j]] <- wg
        }
      }
      
      #### Only RMSE ####
      if(DES.Ensemble == F){
        wg <- 1/RMSE
        wg <- wg/sum(wg, na.rm = T)
        Res[[paste(j, "FT_Ensemble", sep = "_")]] <- as.numeric(as.matrix(List.predict) %*% wg)
        prediction_matrix <- as.matrix(List.predict)
        weighted_variance <-row_weighted_variance(prediction_matrix, wg, method = Variance.method)
        Res[[paste(j, "FT_Ensemble_MaxI", sep = "_")]] <- Res[[paste(j, "FT_Ensemble", sep = "_")]] + 1.96*sqrt(weighted_variance)
        Res[[paste(j, "FT_Ensemble_MinI", sep = "_")]] <- Res[[paste(j, "FT_Ensemble", sep = "_")]] - 1.96*sqrt(weighted_variance)
      }
    }
  }
  Res <- cbind(Keep.plot.x, do.call(cbind, Res))
  
  #### Prep melted data ####
  
  melt <- melt(Res, names(Res)[c(1, grep("MinI", names(Res)), grep("MaxI", names(Res)))])
  names(melt)[c(ncol(melt)-1, ncol(melt))] <- c("Param.clim", "Mean_val")
  
  melt <- melt(melt, c("Param.clim", "Mean_val", names(melt)[c(1, grep("MaxI", names(melt)))]))
  names(melt)[c(ncol(melt)-1, ncol(melt))] <- c("MinI", "Min_val")
  
  melt <- melt(melt, c("Param.clim", "Mean_val", "MinI", "Min_val", Plot.x))
  names(melt)[c(ncol(melt)-1, ncol(melt))] <- c("MaxI", "Max_val")
  
  names(melt)[names(melt) == Plot.x] <- "Plot.x"
  melt$Param.clim <- gsub("_FT_Ensemble", "", melt$Param.clim)
  melt$MaxI <- gsub("_FT_Ensemble_MaxI", "", melt$MaxI)
  melt$MinI <- gsub("_FT_Ensemble_MinI", "", melt$MinI)
  
  melt <- melt[melt$Param.clim == melt$MinI & melt$Param.clim == melt$MaxI,]
  melt$Param.clim <- factor(melt$Param.clim, levels = Keep.clim, ordered = T)
  melt <- melt[!duplicated(melt),]
  
  #### Melt contrib ####
  if(Show.DES.contrib == T & Display.sim == F){
    Res.contrib <- cbind(Keep.plot.x, do.call(cbind, Res.contrib))
    Res.contrib <- melt(Res.contrib, Plot.x)
  }
  
  #### Plot ####
  for(i in 1:length(Keep.clim)){
    melt.i <- melt[melt$Param.clim == Keep.clim[i],]
    if(is.null(Surf.val) == F){Surf.val.j <- subset(Surf.val, select = c(Keep.clim[i]))}
    
    #### Surf val ####
    if(is.null(Surf.val) == F){
      Surf.val.j.i <- Surf.val.j[Name.core,]
      if(is.na(Surf.val.j.i) == T){Surf.val.j.i = NULL}
    }
    
    if(is.null(Surf.val) == F){
      if(is.na(Surf.val.j.i) == F){
        if(is.null(Surf.val.max.age) == T){
          Surf.line <- geom_hline(yintercept = Surf.val.j.i, linetype = "dotdash", linewidth = 0.55, color = "black", alpha = 0.55)
          Surf.dot <- NULL}
        else{
          A = data.frame(x = Limites[1], xend = Surf.val.max.age, yend = Surf.val.j.i, y = Surf.val.j.i)
          Surf.line <- geom_segment(inherit.aes = F, data = A, aes(yend = yend, x = x, xend = xend, y = y), linetype = "dotdash", linewidth = 0.55, color = "black", alpha = 0.55)
          Surf.dot <- geom_point(inherit.aes = F, data = A, aes(x = x, y = y),color = "darkred", alpha = 1, size = Dot.size*1.5, shape = 9)
        }
        
      }
      if(is.na(Surf.val.j.i) == T | is.null(Surf.val.j.i) == T){Surf.line <- NULL; Surf.dot <- NULL}}
    else{Surf.line <- NULL; Surf.dot <- NULL}
    
    if(Show.DES.contrib == T & DES.Ensemble == T & is.null(Surf.dot) == F){
      Surf.dot.dummy <- geom_point(inherit.aes = F, data = A, aes(x = x, y = 0.5),color = "white", alpha = 1, size = Dot.size*1.5, shape = 9)
    }  
    
    #### Breaks ####
    if(is.null(Manual.y.val) == T){
      if(is.null(Surf.val) == F){
        Kmax <- max(melt.i$variable, Surf.val.j.i, na.rm = T)
        Kmin <- min(melt.i$variable, Surf.val.j.i, na.rm = T)}
      else{
        Kmax <- max(melt.i$variable, na.rm = T)
        Kmin <- min(melt.i$variable, na.rm = T)
      }
      
      if(Kmax >= 40){
        Kmin <- round(Kmin, digits = -1)
        Kmin = floor(Kmin / 50) * 50
        Kmax <- round(Kmax, digits = -1)
        Kmax = ceiling(Kmax / 50) * 50
        Kmid <- Kmin + (Kmax - Kmin)*1/3
        Kmid = ceiling(Kmid / 50) * 50
        Kmid2 <- Kmin + (Kmax - Kmin)*2/3
        Kmid2 = ceiling(Kmid2 / 50) * 50
      }
      if(Anomaly == T & is.null(Surf.val) == F){
        Kmid <- Kmin + (Kmax - Kmin)*1/3
        Kmid2 <- Kmin + (Kmax - Kmin)*2/3
      }
      else{
        Kmin <- round(Kmin, digits = 1)
        Kmin = floor(Kmin / 0.5) * 0.5
        Kmax <- round(Kmax, digits = 1)
        Kmax = ceiling(Kmax / 0.5) * 0.5
        Kmid <- Kmin + (Kmax - Kmin)*1/3
        Kmid = ceiling(Kmid / 0.5) * 0.5
        Kmid2 <- Kmin + (Kmax - Kmin)*2/3
        Kmid2 = ceiling(Kmid2 / 0.5) * 0.5
      }
      Range.round <- c(Kmin, Kmid, Kmid2, Kmax)
      if(Kmax >= 40){
        Range.round <- round(Range.round, digits = -1)
      }
      if(Anomaly == T & is.null(Surf.val) == F){
        Range.round <- round(Range.round, digits = 2)
      }
      else{
        Range.round <- round(Range.round, digits = 1)
      }
      
    }
    else{
      Clim.break <- gsub("\\..*", "", names(Manual.y.val))
      Core.break <- gsub(".*\\.", "", names(Manual.y.val))
      Match.clim <- which(Clim.break %in% Keep.clim[i])
      Match.tot <- unique(Match.clim)
      Range.round <- Manual.y.val[[Match.tot]]
    }
    
    My_axis_clim <- data.frame(A = Range.round[c(1,length(Range.round))], B = Limites)
    if(Axis.right == F){Axis.side <- "l"; Axis.side2 <- "left"}
    else{Axis.side <- "r"; Axis.side2 <- "right"}
    
    if(is.null(Cores.lab) == F){
      if(is.null(Manual.y.val) == T){Scale_y <- scale_y_continuous(name = Cores.lab,  breaks = scales::breaks_extended(n = 4))}
      else{Scale_y <- scale_y_continuous(name = Cores.lab, breaks = Range.round, position = Axis.side2)}
    }
    else{Scale_y <- NULL}
    
    #### Plot ####
    p <- ggplot(melt.i, aes(x = Plot.x, y = Mean_val))+
      #### Rectangles ####
    scale_fill_manual(values = Rect.color.scale, guide = "none")+
      Rect.clim +
      geom_vline(xintercept = Manual.vlines, col = "grey30", lty = 2, alpha = 0.7)+
      new_scale_fill()+
      
      #### Items ####
    facet_wrap(vars(Param.clim), scale = "free")+
      geom_ribbon(aes(ymin = Min_val, ymax = Max_val), fill = Col.pol, color = NA, size = 0.15, linetype = "dashed", alpha = 0.1)+
      
      xlab(Plot.x)+
      geom_point(size = Dot.size, color = Col.pol, shape = 1)+
      geom_line(color = Col.pol, linetype = "solid", linewidth = 0.5, alpha = 0.6)+
      scale_x_continuous(name = Title.age, breaks = round(seq(0, Limites[2], by = Select.interv)), expand = c(0.02, 0.02))+
      Scale_y + 
      My_smooth +
      ggtitle(Clim.lab[[i]])+
      guides(colour = guide_legend(nrow = 1))+        # force les legendes a salligner sur une unique ligne
      #### Annotations ####
    new_scale_color()+
      scale_color_manual(values = Rect.color.scale, guide = "none", name = NULL, labels = NULL, breaks = NULL, na.translate = FALSE)+
      ggthemes::geom_rangeframe(data = My_axis_clim, inherit.aes = F, mapping = aes(x = B, y = A), colour = "grey50", sides = Axis.side, size = 0.5) +
      Surf.line + Surf.dot +
      
      #### Theme ####
    theme(plot.background = element_blank(), axis.ticks.x = Ticks.x,
          plot.title = My_title, panel.grid = element_blank(),
          axis.title.y = element_text(size = Legende.Core.name.chart[i], colour = "grey20"),              # taile police du titre de l'axe
          strip.text = element_blank(),
          axis.line.x = Axis.x, axis.line.y = Axis.y,
          panel.background = My_border, strip.background = element_blank(),
          axis.text.x = Axis.lab.x, plot.margin = unit(c(0,0.2,0,0),"cm"),
          strip.text.y = element_text(size = 11.5, angle = 180),
          strip.placement = "outside",
          legend.position = "none")
    
    
    #### Plot contribution DES ####
    if(Show.DES.contrib == T & DES.Ensemble == T){
      p2 <- ggplot(Msim.m, aes(x = Plot.x, y = value, color = variable, fill = variable)) + Surf.dot.dummy +
        geom_vline(xintercept = Manual.vlines, col = "grey30", lty = 2, alpha = 0.7)+
        geom_col(position = "stack") +
        scale_color_manual(values = My_color_models)+
        scale_fill_manual(values = My_color_models)+ #ylab("Weights (%)")+
        scale_x_continuous(name = Title.age, breaks = round(seq(0, Limites[2], by = Select.interv)), limits = Limites, expand = c(0.02, 0.02))+
        scale_y_continuous(name = "Weights (%)", breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1), labels = c(0, 25, 50, 75, 100))+
        theme(plot.background = element_blank(),
              plot.title = element_blank(), panel.grid = element_blank(),
              strip.text = element_blank(),
              axis.title.x = element_blank(), strip.placement = "outside",
              axis.line.y.left = element_line(),
              axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
              panel.background = My_border, strip.background = element_blank(), plot.margin = unit(c(0,0.2,0,0),"cm"),
              legend.position = "top", legend.direction = "horizontal", legend.title = element_blank(),
        )
      
      p <- p2 / p + plot_layout(heights = c(25,75))
    }
    
    Plot.list.col[[i]] <- p
  }
  #### Facet plot clim / cores ####
  Formula.patch <- paste(paste("Plot.list.col[[", seq(length(Keep.clim)), "]]", sep = ""), collapse =  " | ")
  p <- eval(parse(text = Formula.patch))
  
  #### Export plots #### 
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(p, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", limitsize = F)}
    else{ggsave(Save.plot)}}
  
  if(is.null(Plot.Rds) == F){saveRDS(p, Plot.Rds)}
  if(is.null(Save.path) == F){saveRDS(Res, Save.path)}
  if(return.plot == T){return(p)}
  else{return(Res)}
}

Mat.corel.CWT.clim <- function(Mclim, Mtrait, I.confiance = NULL, Display.pval, Disp.R = "pie", Display = "full", ggplot.version = F,
                               Title, Label, Save.path = NULL, return.pick, Bar.pos = "b", my_color, Use.cor = "pairwise.complete.obs", Xlab.pos = NULL,
                               Good.tiles = NULL,  Bad.tiles = NULL, Xlab.rot = 45, Return.slope = F, Print.result = T, No.axis.lab = F,
                               Permutation.test = F, Nb.permutations = 1000, Save.pval = NULL, Bold = NULL, Order.hclust = F,
                               Method.cor, Save.plot = NULL, W, H, Order.by.alphabet, Average){
  #### Settings ####
  Save.tab = T
  if(missing(return.pick)){return.pick = F}
  if(missing(Order.by.alphabet)){Order.by.alphabet = F}
  if(missing(Label)){Label = T}
  if(missing(Average)){Average = T}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Method.cor)){Method.cor = "pearson"}
  if(missing(Display.pval)){Display.pval = "pch"}
  if(missing(Mtrait)){Mtrait = NULL}
  if(missing(Title)){Title = NULL}
  if(missing(my_color)){my_color = colorRampPalette(c("royalblue", "white", "darkorange"))(50)}
  if(Display.pval == "blank"){Subt = paste("*The cells are blanked when ever the p-value is over the", I.confiance, "significant level.")}
  if(Display.pval == "pch"){Subt = paste("*The cells are crossed when ever the p-value is over the", I.confiance, "significant level.")}
  if(Label == F){Subt <- NULL}
  
  #### Calculations ####
  
  if(identical(row.names(Mclim), row.names(Mtrait)) == F){print("****Be carefull !!! The 2 matrixes may be not sorted or similar shape !****")}
  if(is.null(Mtrait) == F){names(Mtrait) <- gsub("TRY_", "", names(Mtrait))}
  
  if(Order.by.alphabet == T){Mtrait <- Mtrait[order(names(Mtrait))]}
  
  Mcor <- cor(Mclim, Mtrait, use = Use.cor, method = Method.cor)
  Mvar <- var(Mclim, Mtrait)
  Mtot <- cbind(Mclim, Mtrait)
  
  if(is.null(I.confiance) == F){
    Resid <- cor.mtest(Mtot, conf.level = I.confiance)
    PV <- Resid$p
  }
  else{PV <- NULL}
  
  if(is.null(PV) == F){
    if(is.null(Mtrait) == F){PV <- PV[(ncol(Mtrait)+1):(ncol(Mtrait)+ncol(Mclim)), 1:ncol(Mtrait)]} 
    row.names(PV) <- row.names(Mcor)
    colnames(PV) <- colnames(Mcor)
  }
  
  if(Return.slope == T){
    Mtrait2 <- Mtrait
    Mclim2 <- Mclim
    slope_mat <- sapply(1:ncol(Mtrait2), function(j) {apply(Mclim2, 2, function(y) coef(lm(y ~ Mtrait2[, j]))[2])})
    colnames(slope_mat) <- colnames(Mcor)
  }
  
  #### Permutations test ####
  if(Permutation.test == T){
    print(paste("p-value calculated with permutation test made on ", Nb.permutations, " permutations.", sep = ""))
    PV.perm <- setNames(data.frame(matrix(NA, nrow = length(Mclim), length(Mtrait))), names(Mtrait))
    row.names(PV.perm) <- names(Mclim) 
    for(i in 1:length(names(Mclim))){
      for(j in 1:length(names(Mtrait))){
        set.seed(0)
        M.perm <- na.omit(cbind(Mclim[i], Mtrait[j]))
        PV.perm[i, j] <- perm.cor.test(M.perm[[1]], M.perm[[2]], num.sim = Nb.permutations)$p.value
      }
    }
    PV <- as.matrix(PV.perm)
  }
  
  #### Save p-values ####
  if(is.null(Save.pval) == F & is.null(PV) == F){
    Path.to.create <- gsub("(.*/).*\\.csv.*","\\1", Save.pval)
    dir.create(file.path(Path.to.create), showWarnings = F)
    write.table(PV, file = Save.pval, row.names=T, col.names = NA, sep=",", dec = ".")
    Save.path.RDS <- gsub("\\.csv", "\\.Rds", Save.pval)
    saveRDS(PV, Save.path.RDS)
  }
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Graphic param ####
  if(is.null(Bold) == F & ggplot.version == T){
    names(Bold) <- c("Var1", "Var2")
    Bold.case <- geom_tile(data = Bold, width = .95, height = .95, linewidth = 1, color = "grey10", fill = NA, alpha = .5)
  }
  else{Bold.case <- NULL}
  
  if(typeof(Bar.pos) == "double" & ggplot.version == F){
    Bar.pos <- "b"
    print("**** Legend position impossible, please choose between t, b, l or r. Default to b.")
  }
  
  #### Highlight some tiles ####
  if(is.null(Good.tiles) == F){
    heat_data <- reshape2::melt(Mcor)
    heat_data$highlight <- apply(heat_data[, c("Var1", "Var2")], 1, function(row) {any(mapply(function(x) all(row == x), Good.tiles))})
    Bold.tiles <- geom_tile(data = subset(heat_data, highlight), mapping = aes(x = Var1, y = Var2), color = "green", size = 1.2, fill = NA)
  }
  else{Bold.tiles <- NULL}
  
  if(is.null(Bad.tiles) == F){
    heat_data <- melt(Mcor)
    heat_data$highlight <- apply(heat_data[, c("Trait1", "Trait2")], 1, function(row) {any(mapply(function(x) all(row == x), Bad.tiles))})
    Bold.tiles.2 <- geom_tile(data = subset(heat_data, highlight), color = "black", size = 1.2, fill = NA)
  }
  else{Bold.tiles.2 <- NULL}
  
  #### Plots (coorplor) ####
  if(ggplot.version == F){
    CP <- corrplot(Mcor, type = Display, col = my_color, tl.pos = Xlab.pos,
                   tl.col="black", tl.srt = Xlab.rot, tl.cex = .7, method = Disp.R, cl.align.text = "l", cl.pos = Bar.pos, 
                   p.mat = PV, sig.level = (1-I.confiance), insig = Display.pval, pch.cex = 2)
    
    title(main = Title, sub = Subt)
  }
  
  #### Plots (ggplot) ####
  if(ggplot.version == T){
    if(Disp.R == "pie"|Disp.R == "circle"){Disp.R.gg <- "circle"}
    if(Disp.R == "number"|Disp.R == "square"){Disp.R.gg <- "square"}
    
    if(Xlab.rot == 0){Xjust <- 0.5}
    if(Xlab.rot != 0){Xjust <- 1}
    
    if(typeof(Bar.pos) != "double"){
      Bar.pos.gg <- "top"
      if(Bar.pos == "b" | Bar.pos == "bottom"){Bar.pos.gg <- "bottom"}
      if(Bar.pos == "n" | Bar.pos == "none"){Bar.pos.gg <- "none"}
    }
    else{Bar.pos.gg <- Bar.pos}
    
    if(No.axis.lab == T){Xlab.opt <- element_blank(); Ylab.opt <- element_blank()}
    else{Xlab.opt <- element_text(hjust = Xjust); Ylab.opt <- element_text()}
    
    CP <- ggcorrplot(Mcor, hc.order = Order.hclust, type = Display, title = Title, colors = c("royalblue", "white", "darkorange"),
                     tl.col="black", tl.srt = Xlab.rot, tl.cex = 10,
                     method = Disp.R.gg,
                     p.mat = PV,
                     sig.level = (1-I.confiance), insig = Display.pval, pch.cex = 5,
                     lab = TRUE)+
      Bold.case+ Bold.tiles.2 + Bold.tiles +
      scale_x_discrete(expand = c(0,0))+
      scale_y_discrete(expand = c(0,0))+
      theme(legend.position = Bar.pos.gg, panel.border = element_rect(NA, "black"),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            axis.text.x = Xlab.opt,
            axis.text.y = Ylab.opt,
            plot.margin=unit(c(0,0,0,0),"cm"), plot.background = element_blank(), panel.background = element_rect(fill = NA, colour = "black"))
    
    print(CP)
  }
  
  #### Average values ####
  if(Average == T){
    type.correl <- gsub(".*ACA", "", Title)
    type.correl <- gsub("\\s.*", "", type.correl)
    
    MeanR2 <- round((Mcor)^2, digits = 2)
    MeanR2 <- cbind(MeanR2, MeansCol = rowMeans(MeanR2, na.rm = T))
    MeanR2 <- rbind(MeanR2, MeansRow = colMeans(MeanR2, na.rm = T))
    
    if(is.null(PV) == F & Print.result == T){print(PV)}
    
    R.mean.ACAV.t.chel <- round(mean(MeanR2["MeansRow", grepl("_chel", colnames(MeanR2))], na.rm = T), digits = 2)
    R.mean.ACAV.t.wc <- round(mean(MeanR2["MeansRow", grepl("_wc", colnames(MeanR2))]), digits = 2)
    R.mean.ACAV.gl.clim <- round(mean(MeanR2[grepl("_gf", rownames(MeanR2)), "MeansCol"]), digits = 2)
    R.mean.ACAV.clim <- round(mean(MeanR2[!grepl("_gf", rownames(MeanR2)), "MeansCol"]), digits = 2)
    
    if(Print.result == T){
      print(paste("R² mean ACA-", type.correl, " vs. Climat-CHELSA:", R.mean.ACAV.t.chel, sep = ""))
      print(paste("R² mean ACA-", type.correl, " vs. Climat-WORLDCLIM:", R.mean.ACAV.t.wc, sep = ""))
      print(paste("R² mean ACA-", type.correl, "-NORMAL vs. Climat :", R.mean.ACAV.clim, sep = ""))
      print(paste("R² mean ACA-", type.correl, "-GAPFILLING vs. Climat :", R.mean.ACAV.gl.clim, sep = ""))
    }
  }
  
  #### Export datas ####
  if(is.null(PV) == F){Mfull <- list(R2 = round(Mcor, digits = 2), Var = round(Mvar, digits = 2),p.val = round(PV, digits = 5))}
  else{Mfull <- list(R2 = round(Mcor, digits = 2), Var = round(Mvar, digits = 2))}
  
  if(Return.slope == T){Mfull <- append(Mfull, list("Slope" = slope_mat))}
  
  if(is.null(Save.path) == F){
    Path.to.create <- gsub("(.*/).*\\.csv.*","\\1", Save.path)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    write.table(Mfull, file = Save.path, row.names=T, col.names=NA, sep=",", dec = ".")}
  
  #### Return ####
  if(is.null(Save.plot) == F){dev.off()}
  if(return.pick == T){return(CP)}
  else{return(Mfull)}
}
