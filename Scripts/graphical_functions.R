Plot.clim.phase <- function(MFT = NULL, MGDGT = NULL, MXRF = NULL, Md18O = NULL, Site.order, Limites, Save.Rds = NULL, Save.plot = NULL, Save.plot.Rds = NULL, W = 500, H = 500, Save.top = NULL,
                            Zone.clim = NULL, Name.zone = NULL, Temp.zone, Diff.MAAT.MAP = T, Meta.data = NULL, Top.rect = F, Top.ratio = .5, Save.merged.data = NULL,
                            Top.smooth = 20, # Top.smooth in [5;30] c'est bien si le signal de base est mal lissé
                            Use.symbol = F, Display.working.plot = F, Keep.top.from = NULL, Stack.proxy = F, Remove.sites = F,
                            Show.proxy.t = T, Detrending = F, Symb.size = .4, Panel.param = .25){
  #### Init Val ####
  # if(missing(Site.order)){Site.order = names(MFT[[1]][-1])}
  if(missing(Temp.zone)){Temp.zone = rep("C", length(Zone.clim))}
  if(is.null(Name.zone) == F){Site.order <- c("Climate Periods", Site.order)}
  library(RColorBrewer)
  
  #### Zone clim graphical settings ####
  if(length(unique(Temp.zone)) == 0){Title.color.scale <- c("grey")}
  if(length(unique(Temp.zone)) == 1){Title.color.scale <- c("grey")}
  if(length(unique(Temp.zone)) == 3){Title.color.scale <- c("#75AADB", "black", "#E76D51")}
  if(length(unique(Temp.zone)) == 2 & "C" %in% unique(Temp.zone) & "W" %in% unique(Temp.zone)){Title.color.scale <- c("#75AADB","#E76D51")}
  if(length(unique(Temp.zone)) == 2 & "G" %in% unique(Temp.zone) & "W" %in% unique(Temp.zone)){Title.color.scale <- c("black", "#E76D51")}
  if(length(unique(Temp.zone)) == 2 & "G" %in% unique(Temp.zone) & "C" %in% unique(Temp.zone)){Title.color.scale <- c("#75AADB", "black")}
  
  Rect.color.scale <- c("grey40", "grey5")
  
  if(is.null(Zone.clim) == F){
    yo = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                    xmax = Zone.clim[seq(2,length(Zone.clim), by=2)], 
                    Temp.col = Temp.zone)
    
    My_climate_zone <- geom_rect(data = yo, inherit.aes = F, 
                                 mapping = aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=+Inf, fill = Temp.col), 
                                 alpha=0.1, color = NA, size = 0.3, linetype = 2)
    
    My_lines <- unique(c(yo$xmin, yo$xmax))
    My_lines <- geom_vline(xintercept = My_lines, col = "black", lty = 2, alpha = 1)
  } 
  else{
    yo = data.frame(xmin = 0, xmax = 0, Temp.col = "black")
    My_climate_zone <- NULL; My_lines <- NULL
  }
  if(is.null(Name.zone) == F){yo2 = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                                               xmax = Zone.clim[seq(2,length(Zone.clim), by=2)],
                                               Temp.col = Temp.zone,
                                               Title.zone = Name.zone,
                                               L1 = 2.5,
                                               Lake = factor("Climate Periods", levels = Site.order))
  My_climate_name <- geom_text(data = yo2, inherit.aes = F,
                               aes(x = (xmax+xmin)/2, y = L1, label = Title.zone, color = Temp.col), size = 3, vjust = 0.5, fontface = "bold")
  
  
  }
  else{
    yo2 = data.frame(xmin = 0, xmax = 0, Temp.col = "", Title.zone = "")
    My_climate_name = NULL}
  
  #### Merge pollen, gdgt and other ####
  if (is.null(MFT) && is.null(MGDGT) && is.null(MXRF) && is.null(Md18O)) {stop("Import models to plot: at least one of MFT, MGDGT, MXRF must be provided.")}
  rename_proxy <- function(df, proxy) {names(df)[-1] <- paste0(names(df)[-1], "", proxy); return(df)}
  
  n_models <- max(
    length(MFT %||% list()),
    length(MGDGT %||% list()),
    length(MXRF %||% list()),
    length(Md18O %||% list())
  )
  
  MModel <- vector("list", n_models)
  
  for (i in seq_len(n_models)) {
    to_merge <- list()
    if (!is.null(MFT) && length(MFT) >= i) {to_merge <- append(to_merge, list(rename_proxy(MFT[[i]], ".Pollen")))}
    if (!is.null(MGDGT) && length(MGDGT) >= i) {to_merge <- append(to_merge, list(rename_proxy(MGDGT[[i]], "brGDGT")))}
    if (!is.null(MXRF) && length(MXRF) >= i) {to_merge <- append(to_merge, list(rename_proxy(MXRF[[i]], ".XRF")))}
    if (!is.null(Md18O) && length(Md18O) >= i) {to_merge <- append(to_merge, list(rename_proxy(Md18O[[i]], ".d18O")))}
    if (length(to_merge) == 1) {MModel[[i]] <- to_merge[[1]]} 
    else {MModel[[i]] <- Reduce(function(x, y) merge(x, y, by = "Age", all = TRUE), to_merge)}
  }
  
  #### Remove manually sites ####
  if(is.null(Remove.sites) == F){
    for(i in 1:length(MModel)){
      MModel[[i]] <- MModel[[i]][,which(!gsub("\\..*", "", names(MModel[[i]])) %in% Remove.sites)]
    }
  }
  
  #### Detrending ####
  if(Detrending == T){
    for(i in 1:length(MModel)){
      Mp <- MModel[[i]]
      time_col <- "Age"
      signal_cols <- setdiff(names(Mp), time_col)
      
      for (col in signal_cols) {
        if(!all(is.na(Mp[[col]]))){
          model <- resid(lm(Mp[[col]] ~ Mp[[time_col]], na.action = na.exclude))
          # model <- 2*(model-min(model, na.rm = T))/(max(model, na.rm = T)-min(model, na.rm = T))-1
          Mp[[col]] <- model
        }
      }
      MModel[[i]] <- Mp
    }
  }
  #### Calcul difference param ####
  if(Diff.MAAT.MAP == T){
    MModel[[3]] <- MModel[[1]] + MModel[[2]]
    MModel[[3]]$Age <- MModel[[1]]$Age
    names(MModel[[3]]) <- gsub("MAAT", "MAAT + MAP", names(MModel[[1]]))}
  
  #### Stack proxies ####
  if(Stack.proxy == T){
    for(i in 1:length(MModel)){
      fill_group_nas <- function(df, group_cols) {
        for (cols in group_cols) {
          group_data <- df[, cols, drop = FALSE]
          mask <- apply(group_data, 1, function(row) any(!is.na(row)))
          sub <- group_data
          sub[!mask, ] <- NA  # keep rows with all NA untouched
          
          fill_down <- function(x) {
            for (i in 2:length(x)) {
              if (is.na(x[i])) x[i] <- x[i - 1]
            }
            x
          }
          fill_up <- function(x) {
            for (i in (length(x) - 1):1) {
              if (is.na(x[i])) x[i] <- x[i + 1]
            }
            x
          }
          
          for (j in seq_along(cols)) {
            x <- sub[[j]]
            x <- fill_down(x)
            x <- fill_up(x)
            sub[[j]] <- x
          }
          
          group_data[mask, ] <- sub[mask, ]
          df[, cols] <- group_data
        }
        return(df)
      }
      
      group_names <- gsub("\\..*", "", names(MModel[[i]]))
      Param <- gsub("\\.", "", stringr::str_match(names(MModel[[i]]), "\\.\\s*(.*?)\\s*\\."))[,1]
      Proxy <- gsub(".*\\.", "", names(MModel[[i]]))
      
      A <- list()
      for(tt in 1:length(unique(group_names))){A[[tt]] <- which(group_names %in% unique(group_names)[tt])}
      A <- Filter(function(x) length(x) >= 2, A)
      MModel[[i]] <- fill_group_nas(MModel[[i]], A)
      
      New.names <- data.frame(Lake = group_names, Proxy = Proxy, Param = Param)
      proxy_list <- aggregate(Proxy ~ Lake, data = New.names, FUN = function(x) unique(as.character(x)))
      param_list <- aggregate(Param ~ Lake, data = New.names, FUN = function(x) unique(as.character(x)))
      combined <- merge(param_list, proxy_list, by = "Lake")
      
      combined$Param <- sapply(combined$Param, paste, collapse = "|")
      combined$Proxy <- sapply(combined$Proxy, paste, collapse = "|")
      combined$New.name <- apply(combined, 1, function(x) paste(x, collapse = "."))
      aggregated <- t(aggregate(t(MModel[[i]]), by = list(Group = New.names$Lake), FUN = mean, na.rm = T))
      colnames(aggregated) <- aggregated[1, ]  # first row has group names
      aggregated <- aggregated[-1, ]           # remove grouping row
      aggregated <- as.data.frame(aggregated)
      aggregated[] <- lapply(aggregated, as.numeric)  # convert all back to numeric
      names(aggregated)[which(names(aggregated) %in% combined$Lake)] <- combined$New.name
      MModel[[i]] <- aggregated
    }
  }
  
  #### Top rectangles ####
  if(Top.rect == T){
    find_wave_tops_df <- function(df, k = Top.smooth, height_ratio = Top.ratio) {
      time_vec <- df[[1]]  # first column is time
      signal_cols <- df[, -1]  # remaining columns are signals
      
      #### My functions ####
      find_wave_tops_edge <- function(y) {
        if (all(is.na(y))) return(list())
        
        if(is.null(k) == F){smooth_y <- stats::filter(y, rep(1/k, k), sides = 2)}
        else{smooth_y <- y}
        
        valid <- which(!is.na(smooth_y))
        y_clean <- smooth_y[valid]
        
        # Detect local maxima (internal)
        maxima_clean <- which(diff(sign(diff(y_clean))) == -2) + 1
        maxima <- valid[maxima_clean]
        
        # Check for peak at the first valid point
        if (length(valid) >= 2 && y_clean[1] > y_clean[2]) {
          maxima <- c(valid[1], maxima)
        }
        
        # Check for peak at the last valid point
        if (length(valid) >= 2 && y_clean[length(y_clean)] > y_clean[length(y_clean) - 1]) {
          maxima <- c(maxima, valid[length(valid)])
        }
        
        # Find top regions
        bounds <- lapply(maxima, function(i) {
          peak_val <- y[i]
          threshold <- peak_val * height_ratio
          
          left <- i
          while (left > 1 && !is.na(y[left - 1]) && y[left - 1] >= threshold) {
            left <- left - 1
          }
          
          right <- i
          while (right < length(y) && !is.na(y[right + 1]) && y[right + 1] >= threshold) {
            right <- right + 1
          }
          
          return(c(start_time = time_vec[left], end_time = time_vec[right]))
        })
        
        return(bounds)
      }
      
      find_wave_tops_slope_expansion <- function(y, time_vec, k = 10, slope_threshold = 0.01, max_expansion = 20) {
        if (all(is.na(y))) return(list())
        
        # Smooth signal
        if (!is.null(k)) {
          smooth_y <- stats::filter(y, rep(1/k, k), sides = 2)
        } else {
          smooth_y <- y
        }
        
        valid <- which(!is.na(smooth_y))
        y_clean <- smooth_y[valid]
        
        # Find local maxima
        maxima_clean <- which(diff(sign(diff(y_clean))) == -2) + 1
        maxima <- valid[maxima_clean]
        
        bounds <- lapply(maxima, function(i) {
          peak_val <- y[i]
          
          # Expand left
          left <- i
          steps_left <- 0
          while (left > 1 && steps_left < max_expansion) {
            slope_left <- y[left] - y[left - 1]
            # if (is.na(slope_left) || abs(slope_left) < slope_threshold) break
            if (is.na(slope_left) || abs(slope_left) > slope_threshold) break
            left <- left - 1
            steps_left <- steps_left + 1
          }
          
          # Expand right
          right <- i
          steps_right <- 0
          while (right < length(y) && steps_right < max_expansion) {
            slope_right <- y[right + 1] - y[right]
            # if (is.na(slope_right) || abs(slope_right) < slope_threshold) break
            if (is.na(slope_right) || abs(slope_right) > slope_threshold) break
            right <- right + 1
            steps_right <- steps_right + 1
          }
          
          return(c(start_time = time_vec[left], end_time = time_vec[right]))
        })
        
        return(bounds)
      }
      
      find_wave_tops_relative_expansion <- function(y, time_vec, k = 10,
                                                    rel_min = 0.001, rel_max = 0.2,
                                                    max_expansion = 50) {
        if (all(is.na(y))) return(list())
        
        # Smoothing
        if (!is.null(k)) {
          smooth_y <- stats::filter(y, rep(1 / k, k), sides = 2)
        } else {
          smooth_y <- y
        }
        
        valid <- which(!is.na(smooth_y))
        y_clean <- smooth_y[valid]
        
        # Find local maxima (internal)
        maxima_clean <- which(diff(sign(diff(y_clean))) == -2) + 1
        maxima <- valid[maxima_clean]
        
        # Check for peak at the first valid point
        if (length(valid) >= 2 && y_clean[1] > y_clean[2]) {
          maxima <- c(valid[1], maxima)
        }
        
        # Check for peak at the last valid point
        if (length(valid) >= 2 && y_clean[length(y_clean)] > y_clean[length(y_clean) - 1]) {
          maxima <- c(maxima, valid[length(valid)])
        }
        
        bounds <- lapply(maxima, function(i) {
          peak_val <- y[i]
          if (is.na(peak_val) || peak_val == 0) return(NULL)
          
          # Expand left
          left <- i
          steps_left <- 0
          while (left > 1 && steps_left < max_expansion) {
            delta <- y[left] - y[left - 1]
            rel_change <- abs(delta / peak_val)
            if (is.na(rel_change) || rel_change < rel_min || rel_change > rel_max) break
            left <- left - 1
            steps_left <- steps_left + 1
          }
          
          # Expand right
          right <- i
          steps_right <- 0
          while (right < length(y) && steps_right < max_expansion && (right + 1) <= length(y)) {
            delta <- y[right + 1] - y[right]
            rel_change <- abs(delta / peak_val)
            if (is.na(rel_change) || rel_change < rel_min || rel_change > rel_max) break
            right <- right + 1
            steps_right <- steps_right + 1
          }
          
          # Force region to full edge if the peak is at the boundary
          if (i == valid[1]) {
            left <- 1  # Extend fully to start
          }
          if (i == valid[length(valid)]) {
            right <- length(y)  # Extend fully to end
          }
          
          return(c(start_time = time_vec[left], end_time = time_vec[right]))
        })
        
        bounds <- bounds[!sapply(bounds, is.null)]
        return(bounds)
      }
      
      find_wave_tops_balanced <- function(y, time_vec, k = 10,
                                          rel_min = 0.001, rel_max = 0.2,
                                          max_expansion = 50) {
        if (all(is.na(y))) return(list())
        
        # Smoothing
        if (!is.null(k)) {
          smooth_y <- stats::filter(y, rep(1 / k, k), sides = 2)
        } else {
          smooth_y <- y
        }
        
        valid <- which(!is.na(smooth_y))
        y_clean <- smooth_y[valid]
        
        # Find local maxima (internal)
        maxima_clean <- which(diff(sign(diff(y_clean))) == -2) + 1
        maxima <- valid[maxima_clean]
        
        # Check for edge peaks
        if (length(valid) >= 2 && y_clean[1] > y_clean[2]) {
          maxima <- c(valid[1], maxima)
        }
        if (length(valid) >= 2 && y_clean[length(y_clean)] > y_clean[length(y_clean) - 1]) {
          maxima <- c(maxima, valid[length(valid)])
        }
        
        bounds <- lapply(maxima, function(i) {
          peak_val <- y[i]
          if (is.na(peak_val) || peak_val == 0) return(NULL)
          
          left <- i
          right <- i
          
          steps <- 0
          done_left <- FALSE
          done_right <- FALSE
          
          while (steps < max_expansion && (!done_left || !done_right)) {
            # Alternate step
            if (steps %% 2 == 0 && !done_left && left > 1) {
              slope <- abs((y[left] - y[left - 1]) / peak_val)
              if (!is.na(slope) && slope >= rel_min && slope <= rel_max) {
                left <- left - 1
              } else {
                done_left <- TRUE
              }
            } else if (steps %% 2 == 1 && !done_right && right < length(y)) {
              slope <- abs((y[right + 1] - y[right]) / peak_val)
              if (!is.na(slope) && slope >= rel_min && slope <= rel_max) {
                right <- right + 1
              } else {
                done_right <- TRUE
              }
            }
            steps <- steps + 1
          }
          
          # Force start/end if it's an edge peak
          if (i == valid[1]) left <- 1
          if (i == valid[length(valid)]) right <- length(y)
          
          return(c(start_time = time_vec[left], end_time = time_vec[right]))
        })
        
        bounds <- bounds[!sapply(bounds, is.null)]
        return(bounds)
      }
      
      merge_overlapping_regions <- function(bounds) {
        if (length(bounds) == 0) {
          return(data.frame(start_time = numeric(0), end_time = numeric(0)))
        }
        
        # Filter out any invalid or incomplete entries
        bounds <- Filter(function(x) !is.null(x) && length(x) == 2 && all(!is.na(x)), bounds)
        if (length(bounds) == 0) {
          return(data.frame(start_time = numeric(0), end_time = numeric(0)))
        }
        
        df <- as.data.frame(do.call(rbind, bounds))
        colnames(df) <- c("start_time", "end_time")
        
        # Remove rows with NA just in case
        df <- df[complete.cases(df), ]
        if (nrow(df) == 0) {
          return(data.frame(start_time = numeric(0), end_time = numeric(0)))
        }
        
        # Sort by start time
        df <- df[order(df$start_time), , drop = FALSE]
        
        merged <- list()
        current <- df[1, , drop = FALSE]
        
        for (i in 2:nrow(df)) {
          next_start <- df[i, "start_time"]
          next_end   <- df[i, "end_time"]
          cur_end    <- current[1, "end_time"]
          
          if (!is.na(next_start) && !is.na(cur_end) && next_start <= cur_end) {
            current[1, "end_time"] <- max(cur_end, next_end, na.rm = TRUE)
          } else {
            merged[[length(merged) + 1]] <- current
            current <- df[i, , drop = FALSE]
          }
        }
        
        merged[[length(merged) + 1]] <- current
        merged_df <- do.call(rbind, merged)
        colnames(merged_df) <- c("start_time", "end_time")
        return(as.data.frame(merged_df))
      }
      
      #### Apply to each signal column ####
      # result <- lapply(signal_cols, find_wave_tops_edge)
      # result <- lapply(signal_cols, function(y) {find_wave_tops_slope_expansion(y, time_vec, k = 10, slope_threshold = 0.02)})
      # result <- lapply(signal_cols, function(y) {find_wave_tops_relative_expansion(y, time_vec, k = 10, rel_min = 0.0001, rel_max = 0.05, max_expansion = 20)})
      result <- lapply(signal_cols, function(y) {find_wave_tops_balanced(y, time_vec, k = 20, rel_min = 0.0001, rel_max = 0.04, max_expansion = 30)})
      
      names(result) <- names(signal_cols)
      result <- lapply(result, function(bounds_list) {
        if (length(bounds_list) == 0) return(NULL) 
        merge_overlapping_regions(bounds_list)})
      
      return(result)
    }
    
    find_wave_tops_df2 <- function(df, k = Top.smooth, height_ratio = Top.ratio) {
      time_vec <- df[[1]]  # first column is time
      signal_cols <- df[, -1]  # remaining columns are signals
      
      adjusted_height_ratio <- function(i, len, base_ratio = 0.3) {
        edge_buffer <- floor(len * 0.1)
        if (i <= edge_buffer || i >= (len - edge_buffer)) {
          return(base_ratio * 0.5)
        }
        return(base_ratio)
      }
      
      merge_overlapping_regions <- function(bounds) {
        if (length(bounds) == 0) {
          return(data.frame(start_time = numeric(0), end_time = numeric(0)))
        }
        
        # Filter out any invalid or incomplete entries
        bounds <- Filter(function(x) !is.null(x) && length(x) == 2 && all(!is.na(x)), bounds)
        if (length(bounds) == 0) {
          return(data.frame(start_time = numeric(0), end_time = numeric(0)))
        }
        
        df <- as.data.frame(do.call(rbind, bounds))
        colnames(df) <- c("start_time", "end_time")
        
        # Remove rows with NA just in case
        df <- df[complete.cases(df), ]
        if (nrow(df) == 0) {
          return(data.frame(start_time = numeric(0), end_time = numeric(0)))
        }
        
        # Sort by start time
        df <- df[order(df$start_time), , drop = FALSE]
        
        merged <- list()
        current <- df[1, , drop = FALSE]
        
        for (i in 2:nrow(df)) {
          next_start <- df[i, "start_time"]
          next_end   <- df[i, "end_time"]
          cur_end    <- current[1, "end_time"]
          
          if (!is.na(next_start) && !is.na(cur_end) && next_start <= cur_end) {
            current[1, "end_time"] <- max(cur_end, next_end, na.rm = TRUE)
          } else {
            merged[[length(merged) + 1]] <- current
            current <- df[i, , drop = FALSE]
          }
        }
        
        merged[[length(merged) + 1]] <- current
        merged_df <- do.call(rbind, merged)
        colnames(merged_df) <- c("start_time", "end_time")
        return(as.data.frame(merged_df))
      }
      
      find_wave_tops_core <- function(y, edge_sensitive = TRUE) {
        if (all(is.na(y))) return(data.frame())
        
        if (!is.null(k)) {
          smooth_y <- stats::filter(y, rep(1 / k, k), sides = 2)
        } else {
          smooth_y <- y
        }
        
        valid <- which(!is.na(smooth_y))
        y_clean <- smooth_y[valid]
        if (length(valid) < 3) return(data.frame())
        
        maxima_clean <- which(diff(sign(diff(y_clean))) == -2) + 1
        maxima <- valid[maxima_clean]
        
        # Add edge peaks if needed
        if (edge_sensitive) {
          if (y_clean[1] > y_clean[2]) maxima <- c(valid[1], maxima)
          if (y_clean[length(y_clean)] > y_clean[length(y_clean) - 1]) {
            maxima <- c(maxima, valid[length(valid)])
          }
        }
        
        # Define top bounds
        bounds <- lapply(maxima, function(i) {
          peak_val <- y[i]
          h_ratio <- adjusted_height_ratio(i, length(y), height_ratio)
          threshold <- peak_val * h_ratio
          
          left <- i
          while (left > 1 && !is.na(y[left - 1]) && y[left - 1] >= threshold) {
            left <- left - 1
          }
          
          right <- i
          while (right < length(y) && !is.na(y[right + 1]) && y[right + 1] >= threshold) {
            right <- right + 1
          }
          
          return(c(time_vec[left], time_vec[right]))
        })
        # print(bounds)
        # Merge overlapping regions
        merged_bounds <- merge_overlapping_regions(bounds)
        # print(merged_bounds)
        return(merged_bounds)
      }
      
      result <- lapply(signal_cols, function(y) find_wave_tops_core(y, edge_sensitive = TRUE))
      names(result) <- names(signal_cols)
      empty_flags <- sapply(result, function(x) is.data.frame(x) && nrow(x) == 0 && ncol(x) == 0)
      result <- result[!empty_flags]
      return(result)
    }
    
    find_wave_tops_by_bottoms <- function(df, k = Top.smooth, min_width = NULL) {
      time_vec <- df[[1]]  # first column is time
      signal_cols <- df[, -1]  # remaining columns are signals
      
      # Internal helper: detect local minima and maxima with smoothing
      find_local_extrema <- function(y, k) {
        if (all(is.na(y))) return(list(tops = numeric(0), bottoms = numeric(0)))
        
        # Smooth
        smooth_y <- if (!is.null(k)) stats::filter(y, rep(1/k, k), sides = 2) else y
        valid <- which(!is.na(smooth_y))
        y_clean <- smooth_y[valid]
        
        # Local maxima and minima in the clean signal
        diff_sign <- sign(diff(y_clean))
        maxima_clean <- which(diff(diff_sign) == -2) + 1
        minima_clean <- which(diff(diff_sign) == 2) + 1
        
        maxima <- valid[maxima_clean]
        minima <- valid[minima_clean]
        
        # Edge checks
        if (length(valid) >= 2) {
          if (y_clean[1] > y_clean[2]) maxima <- c(valid[1], maxima)
          if (y_clean[1] < y_clean[2]) minima <- c(valid[1], minima)
          if (y_clean[length(y_clean)] > y_clean[length(y_clean) - 1]) maxima <- c(maxima, valid[length(valid)])
          if (y_clean[length(y_clean)] < y_clean[length(y_clean) - 1]) minima <- c(minima, valid[length(valid)])
        }
        
        return(list(tops = maxima, bottoms = minima))
      }
      
      # Peak region estimation based on surrounding bottoms
      calculate_peak_regions <- function(tops, bottoms, time_vec) {
        bounds <- list()
        
        for (peak in tops) {
          prev_bottoms <- bottoms[bottoms < peak]
          next_bottoms <- bottoms[bottoms > peak]
          
          if (length(prev_bottoms) == 0 || length(next_bottoms) == 0) next  # skip edge peaks
          
          left <- prev_bottoms[length(prev_bottoms)]
          right <- next_bottoms[1]
          
          # Calculate average width and center on peak
          width <- as.numeric(time_vec[right] - time_vec[left])
          half_width <- width / 20
          
          center_time <- as.numeric(time_vec[peak])
          start_time <- center_time - half_width
          end_time <- center_time + half_width
          
          # Optionally skip narrow peaks
          if (!is.null(min_width) && width < min_width) next
          
          bounds[[length(bounds) + 1]] <- c(start_time = start_time, end_time = end_time)
        }
        
        return(bounds)
      }
      
      merge_overlapping_regions <- function(bounds) {
        if (length(bounds) == 0) {
          return(data.frame(start_time = numeric(0), end_time = numeric(0)))
        }
        
        # Filter out any invalid or incomplete entries
        bounds <- Filter(function(x) !is.null(x) && length(x) == 2 && all(!is.na(x)), bounds)
        if (length(bounds) == 0) {
          return(data.frame(start_time = numeric(0), end_time = numeric(0)))
        }
        
        df <- as.data.frame(do.call(rbind, bounds))
        colnames(df) <- c("start_time", "end_time")
        
        # Remove rows with NA just in case
        df <- df[complete.cases(df), ]
        if (nrow(df) == 0) {
          return(data.frame(start_time = numeric(0), end_time = numeric(0)))
        }
        
        # Sort by start time
        df <- df[order(df$start_time), , drop = FALSE]
        
        merged <- list()
        current <- df[1, , drop = FALSE]
        
        for (i in 2:nrow(df)) {
          next_start <- df[i, "start_time"]
          next_end   <- df[i, "end_time"]
          cur_end    <- current[1, "end_time"]
          
          if (!is.na(next_start) && !is.na(cur_end) && next_start <= cur_end) {
            current[1, "end_time"] <- max(cur_end, next_end, na.rm = TRUE)
          } else {
            merged[[length(merged) + 1]] <- current
            current <- df[i, , drop = FALSE]
          }
        }
        
        merged[[length(merged) + 1]] <- current
        merged_df <- do.call(rbind, merged)
        colnames(merged_df) <- c("start_time", "end_time")
        return(as.data.frame(merged_df))
      }
      
      
      #### Main processing ####
      result <- lapply(signal_cols, function(y) {
        extrema <- find_local_extrema(y, k)
        if (length(extrema$tops) == 0 || length(extrema$bottoms) < 2) return(NULL)
        calculate_peak_regions(extrema$tops, extrema$bottoms, time_vec)
      })
      
      names(result) <- names(signal_cols)
      
      # Merge overlapping peaks per signal
      result <- lapply(result, function(bounds_list) {
        if (length(bounds_list) == 0) return(NULL)
        merge_overlapping_regions(bounds_list)
      })
      
      return(result)
    }
    
    Top.list <- list()
    for(i in 1:length(MModel)){
      # tops_bounds <- find_wave_tops_df(MModel[[i]])
      tops_bounds <- find_wave_tops_by_bottoms(MModel[[i]], k = Top.smooth, min_width = 0)
      
      tops_df <- do.call(rbind,
                         Filter(Negate(is.null), lapply(names(tops_bounds), function(signal) {
                           df <- tops_bounds[[signal]]
                           if (is.null(df) || nrow(df) == 0) return(NULL)
                           data.frame(signal = signal, start_time = df$start, end_time = df$end)})))
      
      rownames(tops_df) <- NULL
      tops_df$Proxy <- gsub(".*\\.", "", tops_df$signal)
      tops_df$Lake <- gsub("\\..*", "", tops_df$signal)
      tops_df <- na.omit(tops_df)
      tops_df <- tops_df[tops_df$start_time != tops_df$end_time,]
      
      Top.list[[i]] <- tops_df
      names(Top.list)[[i]] <- unique(gsub("\\.", "", stringr::str_match(tops_df$signal, "\\.\\s*(.*?)\\s*\\."))[,1])
    }
  }
  else{Top.list <- NULL}
  
  if(is.null(Keep.top.from) == F){
    Val.to.keep <- Top.list[names(Top.list) == Keep.top.from]
    Top.list[names(Top.list) != Keep.top.from] <- Val.to.keep
  }
  
  #### Recup des stats par période ####
  Stats.phase.clim.param <- list()
  if(is.null(My_climate_name) == F){
    Stats.phase.clim <- data.frame(matrix(nrow = nrow(yo2), ncol = (ncol(MModel[[1]])-1), NA))
    row.names(Stats.phase.clim) <- yo2$Title.zone
    
    for(i in 1:length(MModel)){
      Mp <- MModel[[i]]
      colnames(Stats.phase.clim) <- gsub("\\..*", "", names(MModel[[i]])[2:ncol(MModel[[i]])])
      
      for(j in 1:nrow(yo2)){
        Stats.phase.clim[j,] <- colMeans(Mp[which(Mp$Age <= yo2$xmax[j] & Mp$Age >= yo2$xmin[j]),2:ncol(Mp)], na.rm = T)
        
      }
      
      Stats.phase.clim.param[[i]] <- t(Stats.phase.clim)
      if(is.null(Meta.data) == F){
        Keep <- Stats.phase.clim.param[[i]]
        Common.site <- intersect(row.names(Keep), row.names(Meta.data))
        Stats.phase.clim.param[[i]] <- cbind(Meta.data[Common.site,], Keep[Common.site,])
      }
      # print(Stats.phase.clim.param)
    }
    if(Diff.MAAT.MAP == T){names(Stats.phase.clim.param) <- c("MAAT", "MAP", "DeltaParam")}
    if(Diff.MAAT.MAP == F){names(Stats.phase.clim.param) <- c("MAAT", "MAP")}}
  
  #### Data preparations ####
  Mfull.keep <- MModel
  MModel <- melt(MModel, id = "Age")
  MModel$Lake <- gsub("\\..*", "", MModel$variable)
  MModel$Lake <- factor(MModel$Lake, levels = Site.order)
  MModel$Param <- gsub("\\.", "", stringr::str_match(MModel$variable, "\\.\\s*(.*?)\\s*\\."))[,1]
  MModel$Proxy <- gsub(".*\\.", "", MModel$variable)
  MModel$L1 <- 1
  MModel$L1[MModel$Proxy == "Pollen"] <- 1
  MModel$L1[MModel$Proxy == "brGDGT"] <- 2
  MModel$L1[MModel$Proxy == "XRF"] <- 3
  MModel$L1[MModel$Proxy == "d18O"] <- 4
  
  #### Graph settings ####
  if(missing(Limites)){Limites = c(min(MModel$Age), max(MModel$Age))}
  Ylab.size <- c(rep(0, (length(unique(MModel$Param))-1)),9.5)
  A = 1:11
  my_orange = list(brewer.pal(n = 11, "RdBu")[A[-c(4,5,6,7,8)]],#[A[-c(4,5,6,7,8)]],
                   rev(brewer.pal(n = 11, "PuOr")[A[-c(4,5,6,7,8)]]),
                   brewer.pal(n = 11, "Spectral")[A[-c(1,3,4,5,7,8,9,11)]])
  
  if(is.null(Zone.clim) == F){Panel.size2 <- ggh4x::force_panelsizes(rows = c(Panel.param, rep(1, length(unique(MModel$Lake)))))}
  else{Panel.size2 <- NULL}
  
  #### Working plot for climate-tops ####
  if(Display.working.plot == T & is.null(Top.list) == F){
    devices <- dev.list()
    if (is.null(devices) || !"X11cairo" %in% names(devices)) {X11()}
    
    tops_test <- dplyr::bind_rows(lapply(names(Top.list), function(var) {Top.list[[var]] %>% mutate(Param = var)}))
    p <- ggplot(MModel, aes(x = Age, y = value, color = interaction(Lake, Proxy)))+ geom_line()+
      facet_grid(Param ~ Lake, scales = "free")+
      geom_rect(data = tops_test, inherit.aes = F, 
                mapping = aes(xmin = start_time, xmax = end_time, ymin = 1, ymax = 1.01, color = interaction(Lake, Proxy)),
                size = 0.5, linetype = 1, fill = NA)+
      theme(legend.position = "none")
    print(p)
  }
  
  #### Main loop ####
  Plot.list <- list()
  for(i in 1:length(unique(MModel$Param))){
    Param.selected <- unique(MModel$Param)[i]
    Mplot <- MModel[MModel$Param == Param.selected,]
    # print(seq(Limites[1], Limites[2], 50))
    if(Show.proxy.t == T){Alf <- c("(B) ", "(C) " , "(D) ")}
    else{Alf <- c("(A) ", "(B) " , "(C) ")}
    Title.subplot <- paste(Alf[i], Param.selected, sep = " ")
    
    #### Color plot ####
    orange_palette = colorRampPalette(my_orange[[i]])
    my_orange.p = rev(orange_palette(length(seq(-10, 10, by = 1))))
    
    #### Plot top rectangles ####
    if(Top.rect == T){
      tops_df <- Top.list[[i]]
      tops_df$Lake <- factor(tops_df$Lake, levels = levels(Mplot$Lake))
      tops_df$L1 <- 1
      tops_df$L1[tops_df$Proxy == "Pollen"] <- 1
      tops_df$L1[tops_df$Proxy == "brGDGT"] <- 2
      tops_df$L1[tops_df$Proxy == "XRF"] <- 3
      tops_df$L1[tops_df$Proxy == "d18O"] <- 4
      
      My_top_zone <- geom_rect(data = tops_df, inherit.aes = F, 
                               # mapping = aes(xmin = start_time, xmax = end_time), ymin = -Inf, ymax = +Inf,
                               mapping = aes(xmin = start_time, xmax = end_time, ymin = L1-0.5, ymax = L1+0.5),
                               fill = NA, alpha = 0.2, color = "black", size = 0.5, linetype = 1)
    }
    else{My_top_zone <- NULL}
    
    #### Plot ####
    sig <- (Mplot$Age[2] - Mplot$Age[1])/2
    Plot.list[[i]]  <-  ggplot(Mplot)+
      #### Items ####
    geom_rect(mapping = aes(xmin = Age - sig, xmax = Age + sig, ymin = L1-0.5, ymax = L1+0.5, colour = value), na.rm = T)+
      geom_rect(mapping = aes(xmin = Age - sig, xmax = Age + sig, ymin = L1-0.5, ymax = L1+0.5, fill = value), na.rm = T)+
      ggtitle(Title.subplot)+
      facet_grid(Lake ~ ., scales = "free") + Panel.size2 +
      scale_color_gradientn(colours = my_orange.p,  na.value = "gray90")+
      scale_fill_gradientn(colours = my_orange.p,  na.value = "gray90", guide = "none")+
      scale_y_continuous(name = NULL)+
      scale_x_continuous(name = "Time (cal. year BP)", limits = c(Limites[1], Limites[2]), n.breaks = 10, minor_breaks = seq(Limites[1], Limites[2], 50), expand = c(0,0))+
      
      #### Zone clim ####
    new_scale_color()+
      new_scale_fill()+
      scale_fill_manual(values = Rect.color.scale, guide = "none")+
      scale_color_manual(values = Title.color.scale, guide = "none", name = NULL, labels = NULL, breaks = NULL, na.translate = FALSE)+
      
      My_climate_zone +
      My_lines +
      My_climate_name +
      My_top_zone +
      #### Theme ####
    # theme_bw()+
    theme(axis.ticks.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5),
          legend.title = element_blank(),
          axis.line.x = element_line(),
          panel.spacing = unit(0, "cm"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          legend.direction = "horizontal",
          legend.position = "top",
          plot.title.position = "plot",
          legend.justification = c("center"),               # left, top, right, bottom
          legend.text = element_text(size = 8),
          strip.text.y = element_text(angle = 0, hjust = 0, size = Ylab.size[i]),
          plot.margin = unit(c(0,0,0,0), "lines")
    )
  }
  
  #### Save plots and export ####
  Formula.patch <- paste(paste("Plot.list[[", seq(length(unique(MModel$Param))), "]]", sep = ""), collapse = " + ")
  Ptot <- eval(parse(text = Formula.patch))
  
  #### Add plot annotation ####
  Mplot <- Mplot[,4:7]
  Mplot <- Mplot[which(duplicated.data.frame(Mplot)==F),]
  
  if(is.null(My_climate_name) == F){
    Mplot <- rbind(Mplot, c(1, "Climate Periods", NA, NA))}
  
  if(Show.proxy.t == T){
    if(Stack.proxy == T){
      Mplot <- tidyr::separate_rows(Mplot, Proxy, sep = "\\|")
      
      Mplot <- Mplot %>%
        group_by(Lake, Param) %>%
        mutate(Proxy_n = n(), L1 = ifelse(Proxy_n == 1, 1, seq_len(Proxy_n) / Proxy_n))
    }
    if(Use.symbol == T){
      library(ggimage)
      SPol = "Figures/ACA/Trait/Workflow/Symboles/svg/Symbole_pollen.svg"
      SGDGT = "Figures/ACA/Trait/Workflow/Symboles/svg/Symbole_GDGT_2.svg"
      SXRF = "Figures/ACA/Trait/Workflow/Symboles/svg/Symbole_XRF.svg"
      Sd18O = "Figures/ACA/Trait/Workflow/Symboles/svg/Symbole_speleothem.svg"
      
      Mplot$image <- ifelse(Mplot$Proxy == "Pollen", SPol,
                            ifelse(Mplot$Proxy == "brGDGT", SGDGT, 
                                   ifelse(Mplot$Proxy == "XRF", SXRF, Sd18O)))
      Dots.proxy <- geom_image(aes(image = image), size = Symb.size)
    }
    else{Dots.proxy <- geom_point(mapping = aes(colour = Proxy, shape = Proxy), size = 3)}
    
    Padd  <-  ggplot(Mplot, mapping = aes(x = 1, y = L1)) + Dots.proxy +
      #### Items ####
    facet_grid(Lake ~ ., scales = "free") + Panel.size2 +
      scale_y_discrete(name = NULL)+
      scale_x_continuous(name = NULL)+
      ggtitle("(A) Proxy")+
      #### Theme ####
    theme_bw()+
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.title = element_blank(),
            panel.spacing = unit(0, "cm"),
            plot.title = element_text(face = "bold"),
            plot.title.position = "plot",
            strip.background = element_blank(),
            panel.border = element_blank(),
            panel.grid = element_blank(),
            legend.direction = "horizontal",
            legend.position = "top",
            legend.justification = c("center"),               # left, top, right, bottom
            legend.text = element_text(size = 8),
            strip.text.y = element_blank(),
            plot.margin = unit(c(0,0,0,0), "lines")
            
      )
    
    
    Ptot <- Padd + Ptot + plot_layout(nrow = 1, widths = c(1, 20))
  }
  
  #### Export plots ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(Ptot, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", limitsize = F)}
    else{ggsave(Save.plot)}}
  
  if(is.null(Save.merged.data) == F){saveRDS(Mfull.keep, Save.merged.data)}
  
  if(is.null(Save.plot.Rds) == F){saveRDS(Ptot, Save.plot.Rds)}
  
  if(is.null(Save.top) == F){
    Top.list <- lapply(Top.list, function(x) x[-c(1)])
    saveRDS(Top.list, Save.top)
  }
  if(is.null(Save.Rds) == F){saveRDS(Stats.phase.clim.param, Save.Rds)}
  # print(Stats.phase.clim.param)
  return(Stats.phase.clim.param)
  
}

Plot.multiprox <- function(Mlake, Order.alpha, Age.depth, Limites, Title, Smooth.param, Param.clim = NULL, Zone.clim.box = T, Show.x.axis = T, 
                           Reverse.y.list = c("dO18", "Delta"), CRUTS = NULL, CRUTS.display.panel = NULL,
                           Axis.x.begin.zero = T,
                           Scale.y.free = F, Clean.labels = F, Add.parenthesis = F, Add.CI = F, Alpha.CI = .2, Alpha.zone.clim = .1, Rotate.panel.tit = F, Size.panel.tit = 10,
                           Show.line = T, Facet.tit.color = F, Show.smooth = T, Smooth.color = F, Smooth.SD = F, Panel.tit.just = "center",
                           Manual.y.val = NULL, Surf.val = NULL, Temp.zone, Color.by.proxy = T, Manual.sort = F, Panel.space.y = 0, Alternate.y.pos = T, Detrending = F, Bar.size = 2,
                           Name.zone = NULL, Zone.clim = NULL, Event.clim = NULL, Zone.rotation = 0, Leg.pos, Select.interv = 1000, Title.age = "Time (cal. year BP)", H, W, Save.plot){
  #### Init val ####
  if(missing(Order.alpha)){Order.alpha = T}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Limites)){Limites = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Smooth.param)){Smooth.param = 2}
  if(missing(Title)){Title = NULL}
  if(missing(Leg.pos)){Leg.pos = "rigth"}
  if(missing(Temp.zone)){Temp.zone = rep("C", length(Zone.clim))}
  
  #### Temp zone ####
  if(any(unique(grepl("#", Temp.zone))) == T){
    print("Manual color scale for climate zone activated.")
    Rect.color.scale <- Temp.zone
    names(Rect.color.scale) <- Rect.color.scale
  }
  else{
    Rect.color.scale <- c("W" = "#E76D51",
                          "C" = "#75AADB",
                          "G" = "grey40",
                          "D" = "#c67f05",
                          "Wt" = "#004266")
    Rect.color.scale <- Rect.color.scale[unique(Temp.zone)]
  }
  
  if(is.null(Zone.clim) == F){
    if(Zone.clim.box == T){Box.col <- "grey"}
    else{Box.col <- NA}
    
    Rect_scale <- scale_fill_manual(values = Rect.color.scale)
    Rect_scale2 <- scale_color_manual(values = Rect.color.scale)
    
    yo = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                    xmax = Zone.clim[seq(2,length(Zone.clim), by=2)], 
                    Temp.col = Temp.zone)
    Zone.clim.to.plot <- geom_rect(data = yo, inherit.aes = F, na.rm = T, show.legend = F,
                                   mapping = aes(ymin = -Inf, ymax= +Inf, xmin = xmin, xmax = xmax, fill = Temp.col),
                                   alpha = Alpha.zone.clim, color = Box.col, linewidth = 0.3, linetype = 2)}
  
  
  else{Zone.clim.to.plot <- NULL; Rect_scale <- NULL; Rect_scale2 <- NULL}
  
  if(is.null(Name.zone) == F){
    yo2 = data.frame(xmin = Zone.clim[seq(1,length(Zone.clim), by=2)], 
                     xmax = Zone.clim[seq(2,length(Zone.clim), by=2)],
                     Temp.col = Temp.zone,
                     Title.zone = Name.zone,
                     variable = "Height[CWM]")
    Zone.txt1 <- geom_text(data = yo2, inherit.aes = F, aes(x = (xmax+xmin)/2, y = 0, label = Title.zone, color = Temp.col), size = 3, angle = Zone.rotation, vjust = 0.5, hjust = 0.5, fontface = "bold")
  }
  else{Zone.txt1 <- NULL}
  
  if(is.null(Event.clim) == F){Event.clim <- geom_vline(xintercept = Event.clim, color = "grey70", linetype = "dashed")}
  
  #### Detrending ####
  if(Detrending == T){
    for(i in 1:length(Mlake)){
      df <- Mlake[[i]]
      signal_col <- setdiff(names(df), Age.depth)
      detrended <- resid(lm(df[[signal_col]] ~ df[[Age.depth]], na.action = na.exclude))
      # loess_fit <- loess(df[[signal_col]] ~ df[[Age.depth]], span = 0.3)
      # detrended <- df[[signal_col]] - predict(loess_fit)
      
      df[[signal_col]] <- detrended
      Mlake[[i]] <- df
    }
  }
  
  #### Melt ####
  Lake.list <- names(Mlake)
  
  if(Add.CI == F){
    Mmelt <- data.frame(Age.depth = NA, variable = NA, value = NA, lake = NA, lake.lab = NA, Proxy = NA)
  }
  else{Mmelt <- data.frame(Age.depth = NA, variable = NA, value = NA, lake = NA, lake.lab = NA, Proxy = NA, Max = NA, Min = NA)}
  
  for(i in 1:length(Mlake)){
    if(Add.CI == F){A <- melt(Mlake[[i]], id = Age.depth)}
    else{
      A <- melt(Mlake[[i]], id = c(Age.depth, names(Mlake[[i]])[grep("Min|Max", names(Mlake[[i]]))]))
      names(A)[grep("Min", names(A))] <- "Min"
      names(A)[grep("Max", names(A))] <- "Max"
    }
    
    names(A)[1] <- "Age.depth"
    A$lake <- names(Mlake[i])
    A$lake.lab <- gsub("\\..*", "", names(Mlake[i]))
    A$Proxy <- gsub(".*\\.", "", names(Mlake[i]))
    Mmelt <- rbind(Mmelt, A)
  }
  
  #### Clean Lake, proxy and method names ####
  Mmelt$Proxy <- sub("^[^.]*\\.", "", Mmelt$lake)
  Mmelt$lake <- sub("\\.", "\n", Mmelt$lake)
  # print(Mmelt)
  
  #### Manual order of facets ####
  if(Manual.sort == T){
    Mmelt$lake <- factor(Mmelt$lake, unique(Mmelt$lake), ordered = T)
  }
  else{Mmelt$lake <- factor(Mmelt$lake, ordered = T)}
  
  #### Add CI / lines ####
  if(Add.CI == F){CI <- NULL}
  else{CI <- geom_ribbon(aes(ymin = Min, ymax = Max, fill = Proxy), color = NA, alpha = Alpha.CI, show.legend = F)}
  
  if(Show.line == T){My_lines <- geom_line()}
  else{My_lines <- NULL}
  
  #### Add surface values ####
  if(is.null(Surf.val) == F){
    MM <- unique(Mmelt[c("lake", "lake.lab")])
    # MM <- unique(Mmelt[c("variable", "lake", "lake.lab")])
    MM <- na.omit(MM)
    MM$surf.val <- Surf.val[match(MM$lake.lab, Surf.val$Lake), Param.clim]
    
    # #### Test pour faire MAAT et MAF ensemble...#####
    # if(length(Param.clim) == 1){MM$surf.val <- Surf.val[match(MM$lake.lab, Surf.val$Lake), Param.clim]}
    # else{
    #   MM$surf.val <- NA
    #   for(i in Param.clim){
    #     print(i)
    #     # Param.T <- grepl(i, MM$variable)
    #     print(MM$surf.val[grepl(i, MM$variable)])
    #     MM$surf.val[grepl(i, MM$variable)] <- Surf.val[match(MM$lake.lab[grepl(i, MM$variable)], Surf.val$Lake), i]
    #       
    #   }
    # }
    # 
    # print(MM)
    
    Surf_lines <- geom_hline(data = MM, mapping = aes(yintercept = surf.val), col = "grey30", lty = 2, alpha = 0.7)
    
    # print(Surf.val)
  }
  else{Surf_lines <- NULL}
  
  #### Add smooth ####
  if(Show.smooth == T){
    if(Smooth.color == F){
      My_smooth <- geom_smooth(method = "loess", se = F, fullrange = F, level = 0.95, linetype="solid", formula = "y ~ x",
                               linewidth = .7, show.legend = F,  span = Smooth.param, alpha = 0.2, aes(color = Loess))
    }
    else{
      My_smooth <- geom_smooth(method = "loess", se = Smooth.SD, fullrange = F, level = 0.95, linetype="solid", formula = "y ~ x",
                               linewidth = .7, show.legend = F,  span = Smooth.param, alpha = 0.2, aes(color = Proxy, fill = Proxy))
      
    }
  }
  else{My_smooth <- NULL}
  
  #### Color scales ####
  if(Color.by.proxy == T){
    Mmelt$Proxy <- gsub("\n.*", "", Mmelt$Proxy)
    My_color <- c("Pollen" = "#BEA33A",
                  "brGDGT" = "#CF5612", "GDGT" = "#CF5612",
                  "XRF" = "#35978F", "Lake_level" = "#485BB3",
                  "Speleothem" = "#00BA38",
                  "Normal" = "grey40",
                  "Summer" = "darkred",
                  "Spring" = "darkorange",
                  "Winter" = "royalblue",
                  "Annual" = "darkgreen",
                  "Loess" = "#755841", "Paleodust" = "#755841", "NPP" = "#755841",
                  "Tree rings" = "brown",
                  "n-alkane" = "darkgreen", "d18O" = "darkgreen", "d13C" = "darkgreen",
                  "TraCE-21Ka" = "#619CFF",
                  "CRU-TS" = "#619CFF", "HYDE" = "#619CFF",
                  "Ice core" = "#619CFF"
    )
    My_color_scale <- scale_color_manual(values = My_color)
    My_fill_scale <- scale_fill_manual(values = My_color)
    
    Color.facet <- My_color
    Color.facet[which(names(Color.facet) == "Pollen")] <- "#7b6514ff"
    Color.facet[which(names(Color.facet) %in% c("GDGT", "brGDGT"))] <- "#b05e30ff"
    Color.facet[which(names(Color.facet) %in% c("Loess", "Paleodust"))] <- "#593a21ff"
    Color.facet[which(names(Color.facet) == "Speleothem")] <- "#00952dff"
    Color.facet[which(names(Color.facet) == "Lake_level")] <- "#3e51a5ff"
    
  }
  else{My_color_scale <- NULL; My_fill_scale <- NULL}
  
  #### Param graphs ####
  if(is.null(Title) == F){
    if(is.null(Name.zone) == T){Title.1 <- ggtitle(Title)}
    else{Title.2 <- ggtitle(Title); Title.1 <- NULL}
    
  }
  else{Title.1 <- NULL; Title.2 <- NULL}
  
  if(is.null(Limites) == F){Mmelt <- Mmelt[Mmelt$Age.depth > Limites[1] & Mmelt$Age.depth < Limites[2],]}
  
  if(Axis.x.begin.zero == T){Y.zero = 0}
  else{Y.zero <- Limites[1]}
  
  if(Order.alpha == F){
    Mmelt <- Mmelt[-c(1),]
    Mmelt$lake <- factor(Mmelt$lake, levels = names(Mlake), ordered = T)
  }
  
  Mmelt <- Mmelt[!is.na(Mmelt$Age.depth),]  
  
  #### Settings shape ####
  Mmelt$Shape <- "Normal"
  Laskar <- c("Summer", "Winter", "Spring", "Automn", "Annual")
  Mmelt$Shape[Mmelt$Proxy %in% Laskar] <- NA
  
  Mmelt$Loess <- "Normal"
  Mmelt$Loess[Mmelt$Proxy %in% Laskar[1]] <- Laskar[1]
  Mmelt$Loess[Mmelt$Proxy %in% Laskar[2]] <- Laskar[2]
  Mmelt$Loess[Mmelt$Proxy %in% Laskar[3]] <- Laskar[3]
  
  #### Graphical settings ####
  if(Show.x.axis == T){
    Axe.x.param <- element_text(angle = 45, hjust = 1, vjust = 1)
    Axe.x.title <- element_text()
    Axe.x.line <- element_line()
    Axe.x.ticks <- element_line()
  }
  else{
    Axe.x.title <- element_blank()
    Axe.x.param <- element_blank()
    Axe.x.line <- element_blank()
    Axe.x.ticks <- element_blank()
  }
  
  if(Rotate.panel.tit == T){Ang = 90}
  else{Ang = 0}
  
  if(Panel.tit.just == "center"){yhjust <- 0.5}
  if(Panel.tit.just == "left"){yhjust <- 0}
  if(Panel.tit.just == "right"){yhjust <- 1}
  
  #### Alternate left/right axis (or not) ####
  if(Alternate.y.pos == T){
    if(Scale.y.free == F){
      Axis.matrix <- list(scale_y_continuous(), scale_y_continuous(position = "right"))
      Axis.matrix <- Axis.matrix[1:length(unique(Mmelt$lake)) %% 2 +1]
    }
    else{
      if(Add.CI == F){
        limits_df <- Mmelt %>%
          dplyr::group_by(lake) %>%
          dplyr::summarise(ymin = min(value, na.rm = TRUE), ymax = max(value, na.rm = TRUE))}
      else{
        limits_df <- Mmelt %>%
          dplyr::group_by(lake) %>%
          dplyr::summarise(ymin = min(min(value, na.rm = TRUE), min(Min, na.rm = TRUE)),
                           ymax = max(max(value, na.rm = TRUE), max(Max, na.rm = TRUE)))
      }
      
      Axis.matrix <- lapply(seq_len(nrow(limits_df)), function(i) {
        lims <- c(limits_df$ymin[i], limits_df$ymax[i])
        if (i %% 2 == 1) {scale_y_continuous(position = "left", limits = lims)} 
        else {scale_y_continuous(position = "right", limits = lims)}})
    }  
  }
  else{Axis.matrix <- rep(list(scale_y_continuous()), length(unique(Mmelt$lake)))}
  
  if(is.null(Reverse.y.list) == F){
    matches <- grep(paste(Reverse.y.list, collapse = "|"), levels(Mmelt$lake))
    
    Axis.matrix[matches[matches %% 2 == F]] <- list(scale_y_reverse())
    Axis.matrix[matches[matches %% 2 == T]] <- list(scale_y_reverse(position = "right"))
  }
  # if(length(grep("dO18", levels(Mmelt$lake))) >= 1){
  #   Axis.matrix[grep("dO18", levels(Mmelt$lake))[grep("dO18", levels(Mmelt$lake)) %% 2 == F]] <- list(scale_y_reverse())
  #   Axis.matrix[grep("dO18", levels(Mmelt$lake))[grep("dO18", levels(Mmelt$lake)) %% 2 == T]] <- list(scale_y_reverse(position = "right"))
  # }
  # 
  # if(any(grepl("5Me", levels(Mmelt$lake)) & grepl("6Me", levels(Mmelt$lake)))){
  #   Pos.Delta <- which(grepl("5Me", levels(Mmelt$lake)) & grepl("6Me", levels(Mmelt$lake)))
  #   Axis.matrix[Pos.Delta[Pos.Delta %% 2 == F]] <- list(scale_y_reverse())
  #   Axis.matrix[Pos.Delta[Pos.Delta %% 2 == T]] <- list(scale_y_reverse(position = "right"))
  # }
  
  #### Labeller (facet_wrap) ####
  if(Clean.labels == T){
    lake_labeller <- function(x) {
      sapply(x, function(lbl) {
        lbl <- gsub(" ", "~", lbl)
        lbl <- gsub("d18O", "delta^{18}*O", lbl)
        lbl <- gsub("d13C", "delta^{13}*C", lbl)
        lbl <- gsub("Lake_level", "Lake-level", lbl)
        lbl <- gsub("TraCE-21Ka", "TraCE*'-'*21*Ka", lbl)
        lbl <- gsub("n-alkane", "italic(n)*'-alkane'", lbl)
        
        parts <- strsplit(lbl, "\n")[[1]]
        if(length(parts) == 2){
          if(!Add.parenthesis){paste0("atop(", parts[1], ",", parts[2], ")")} 
          else{paste0("atop(", parts[1], ",(", parts[2], "))")}} 
        else{lbl}})}
    
    if(Facet.tit.color == T){
      Color.facet <- Color.facet[match(Mmelt$Proxy[!duplicated(Mmelt$lake)], names(Color.facet))]
      
      My_facet <- ggh4x::facet_wrap2(
        vars(lake), ncol = 1, scales = "free_y", strip.position = "left", labeller = as_labeller(lake_labeller, default = label_parsed),
        strip = ggh4x::strip_themed(text_y = ggh4x::elem_list_text(face = "bold", color = Color.facet, angle = Ang, vjust = .5, hjust = .5)))
    }
    else{
      My_facet <- facet_wrap(vars(lake), ncol = 1, scales = "free_y", strip.position = "left", labeller = as_labeller(lake_labeller, default = label_parsed)) 
    }
  }
  else{
    if(Facet.tit.color == T){
      Color.facet <- Color.facet[match(Mmelt$Proxy[!duplicated(Mmelt$lake)], names(Color.facet))]
      My_facet <- ggh4x::facet_wrap2(vars(lake), ncol = 1, scales = "free_y", strip.position = "left",
                                     strip = ggh4x::strip_themed(text_y = ggh4x::elem_list_text(color = Color.facet, angle = Ang, hjust = 0.5, vjust = 5)))}
    else{My_facet <- facet_wrap(vars(lake), ncol = 1, scales = "free_y", strip.position = "left")}}
  
  #### CRUTS v4 plot ####
  if(is.null(CRUTS) == F){
    CRUTS.i <- CRUTS[c(Age.depth, Param.clim)]
    names(CRUTS.i) <- c("Age.depth", "value")
    
    All.panel <- levels(Mmelt$lake)
    if(is.null(CRUTS.display.panel) == F){All.panel <- All.panel[CRUTS.display.panel]}
    
    CRUTS.i.j <- setNames(data.frame(NA, NA, NA), c("Age.depth", "value", "lake"))
    for(i in 1:length(All.panel)){CRUTS.i.j <- rbind(CRUTS.i.j, cbind(CRUTS.i, lake = All.panel[i]))}
    CRUTS.i.j <- CRUTS.i.j[rowSums(is.na(CRUTS.i.j)) < ncol(CRUTS.i.j), ]
    CRUTS.i.j$lake <- as.factor(CRUTS.i.j$lake)
    CRUTS.line <- geom_line(inherit.aes = F, data = CRUTS.i.j, mapping = aes(x = Age.depth, y = value), color = "grey30", linetype = "solid", linewidth = 0.3, alpha = 0.6)
  }
  else{CRUTS.line <- NULL}
  
  #### Bars for Lake Level ####
  if("Lake_level" %in% unique(Mmelt$Proxy)){
    Mmelt.bar <- Mmelt
    Mmelt.bar$value[Mmelt$Proxy != "Lake_level"] <- NA
    Mmelt$value[Mmelt$Proxy == "Lake_level"] <- NA
    My_bar <- geom_segment(inherit.aes = F, data = Mmelt.bar, aes(yend = 0, y = value+0.1, xend = Age.depth, x = Age.depth, color = Proxy), alpha = .5, linewidth = Bar.size)
  }
  else{My_bar <-  NULL}
  
  #### Plot ggplot ####
  Ptot <- ggplot(Mmelt, aes(x = Age.depth, y = value, color = Proxy, shape = Shape))+ 
    Zone.clim.to.plot +
    Rect_scale +
    Event.clim +
    new_scale_fill()+
    # scale_x_continuous(name = Title.age, breaks = round(seq(Y.zero, Limites[2], by = Select.interv)), expand = c(0.02, 0.02))+
    scale_x_continuous(name = Title.age, breaks = round(seq(Y.zero, Limites[2], by = Select.interv)), expand = c(0.02, 0.02), limits = Limites)+
    CI + My_bar + My_facet + My_lines + My_color_scale+ My_fill_scale+ Surf_lines+ CRUTS.line +
    ggh4x::facetted_pos_scales(y = Axis.matrix) +
    ylab("")+
    xlab(Age.depth)+
    geom_point()+ Title.1 +
    My_smooth +
    coord_cartesian(clip = "off") +
    scale_shape_manual(values = 16, guide = "none")+
    theme(
      strip.placement = "outside", legend.position = Leg.pos,
      panel.spacing.y = unit(Panel.space.y, "mm"),
      axis.text.x = Axe.x.param,
      axis.line.y = element_line(),
      axis.line.x = Axe.x.line,
      axis.title.x = Axe.x.title,
      axis.ticks.x = Axe.x.ticks,
      strip.text.y.left = element_text(angle = Ang, size = Size.panel.tit, hjust = yhjust),
      panel.grid = element_blank(),
      panel.background = element_blank(), plot.background = element_blank(), legend.background = element_blank(),
      plot.margin = unit(c(0,0,0,0), 'pt'),
      strip.background = element_blank()
    )
  
  #### Plot Temp zones ####
  if(is.null(Name.zone) == F){
    M.annot <- Mmelt[c(which.min(Mmelt$Age.depth), which.max(Mmelt$Age.depth)),]
    M.annot$value <- 0
    
    Plot.annot <- ggplot(M.annot, aes(y = value, x = Age.depth))+
      geom_point(colour = "white", alpha = 0.01)+
      scale_x_continuous(name = Title.age, breaks = round(seq(0, Limites[2], by = Select.interv)), expand = c(0.02, 0.02))+
      Title.2 + Zone.txt1 + Rect_scale2 +
      new_scale_color()+
      scale_fill_manual(values = Rect.color.scale, guide = "none")+
      scale_color_manual(values = Rect.color.scale, guide = "none", name = NULL, labels = NULL, breaks = NULL, na.translate = FALSE)+
      theme(legend.position = "none", plot.background = element_blank(), panel.background = element_blank(),
            strip.background = element_blank(), panel.grid = element_blank(), legend.background = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = unit(c(0,0,0,0), 'pt'),
            axis.text = element_blank()
      )
    
    Ptot <- Plot.annot/Ptot +  plot_layout(heights = c(2, 98)) & theme(panel.grid = element_blank(), legend.background = element_blank())
    # Ptot <- Plot.annot/Ptot + plot_layout(heights = c(2, 98, 98))
    # Ptot <- Plot.annot
    
  }
  
  #### Save plot ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(Ptot, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.plot)}}
  return(Ptot)
  
}

Map.panel.climat <- function(List.df = NULL, Timespan = NULL, Period.lab = NULL, Param.clim = c("MAAT", "MAP"), Meta.data = NULL, Emphase.proxy = NULL, Elevation.range = NULL,
                             Color.by.rank = T, Add.emptysite = T, Map.dot.size = 2.7, Title.size = 12, Dot.alpha = 1, Display.plot = T, Use.Meta.data.4.empty.site = F, Show.stats = F,
                             Empty.site.shape = 16, Add.DEM = F, DEM.alpha = .8, Remove.overlap = F, Overlap.seuil = 60000, Show.label = T, Altitudinal.gradient = F,
                             Map.leg.pos = c(0.002, 0.999), Map.limits = c(42, 123, 25, 63.5), Save.plot = NULL, H = 850, W = 800) {
  
  #### Select coordinates of the data ####
  Coords_used <- Meta.data[Meta.data$ID %in% names(List.df), c("ID", "Latitude", "Longitude", "Altitude")]
  Coords_used$Alt_rank <- rank(Coords_used$Altitude, ties.method = "first")
  
  #### Overlapping dots ####
  if(Remove.overlap == T){
    library(sf)
    Sites_df_plot <- Coords_used
    Sites_df_plot <- Sites_df_plot[!is.na(Sites_df_plot$Latitude), ]
    Sites_df_plot <- Sites_df_plot[!is.na(Sites_df_plot$Longitude), ]
    
    sites_sf <- st_as_sf(
      Sites_df_plot,
      coords = c("Longitude", "Latitude"),
      crs = 4326
    )
    
    sites_proj <- st_transform(sites_sf, 3857)  # Web Mercator (meters)
    coords <- st_coordinates(sites_proj)
    
    dot_radius <- Overlap.seuil   # 200 km radius for big dots
    max_iter   <- 999       # max iterations for overlap resolution
    seed <- 123
    
    sites_sf <- st_as_sf(Sites_df_plot, coords = c("Longitude", "Latitude"), crs = 4326)
    sites_proj <- st_transform(sites_sf, 3857)
    coords <- st_coordinates(sites_proj)
    coords_new <- coords
    set.seed(seed)
    
    for (iter in 1:max_iter) {
      moved <- FALSE
      n <- nrow(coords_new)
      for (a in 1:(n-1)) {
        for (b in (a+1):n) {
          dx <- coords_new[b,1] - coords_new[a,1]
          dy <- coords_new[b,2] - coords_new[a,2]
          dist <- sqrt(dx^2 + dy^2)
          
          # Handle exact same coordinates
          if (dist == 0) {
            angle <- runif(1, 0, 2*pi)
            dx <- cos(angle) * 1e-6
            dy <- sin(angle) * 1e-6
            dist <- sqrt(dx^2 + dy^2)
          }
          
          overlap <- 2*dot_radius - dist
          if (!is.na(overlap) && overlap > 0) {
            # push both points apart
            shift_x <- (dx / dist) * (overlap / 2)
            shift_y <- (dy / dist) * (overlap / 2)
            coords_new[a,1] <- coords_new[a,1] - shift_x
            coords_new[a,2] <- coords_new[a,2] - shift_y
            coords_new[b,1] <- coords_new[b,1] + shift_x
            coords_new[b,2] <- coords_new[b,2] + shift_y
            moved <- TRUE
          }
        }
      }
      if (!moved) break
    }
    
    sites_proj$geometry <- st_sfc(lapply(1:nrow(coords_new), function(i) {st_point(coords_new[i, ])}), crs = st_crs(sites_proj))
    
    sites_proj$x <- coords_new[,1]
    sites_proj$y <- coords_new[,2]
    sites_latlon <- st_transform(sites_proj, 4326)
    coords_latlon <- st_coordinates(sites_latlon)
    sites_latlon$Longitude <- coords_latlon[,1]
    sites_latlon$Latitude  <- coords_latlon[,2]
    Coords_used <- sites_latlon[c("ID", "Latitude", "Longitude", "Altitude", "Alt_rank")]
  }
  
  #### Bin calculation and ranking ####
  List.save <- list()
  if(is.null(Period.lab) == T){Period.lab <- paste("P", Timespan[-length(Timespan)], sep = "")}
  
  for(i in Param.clim){
    #### Bin by period ####
    Clim.bin <- function(df, var, Timespan, Period.lab){
      df <- df[df$Age >= Timespan[1],]
      df <- df[df$Age <= Timespan[length(Timespan)],]
      df$Groups <- base::cut(x = df$Age, breaks = Timespan, labels = Period.lab, right = F)
      df_summary <- df %>%
        group_by(Groups) %>%
        summarise(Mean = mean(.data[[var]], na.rm = TRUE)) %>%
        tidyr::pivot_wider(names_from = Groups, values_from = Mean)
      df_summary <- as.data.frame(df_summary)
      return(df_summary)
    }
    Tab.clim <- lapply(List.df, Clim.bin, var = i, Timespan = Timespan, Period.lab = Period.lab)
    Tab.clim <- lapply(names(Tab.clim), function(nm) {df <- Tab.clim[[nm]]; df$Site <- nm; return(df)})
    Tab.clim <- bind_rows(Tab.clim)
    Tab.clim <- Tab.clim[which(names(Tab.clim) != "NA")]
    
    #### Rank of param.clim by period ####
    rank_labels <- paste("V", seq(1:length(Period.lab)), sep = "")
    rank_periods <- function(df, period_labels, rank_labels) {
      df %>%
        rowwise() %>%
        mutate(across(all_of(period_labels), ~ {
          vals <- c_across(all_of(period_labels))
          if (any(is.na(vals))) {
            if (is.na(cur_data()[[cur_column()]])) {
              return(NA_integer_)
            }
          }
          rank(vals, ties.method = "first")[which(names(cur_data()) == cur_column())]
        }, .names = "{.col}_order")) %>%
        ungroup() %>%
        mutate(across(ends_with("_order"), ~ ifelse(!is.na(.), rank_labels[.], NA_character_)))
    }
    if(Color.by.rank == T){
      Tab.clim <- rank_periods(Tab.clim, Period.lab, rank_labels)
      Tab.clim <- Tab.clim[!names(Tab.clim) %in% Period.lab]
      Tab.clim <- reshape2::melt(Tab.clim, id = names(Tab.clim)[!names(Tab.clim) %in% paste(Period.lab, "order", sep = "_")])
      Tab.clim$variable <- gsub("_order", "", Tab.clim$variable) 
    }
    else{
      Tab.clim <- melt(Tab.clim, id = names(Tab.clim)[!names(Tab.clim) %in% Period.lab])
    }
    
    Tab.clim <- cbind(Tab.clim, Coords_used[match(Tab.clim$Site, Coords_used$ID), c("Latitude", "Longitude", "Altitude", "Alt_rank")])
    Tab.clim$Lab.clim <- i
    row.names(Tab.clim) <- NULL
    List.save[[which(Param.clim %in% i)]] <- Tab.clim
  }
  Tab.clim <- bind_rows(List.save)
  Tab.clim$variable <- factor(Tab.clim$variable, levels = Period.lab)
  
  #### Graphical param ####
  if(Display.plot == T & exists("ACA.bo") == F){source("Scripts/Load_maps.R")}
  
  if(Color.by.rank == T){
    Tab.clim$value[which(Tab.clim$Lab.clim == "MAAT" & Tab.clim$value == "V1")] <- "Cold"
    Tab.clim$value[which(Tab.clim$Lab.clim == "MAAT" & Tab.clim$value == "V2")] <- "Mild"
    Tab.clim$value[which(Tab.clim$Lab.clim == "MAAT" & Tab.clim$value == "V3")] <- "Warm"
    Tab.clim$value[which(Tab.clim$Lab.clim == "MAP" & Tab.clim$value == "V1")] <- "Arid"
    Tab.clim$value[which(Tab.clim$Lab.clim == "MAP" & Tab.clim$value == "V2")] <- "Semi-arid"
    Tab.clim$value[which(Tab.clim$Lab.clim == "MAP" & Tab.clim$value == "V3")] <- "Sub-humid"
    
    my_col_group <- c("Cold" = "#2c7bb6", "Mild" = "#ffffbf", "Warm" = "#d7191c",
                      "Arid" = "#89520B", "Semi-arid" = "#DFC27E", "Sub-humid" = "#80CDC1")
    # Col_scale <- scale_color_manual(values = my_col_group, na.translate = F)
    Fill_scale.MAP <- scale_fill_manual(values = my_col_group, na.translate = F)
    Fill_scale.MAAT <- Fill_scale.MAP
  }
  else{
    Fill_scale.MAP <- scale_fill_gradientn(colors = c("#89520B", "#DFC27E", "#80CDC1"))
    Fill_scale.MAAT <- scale_fill_gradientn(colors = c("#2c7bb6", "#ffffbf", "#d7191c"))
  }
  
  My_labels <- labeller(variable = c(
    "Late Hol." = "Late Holocene",
    "Mid-Hol."  = "Middle Holocene",
    "Early Hol."= "Early Holocene"),
    Lab.clim = c("MAAT" = "(A) Temperature", "MAP" = "(B) Precipitation"))
  
  label_longitude <- function(x) {paste0(abs(x), "°", ifelse(x < 0, "W", "E"))}
  label_latitude <- function(y) {paste0(abs(y), "°", ifelse(y < 0, "S", "N"))}
  
  Tab.clim.MAAT <- Tab.clim[Tab.clim$Lab.clim == "MAAT",]
  Tab.clim.MAP <- Tab.clim[Tab.clim$Lab.clim == "MAP",]
  
  #### Display DEM ####
  if(Display.plot == T & Add.DEM == T){
    Path.world.DEM = paste(Path.SIG, "SIG/MNT/WorldClim2_MNT/wc2.1_10m_elev.tif", sep = "")
    DEM = raster::raster(Path.world.DEM)
    # browser()
    DEM <- raster::crop(DEM, raster::extent(Map.limits[1]-3, Map.limits[2]+3, Map.limits[3]-1, Map.limits[4]+1))
    DEM_lowres <- raster::aggregate(DEM, fact = 1, fun = mean, na.rm = TRUE)
    
    breaks <- c(-Inf, 0, 250, 500, 1000, 2000, 3000, 4000, 5000, Inf)
    
    DEM_disc <- raster::cut(
      DEM_lowres,
      breaks = breaks,
      include.lowest = TRUE
    )
    
    DEM_df <- raster::as.data.frame(DEM_disc, xy = TRUE, na.rm = TRUE)
    colnames(DEM_df) <- c("x", "y", "elev_class")
    
    My_DEM <-  geom_raster(data = DEM_df, aes(x = x, y = y, fill = elev_class), alpha = DEM.alpha)
    DEM_scale <- scale_fill_gradientn(colors = c(NA, "grey90", "grey80", "grey60", "tan3", "sienna4", "white"), guide = "none")
  }
  else{My_DEM <- NULL; DEM_scale <- NULL}
  
  #### Empty sites ####
  if(Add.emptysite == T){
    if(Use.Meta.data.4.empty.site == F){df.empty <- Meta.data.ACA}
    else{df.empty <- Meta.data}
    
    df.empty <- df.empty[!grepl("\\|TRA$", df.empty$ID),]
    df.empty <- df.empty[!df.empty$ID %in% unique(Tab.clim.MAAT$Site),]
    # browser()
    Empty.site <- geom_point(data = df.empty, mapping = aes(x = Longitude, y = Latitude),  alpha = 1, size = Map.dot.size*.2, shape = Empty.site.shape, color = "grey10")
  }
  else{Empty.site <- NULL}
  
  #### Emphase shape ####
  if(is.null(Emphase.proxy) == F){
    My_dots_MAAT <- geom_point(data = Tab.clim.MAAT, mapping = aes(x = Longitude, y = Latitude,  fill = value),  alpha = Dot.alpha, size = Map.dot.size, shape = 21, color = "black", na.rm = T)
    My_dots_MAP <- geom_point(data = Tab.clim.MAP, mapping = aes(x = Longitude, y = Latitude,  fill = value),  alpha = Dot.alpha, size = Map.dot.size, shape = 21, color = "black", na.rm = T)
  }
  else{
    My_dots_MAAT <- geom_point(data = Tab.clim.MAAT, mapping = aes(x = Longitude, y = Latitude,  fill = value),  alpha = Dot.alpha, size = Map.dot.size, shape = 21, color = "black", na.rm = T)
    # My_dots_MAAT <- geom_hex(bins = 60, data = Tab.clim.MAAT, mapping = aes(x = Longitude, y = Latitude,  fill = value, color = value, group = Site))
    My_dots_MAP <- geom_point(data = Tab.clim.MAP, mapping = aes(x = Longitude, y = Latitude,  fill = value),  alpha = Dot.alpha, size = Map.dot.size, shape = 21, color = "black", na.rm = T)
    # My_dots_MAP <- geom_label(data = Tab.clim.MAP, mapping = aes(x = Longitude, y = Latitude, label = Site), size = 1, color = "black", na.rm = T)
  }
  
  #### Show labels ####
  if(Show.label == T){Site.labs <- geom_text_repel(data = Coords_used, aes(label = ID, x = Longitude, y = Latitude), direction = "y", 
                                                   hjust = -0.5, #vjust = 2,
                                                   # force = 100, min.segment.length = 3000,
                                                   # na.rm = T, xlim = max(merged_df$value),
                                                   segment.square = F, segment.shape = 1, segment.linetype = 1, segment.size = 0.1, segment.colour = "black",
                                                   segment.curvature = 0.2, segment.inflect = F, min.segment.length = 0.5,
                                                   size = 2, parse = F)}
  else{Site.labs <- NULL}
  
  #### Display stats ####
  if(Show.stats == T){
    Nb.site <- length(unique(Tab.clim.MAAT$Site))
    
    Stats.MAAT <- Tab.clim.MAAT[c("variable", "value")]
    Stats.MAAT <- group_by(Stats.MAAT, variable)
    Stats.MAAT <- dplyr::count(Stats.MAAT, value)
    Stats.MAAT$perc <- Stats.MAAT$n/Nb.site*100
    
    Stats.MAP <- Tab.clim.MAP[c("variable", "value")]
    Stats.MAP <- group_by(Stats.MAP, variable)
    Stats.MAP <- dplyr::count(Stats.MAP, value)
    Stats.MAP$perc <- Stats.MAP$n/Nb.site*100
    
    print("MAAT")
    print(Stats.MAAT)
    print("MAP")
    print(Stats.MAP)
  }
  
  #### Plot maps ####
  if(Display.plot == T){
    pmap.MAAT <- map_empty_world +
      new_scale_color() + new_scale_fill() +
      My_DEM + DEM_scale +
      new_scale_color() + new_scale_fill() +
      geom_polygon(data = World_map, aes(x=long, y=lat, group = group), alpha = 1, fill = NA, color = "grey30", size = 0.3)+
      geom_polygon(data = ACA.bo.proj, aes(x=long, y=lat), colour="black", alpha = 1, fill = NA, size = 0.5)+
      coord_quickmap(xlim = c(Map.limits[1], Map.limits[2]), ylim = c(Map.limits[3], Map.limits[4])) +
      Empty.site +
      new_scale_fill()+ Fill_scale.MAAT + 
      My_dots_MAAT + Site.labs +
      facet_grid(variable ~ Lab.clim, labeller = My_labels)+
      guides(fill = guide_legend(nrow = 1))+
      theme(strip.text.y = element_blank())
    
    pmap.MAP <- map_empty_world +
      new_scale_color() + new_scale_fill() +
      My_DEM + DEM_scale +
      new_scale_color() + new_scale_fill() +
      geom_polygon(data = World_map, aes(x=long, y=lat, group = group), alpha = 1, fill = NA, color = "grey30", size = 0.3)+
      geom_polygon(data = ACA.bo.proj, aes(x=long, y=lat), colour="black", alpha = 1, fill = NA, size = 0.5)+
      coord_quickmap(xlim = c(Map.limits[1], Map.limits[2]), ylim = c(Map.limits[3], Map.limits[4])) +
      Empty.site +
      new_scale_fill()+ Fill_scale.MAP + 
      My_dots_MAP + Site.labs +
      facet_grid(variable ~ Lab.clim, labeller = My_labels)+
      guides(fill = guide_legend(nrow = 1))+
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    p <- pmap.MAAT + pmap.MAP & theme(strip.background = element_blank(), legend.title = element_blank(), strip.text = element_text(size = Title.size, face = "bold")) &
      theme(
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "#7BACCF"),
        panel.border = element_rect(colour = "black", linewidth = 1),
        legend.position = Map.leg.pos,
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 1), color = "black"),
        legend.box.background = element_blank()
      ) &
      scale_x_continuous(labels = label_longitude) &
      scale_y_continuous(labels = label_latitude)
  }
  else{p <- NULL}
  
  #### Plot altitudinal gradient ####
  if(Altitudinal.gradient == T){
    Show.alt.boxplot = T
    Add.Tukey = F
    My_width <- c(4/5, 1/5)
    if(is.null(Elevation.range) == T){Elevation.range <- c(round(min(Tab.clim.MAAT$Altitude, na.rm = T), digits = 0), round(max(Tab.clim.MAAT$Altitude, na.rm = T), digits = 0))}
    
    Tab.clim.MAAT <- Tab.clim.MAAT[!is.na(Tab.clim.MAAT$Altitude),]
    Tab.clim.MAP <- Tab.clim.MAP[!is.na(Tab.clim.MAP$Altitude),]
    
    
    #### Compact letters ####
    anova_tukey_letters_by_panel <- function(data, response, group, panel, alpha = 0.05) {
      library(multcompView)
      
      #### Data cleaning ####
      data[[group]] <- as.factor(data[[group]])
      data[[panel]] <- as.factor(data[[panel]])
      
      panels <- levels(data[[panel]])
      out <- list()
      fml <- as.formula(paste(response, "~", group))
      data[[group]] <- gsub("-", "_", data[[group]])
      
      #### Loop on time period ####
      for (p in panels) {
        sub <- data[data[[panel]] == p, ]
        ng <- length(unique(sub[[group]]))
        
        if (ng < 2) {
          out[[p]] <- data.frame(
            panel  = p,
            group  = unique(sub[[group]]),
            letters = NA_character_, row.names = NULL)
          next
        }
        
        print(fml)
        model <- aov(fml, data = sub)
        print(summary(model))
        tukey <- TukeyHSD(model, conf.level = 1 - alpha)
        # print(tukey)
        cld   <- multcompView::multcompLetters4(model, tukey)
        
        letters <- cld[[1]]$Letters
        
        out[[p]] <- data.frame(panel  = p, group  = names(letters), letters = letters, row.names = NULL)
      }
      
      #### Merge all time period ####
      res <- do.call(rbind, out)
      colnames(res)[1:2] <- c(panel, group)
      
      res[[group]] <- gsub("_", "-", res[[group]])
      data[[group]] <- gsub("_", "-", data[[group]])
      
      #### ggplot Letters ####
      y_pos <- aggregate(
        Altitude ~ variable + value,
        data = data,
        FUN = function(x) max(x, na.rm = T)# * 1.05
        # FUN = function(x) {
        #   q3  <- quantile(x, 0.75, na.rm = TRUE)
        #   iqr <- IQR(x, na.rm = TRUE)
        #   iqr <- q3 + iqr   # 5% of IQR above Q3
        #   return(iqr)
        #   }
      )
      # res <- merge(res, y_pos, by = c("variable", "value"), all.x = T)
      res$y_pos <- max(data[[response]]) + 0.07*max(data[[response]])
      
      my_letters <- geom_text(data = res, aes(y = y_pos, label = letters), size = 3)
      return(my_letters)
    }
    
    if(Add.Tukey == T){
      df.tuk.MAAT <- Tab.clim.MAAT[c("Site", "variable", "value", "Altitude")]
      df.tuk.MAP <- Tab.clim.MAP[c("Site", "variable", "value", "Altitude")]
      
      letters_MAAT <- anova_tukey_letters_by_panel(data = df.tuk.MAAT, response = "Altitude", group = "value", panel = "variable")
      letters_MAP <- anova_tukey_letters_by_panel(data = df.tuk.MAP, response = "Altitude", group = "value", panel = "variable")
    }
    else{letters_MAP <- NULL; letters_MAAT <- NULL}
    
    #### MAAT ####
    palt.MAAT <- ggplot(data = Tab.clim.MAAT, mapping = aes(x = Alt_rank, y = Altitude)) +
      # palt.MAAT <- ggplot(data = Tab.clim.MAAT, mapping = aes(x = Longitude, y = Altitude)) +
      geom_area(fill = "grey90", color = "grey50", size = .2)+
      geom_point(aes(fill = value),  alpha = Dot.alpha, size = Map.dot.size, shape = 21, color = "black", na.rm = T)+
      # coord_quickmap(xlim = c(Map.limits[1], Map.limits[2]), ylim = c(Map.limits[3], Map.limits[4])) +
      # Empty.site +
      Fill_scale.MAAT + ylim(Elevation.range) +
      Site.labs + xlab("Elevation ranking (#)")+ ylab("Elevation (m a.s.l.)")+
      facet_grid(variable ~ Lab.clim, labeller = My_labels)+
      guides(fill = guide_legend(nrow = 1))+
      theme(strip.text.y = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            plot.margin = unit(c(0,0,0,0), "lines"),
            legend.position = Map.leg.pos, legend.title = element_blank())
    
    #### MAAT boxplots ####
    if(Show.alt.boxplot == T){
      palt.MAAT.dist <- ggplot(data = Tab.clim.MAAT, mapping = aes(x = value, y = Altitude, fill = value)) +
        Fill_scale.MAAT +
        geom_violin(alpha = .8, fill = "grey90", na.rm = F, show.legend = F, trim = F, scale = "width", color = NA)+
        geom_boxplot(outlier.colour = "red", outlier.shape = NA, alpha = 0.6, na.rm = T, notch = F)+
        letters_MAAT + 
        facet_grid(variable ~ Lab.clim, labeller = My_labels)+
        guides(fill = guide_legend(nrow = 1))+ ylim(Elevation.range) +
        theme(strip.text = element_blank(), axis.title = element_blank(), 
              plot.margin = unit(c(0,0,0,0), "lines"),
              axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
      
      palt.MAAT <- wrap_plots(palt.MAAT, palt.MAAT.dist, widths = My_width)
    }
    
    #### MAP ####
    palt.MAP <- ggplot(data = Tab.clim.MAP, mapping = aes(x = Alt_rank, y = Altitude)) +
      # palt.MAP <- ggplot(data = Tab.clim.MAP, mapping = aes(x = Longitude, y = Altitude)) +
      geom_area(fill = "grey90", color = "grey50", size = .2)+
      geom_point(aes(fill = value),  alpha = Dot.alpha, size = Map.dot.size, shape = 21, color = "black", na.rm = T)+
      # coord_quickmap(xlim = c(Map.limits[1], Map.limits[2]), ylim = c(Map.limits[3], Map.limits[4])) +
      Fill_scale.MAP + ylim(Elevation.range) +
      Site.labs + xlab("Elevation ranking (#)")+ ylab("Elevation (m a.s.l.)")+
      facet_grid(variable ~ Lab.clim, labeller = My_labels)+
      guides(fill = guide_legend(nrow = 1))+
      theme(strip.text.y = element_blank(), legend.position = Map.leg.pos, plot.margin = unit(c(0,0,0,0), "lines"),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
            axis.title.y = element_blank(), legend.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
    
    #### MAP boxplots ####
    if(Show.alt.boxplot == T){
      palt.MAP.dist <- ggplot(data = Tab.clim.MAP, mapping = aes(x = value, y = Altitude, fill = value)) +
        Fill_scale.MAP +
        geom_violin(alpha = .8, fill = "grey90", na.rm = F, show.legend = F, trim = F, scale = "width", color = NA)+
        geom_boxplot(outlier.colour = "red", outlier.shape = NA, alpha = 0.6, na.rm = T, notch = F)+
        letters_MAP +
        facet_grid(variable ~ Lab.clim, labeller = My_labels)+
        guides(fill = guide_legend(nrow = 1))+ ylim(Elevation.range) +
        theme(strip.text.x = element_blank(), axis.title = element_blank(), 
              plot.margin = unit(c(0,0,0,0), "lines"),
              axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")
      
      palt.MAP <- wrap_plots(palt.MAP, palt.MAP.dist, widths = My_width)
    }
    
    
    p <- (palt.MAAT | palt.MAP) &
      theme(
        plot.background = element_blank(), strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 1), color = "black"),
        legend.box.background = element_blank()
      ) #&
    # scale_x_continuous(labels = label_longitude) &
  }
  
  #### Save + export ####
  if(is.null(Save.plot) == F & is.null(p) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(p, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", limitsize = F)}
    else{ggsave(Save.plot)}}
  
  return(Tab.clim)
}

Plot.stack.groups <- function(List.clim = NULL, Meta.data = NULL, Group = NULL, Add.vlines = NULL, Facet.by.group = F, Subgroup = NULL, Color.by.clim = F,
                              Add.count.per.bin = F, Return.plot = F, Plot.2.merge = F, Show.strip = T, Display.map = F, Add.smooth = F, Add.emptysite = T, Width.pour.map = NULL,
                              Map.dot.size = 2.7, Map.leg.pos = c(0.003, 0.997), Map.limits = c(40, 130, 25, 70), Bar.size = 1.5, Bar.alpha = .3, Smooth.param = .4, CI.alpha = .4,
                              Title.size = 14, Select.interv = 1000, Title.age = "Time (cal. year BP)", Limits = NULL, Remove.clim = NULL, Param.clim = c("MAAT", "MAP"),
                              Time.lab.rot = 45, Add.error.zone = F, Error.alpha = .3, Empty.site.from.Meta.data = F, Remove.unused.facet = F,
                              Scaling = T, Windows.length = 400, Anomaly = F, Display.CI = T, Leg.ncol = 1, CI.perc = 95, Facet.order = NULL, Remove.extrems = F, Name.param.clim = NULL,
                              H = NULL, W = NULL, Save.plot = NULL, Save.plot.Rds = NULL, Save.data.Rds = NULL){
  #### Groups cores ####
  if(any(names(Meta.data) %in% Group) == F){Group <- NULL}
  if(any(names(Meta.data) %in% Subgroup) == F){Subgroup <- NULL}
  
  if(is.null(Group) == F){
    Meta.data <- Meta.data[Meta.data$ID %in% names(List.clim),]
    Meta.data.gp <- split(Meta.data, Meta.data[Group])
    if(is.null(Subgroup) == T){
      print(paste("**** Stack.plot by ", Group, ". ****", sep = ""))
      Site.gp <- lapply(Meta.data.gp, function(x) unlist(x["ID"]))
    }
    else{
      print(paste("**** Stack.plot by ", Group, " and ", Subgroup, ". ****", sep = ""))
      My_subgroups <- unique(Meta.data[[Subgroup]])
      Meta.data <- Meta.data[Meta.data$ID %in% names(List.clim),]
      Meta.data.gp <- split(Meta.data, interaction(Meta.data[[Group]], Meta.data[[Subgroup]], sep = "|"))
      Site.gp <- lapply(Meta.data.gp, function(x) unlist(x["ID"]))
      Site.gp <- Filter(function(x) length(x) != 0, Site.gp)
      
    }
    Grouped.List <- lapply(Site.gp, function(names_vec) {List.clim[names_vec]})
    Res <- lapply(Grouped.List, function(x){Stacking.quantif(x, Binning = T, Scaling = Scaling, Add.bin.count = T, Limits = Limits, Anomaly = Anomaly, Windows.length = Windows.length, 
                                                             Stacking = T, Keep.empty.bin = F, Bin.sd = T, Plot.x = "Age", Keep.clim = Param.clim)})
    
    
    # if(Scaling == T){Res <- lapply(Res, function(x) cbind(Age = x$Age, data.frame(apply(x[c(2,3)], 2, Scaling.fun)), SD_MAAT = x$SD_MAAT, SD_MAP = x$SD_MAP, N = x$N))}
    # Scaling = F
    
    Res <- lapply(names(Res), function(name) {df <- Res[[name]]; df$Group <- name; return(df)})
  }
  else{Res <- Stacking.quantif(List.clim, Binning = T, Add.bin.count = T, Scaling = Scaling, Limits = Limits, Windows.length = Windows.length, Anomaly = Anomaly, Keep.empty.bin = F, Bin.sd = T, Plot.x = "Age", Keep.clim = Param.clim)}
  Res <- bind_rows(Res)
  
  if(is.null(Subgroup) == F){
    Res$Subgroup <- gsub(".*\\|", "", Res$Group)
    Res$Group <- gsub("\\|.*", "", Res$Group)
  }
  
  Res.save <- Res
  if(Remove.extrems == T){
    Res <- Res[!seq(1:nrow(Res)) %in% c(which.min(Res$Age), which.max(Res$Age)),]
    Res <- Res[!seq(1:nrow(Res)) %in% c(which.min(Res$Age), which.max(Res$Age)),]}
  
  if(is.null(Remove.clim) == F){Res <- Res[!grepl(Remove.clim, names(Res))]}
  
  #### Melt + Confidence interval ####
  Keep.unmelt <- c("Age", "N")
  if(is.null(Group) == F){Keep.unmelt <- c(Keep.unmelt, "Group")}
  if(is.null(Subgroup) == F){Keep.unmelt <- c(Keep.unmelt, "Subgroup")}
  if(any(grepl("^SD_", names(Res))) == T){
    df_means <- reshape2::melt(Res, id.vars = Keep.unmelt, measure.vars = names(Res)[!names(Res) %in% Keep.unmelt & !grepl("^SD_", names(Res))], variable.name = c("Param.clim"), value.name = "Val.clim")
    df_sds <- reshape2::melt(Res, id.vars = Keep.unmelt, measure.vars = names(Res)[grepl("^SD_", names(Res))], variable.name = "Param.clim", value.name = "Val.SD")
    df_sds$Param.clim <- gsub("SD_", "", df_sds$Param.clim)
    Res <- left_join(df_means, df_sds, by = c(Keep.unmelt, "Param.clim"))
    
    # Res$se <- Res$Val.SD #/ sqrt(Res$N)  # Standard error
    Res$se <- Res$Val.SD / sqrt(Res$N)  # Standard error
    # Res$t_crit <- qt(0.975, Res$N - 1) # Calculate t critical value for each row
    Res$t_crit <- qt((1 + CI.perc/100)/2, pmax(Res$N - 1, 1))  # avoid df=0
    # Res$t_crit <- ifelse(Res$N <= 3, 3, Res$t_crit)
    # Res$note <- ifelse(Res$N <= 3, "N<3: CI unstable", "")
    Res$CI_lower <- Res$Val.clim - Res$t_crit * Res$se
    Res$CI_upper <- Res$Val.clim + Res$t_crit * Res$se
    # Res$CI_lower <- Res$Val.clim - 1.86*Res$se# /sqrt(Res$N)
    # Res$CI_upper <- Res$Val.clim + 1.86*Res$se# /sqrt(Res$N)
    CI <- geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), color = NA, size = 0.15, linetype = "dashed", alpha = CI.alpha)
  }
  else{
    Res <- melt(Res, id = Keep.unmelt)
    names(Res)[c(ncol(Res)-1,ncol(Res))] <- c("Param.clim", "Val.clim")
    CI <- NULL
  }
  
  #### Graphical settings ####
  if(Add.count.per.bin == T){Tit.x <- element_blank(); Ticks.x <- element_blank(); Lab.x <- element_blank()}
  else{Tit.x <- element_text(); Ticks.x <- element_line(); Lab.x <- element_text(angle = Time.lab.rot, hjust = 1, vjust = 1)}
  
  if(Plot.2.merge == T){
    Tit.x <- element_blank(); Lab.x <- element_blank()
    Tit.x2 <- element_blank(); Lab.x2 <- element_blank()
  }
  else{Tit.x2 <- element_text(); Lab.x2 <- element_text()}
  
  Strip.x <- element_text(size = Title.size, face = "bold")
  if(Show.strip == T & Display.map == F){
    Strip.y <- element_text(size = Title.size, face = "bold")
    Strip.map <- element_blank()}
  else{
    # if(Show.strip == F){
    # Strip.x <- element_blank()
    Strip.y <- element_blank()
    Strip.map <- element_text(size = Title.size, face = "bold")
  }
  
  if(Display.map == T){Leg.pos <- "none"}
  else{Leg.pos <- "bottom"}
  
  My_glob_theme <- theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid = element_blank(), legend.key = element_blank(),
                         strip.background = element_blank(), plot.margin = unit(c(0,0,0,0), "pt"))
  
  if(Add.smooth == T){My_smooth <- geom_smooth(method = "loess", se = T, alpha = .05, fullrange = F, span = Smooth.param, na.rm = T, linewidth = .8, formula = 'y ~ x')}
  else{My_smooth <- NULL}
  
  if(is.null(Add.vlines) == F){
    # My_vlines <- geom_vline(xintercept = Add.vlines, linetype = "dashed", color = "grey40")
    My_vlines <- geom_vline(xintercept = Add.vlines, col = "grey30", lty = 2, alpha = 0.7)
  }
  else{My_vlines <- NULL}
  
  if(is.null(Limits)){Limits <- c(round(min(Res$Age, na.rm = T), digits = -1), round(max(Res$Age, na.rm = T), digits = -1))}
  
  #### Color settings ####
  if(is.null(Group) == F){
    Leg.tit <- paste("Grouped by", gsub("_", " ", Group))
    my_col_group <- c("1. Data" = "royalblue", "2. Model" = "#b37878ff", "Error" = "grey10",
                      "ECA" = "darkorange", "WCA" = "royalblue", 
                      "High-elevation" = "#AB5B33", "Low-elevation" = "#86AE61", "Mid-elevation" = "#F8AA00",
                      "Pollen" = "#BEA33A", "GDGT" = "#CF5612", "Speleothem" = "#47BA62", "Lake_level" = "royalblue", "n-alkane" = "#86AE61",
                      "Loess" = "#755841", "Paleodust" = "#755841", "d18O (Lake)" = "#47BA62", 
                      "Dust (Loess)" = "#755841", "Dust (Lake)" = "#755841",
                      "ACA" = "#A74639", "Monsoon" = "#22A767",
                      
                      "MAAT" = "darkred", "MAP" = "royalblue", "2. Model" = "#b37878ff", "3. Data" = "pink",
                      
                      "Iran" = "#219F49", "Caucasus" = "#F8AA00", "Caspian" = "#10B0C8", "Caspian-Aral" = "#10B0C8", "West. TSH" = "#A250A4", "East. TSH" = "#DF2E26",
                      
                      "1. Data.MAAT" = "darkred", "1. Data.MAP" = "darkblue", "2. Model.MAAT" = "#b37878ff", "2. Model.MAP" = "#646499ff",
                      "1. Data.Temperature" = "darkred", "1. Data.Moisture" = "darkblue", "2. Model.Temperature" = "#b37878ff", "2. Model.Moisture" = "#646499ff"
    )
    
    Col_scale <- scale_color_manual(values = my_col_group, name = Leg.tit)
    Fill_scale <- scale_fill_manual(values = my_col_group, name = Leg.tit)
    
  }
  else{Leg.tit <- "Full sites"; Col_scale <- NULL; Fill_scale <- NULL}
  
  #### Facet param ####
  if(is.null(Facet.order) == F){
    Res$Group <- factor(Res$Group, levels = Facet.order)
    if(Remove.unused.facet == T){Res <- Res[Res$Group %in% levels(Res$Group),]}
  }
  
  facet_map <- NULL
  if(Scaling == F){
    if(Facet.by.group == F){
      my_facet <- facet_wrap(. ~ Param.clim, scale = "free_y")}
    else{
      my_facet <- facet_wrap(Group ~ Param.clim, scale = "free_y", ncol = 2)
      # Strip.x <- element_blank()
    }}
  else{
    if(Facet.by.group == F){
      my_facet <- facet_grid(. ~ Param.clim)}
    else{
      my_facet <- facet_grid(Group ~ Param.clim)
      facet_map <- facet_grid(Group ~ .)
      if(is.null(Facet.order) == F){
        Meta.data[[Group]] <- factor(Meta.data[[Group]], levels = Facet.order)
        if(Remove.unused.facet == T){Meta.data <- Meta.data[Meta.data[[Group]] %in% Facet.order,]}
        
      }
      Map.leg.pos <- "none"}}
  
  if(Display.CI == T){CI <- CI}
  else{CI <- NULL}
  
  #### Rename param.clim ####
  if(is.null(Name.param.clim) == F){
    for(i in 1:length(unique(Res$Param.clim))){
      Res$Param.clim[Res$Param.clim == unique(Res$Param.clim)[i]] <- Name.param.clim[i]}
    Res$Param.clim <- factor(Res$Param.clim, levels = Name.param.clim)
  }
  
  #### Subgroup param ####
  if(is.null(Subgroup) == T){My_plot <- ggplot(Res, aes(x = Age, y = Val.clim, color = Group, fill = Group))}
  else{
    Res$Subgroup <- factor(Res$Subgroup, levels = rev(unique(Res$Subgroup)))
    My_plot <- ggplot(Res, aes(x = Age, y = Val.clim, color = Subgroup, fill = Subgroup))
    
    if(Color.by.clim == T){
      My_plot <- ggplot(Res, aes(x = Age, y = Val.clim, color = interaction(Subgroup, Param.clim), fill = interaction(Subgroup, Param.clim)))
    }
  }
  
  #### Add error zone ####
  if(Add.error.zone == T){
    df_wide <- Res %>%
      group_by(Group, Param.clim, Age, Subgroup) %>%
      summarise(Val.clim = mean(Val.clim), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Subgroup, values_from = Val.clim)
    
    names(df_wide)[names(df_wide) == unique(Res$Subgroup)[1]] <- "A"
    names(df_wide)[names(df_wide) == unique(Res$Subgroup)[2]] <- "B"
    df_wide$Val.clim <- 0
    df_wide$Fill <- "Error"
    Error.zn <- geom_ribbon(data = df_wide, aes(x = Age, ymin = pmin(A, B, na.rm = T), ymax = pmax(A, B, na.rm = T), fill = Fill), color = NA, alpha = Error.alpha)
  }
  else{Error.zn <- NULL}
  
  #### Plot results ####
  p <- My_plot + Error.zn +
    geom_point() + geom_line() +
    my_facet + My_vlines +
    CI + My_smooth + Fill_scale + Col_scale + ylab("z-scores")+
    scale_x_continuous(name = Title.age, breaks = round(seq(Limits[1], Limits[2], by = Select.interv)), expand = c(0.02, 0.02))+
    guides(color = guide_legend(ncol = Leg.ncol), fill = guide_legend(ncol = Leg.ncol))+
    theme(axis.title.x = Tit.x, axis.ticks.x = Ticks.x, axis.text.x = Lab.x, strip.text.x = Strip.x, strip.text.y = Strip.y, legend.position = Leg.pos)
  
  #### Plot counts ####
  p2 <- ggplot(Res, aes(x = Age, y = N, color = Group)) + geom_segment(aes(yend = 0, xend = Age), alpha = Bar.alpha, linewidth = Bar.size) + 
    facet_grid(. ~ Param.clim) + 
    Col_scale + xlab("Time (cal. year BP)") +
    scale_y_continuous(expand = c(0.2, 0.5))+
    guides(color = guide_legend(ncol = Leg.ncol), fill = guide_legend(ncol = Leg.ncol))+
    theme(legend.position = "none", strip.text = element_blank(), axis.title.x = Tit.x2, axis.text.x = Lab.x2)
  
  if(Add.count.per.bin == T){p <- p / p2 + plot_layout(heights = c(8, 2)) & My_glob_theme}
  else{p <- p + My_glob_theme}
  
  #### Plot Maps ####
  if(Display.map == T){
    if(exists("ACA.bo") == F){source("Scripts/Load_maps.R")}
    label_longitude <- function(x) {paste0(abs(x), "°", ifelse(x < 0, "W", "E"))}
    label_latitude <- function(y) {paste0(abs(y), "°", ifelse(y < 0, "S", "N"))}
    
    if(Add.emptysite == T){
      if(Empty.site.from.Meta.data == T){
        Empty.site <- geom_point(data = Meta.data, mapping = aes(x = Longitude, y = Latitude),  alpha = 1, size = Map.dot.size*.2, shape = 16, color = "grey10")
        
      }
      else{
        if(exists("Meta.data.ACA") == T){
          Empty.site <- geom_point(data = Meta.data.ACA, mapping = aes(x = Longitude, y = Latitude),  alpha = 1, size = Map.dot.size*.2, shape = 16, color = "grey10")
        }
        else{print("Impossible to plot the empty dots, the 'Meta.data.ACA' file not found.")
          Empty.site <- NULL}
      }
    }
    else{Empty.site <- NULL}
    
    names(Meta.data)[names(Meta.data) == Group] <- "Group"
    pmap.groups <- map_empty +
      geom_polygon(data = ACA.bo.proj, aes(x=long, y=lat), colour="black", alpha = 1, fill = NA, size = 0.5)+
      coord_quickmap(xlim = c(Map.limits[1], Map.limits[2]), ylim = c(Map.limits[3], Map.limits[4])) +
      new_scale_colour()+
      Empty.site + facet_map +
      guides(color = guide_legend(ncol = Leg.ncol), fill = guide_legend(ncol = Leg.ncol))+
      geom_point(data = Meta.data, mapping = aes(x = Longitude, y = Latitude,  fill = Group),  alpha = 1, size = Map.dot.size, shape = 21, color = "black")+
      Fill_scale + Col_scale +
      scale_x_continuous(labels = label_longitude) +
      scale_y_continuous(labels = label_latitude) +
      theme(plot.background = element_blank(), panel.grid = element_blank(), panel.border = element_rect(colour = "black", linewidth = 1),
            legend.position = Map.leg.pos, legend.justification = c(0, 1), strip.background = element_blank(),
            panel.background = element_rect(fill = "#7BACCF"), axis.text.x = Lab.x2,
            strip.text = Strip.map, 
            legend.background = element_rect(fill = alpha("white", 1), color = "black"), legend.box.background = element_blank(),
      )
    
    if(is.null(Width.pour.map) == F){p <- (p | pmap.groups) + plot_layout(widths = c((100-Width.pour.map), Width.pour.map))}
    else{p <- (p | pmap.groups)}
  }
  
  
  #### Return ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(p, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", limitsize = F)}
    else{ggsave(Save.plot)}}
  if(is.null(Save.plot.Rds) == F){saveRDS(p, Save.plot.Rds)}
  
  if(is.null(Save.data.Rds) == F){saveRDS(Res.save, Save.data.Rds)}
  if(Return.plot == T){return(p)}
  else{return(Res.save)}
}

Compar.stack.inso <- function(Mdata = NULL, Minso = NULL, List.keep = NULL, Show.strip = T, Add.vlines = NULL, Remove.ages = NULL,
                              Rescale.reconstruction = T, Add.CI = T, CI.perc = 95, CI.alpha = .4, Facet.order = NULL, Leg.ncol = 1, Plot.2.merge = F, Leg.pos = "right",
                              Title.size = 14, Select.interv = 1000, Title.age = "Time (cal. year BP)", Limits = NULL, 
                              Time.lab.rot = 45, Add.error.zone = F, Error.alpha = .3, Line.inso.width = 2, Reverse.val = NULL,
                              Save.plot = NULL, H = 500, W = 500, Save.data.Rds = NULL, Save.plot.Rds = NULL, Return.plot = T) {
  #### Color settings ####
  Leg.tit <- "Orbital insolation"
  my_col_group <- c("1. Data" = "royalblue", "2. Model" = "#b37878ff", "Error" = "grey10",
                    "ECA" = "darkorange", "WCA" = "royalblue", "Winter" = "royalblue", "Spring" = "chartreuse3", "Summer" = "darkred", "Autumn" ="#b37878ff",
                    "High-elevation" = "#AB5B33", "Low-elevation" = "#86AE61", 
                    "Pollen" = "#BEA33A", "GDGT" = "#CF5612", "Speleothem" = "#47BA62", "Lake_level" = "royalblue", "Loess" = "brown4", "n-alkane" = "#86AE61",
                    "MAAT" = "darkred", "MAP" = "royalblue", "2. Model" = "#b37878ff", "3. Data" = "pink",
                    "Iran" = "#219F49", "Caucasus" = "#F8AA00", "Caspian" = "#10B0C8", "West. TSH" = "#A250A4", "East. TSH" = "#DF2E26",
                    "1. Data.MAAT" = "darkred", "1. Data.MAP" = "darkblue", "2. Model.MAAT" = "#b37878ff", "2. Model.MAP" = "#646499ff",
                    "1. Data.Temperature" = "darkred", "1. Data.Moisture" = "darkblue", "2. Model.Temperature" = "#b37878ff", "2. Model.Moisture" = "#646499ff")
  
  # my_col_group <- my_col_group[c(names(Minso), names(Mdata))]
  Col_scale <- scale_color_manual(values = my_col_group, name = Leg.tit)
  Fill_scale <- scale_fill_manual(values = my_col_group, name = Leg.tit)
  
  #### Clean and merge data ####
  Minso[c(2:ncol(Minso))] <- as.data.frame(lapply(Minso[c(2:ncol(Minso))], Scaling.fun))
  if(is.null(Reverse.val) == F){Minso[names(Minso) %in% Reverse.val] <- - Minso[names(Minso) %in% Reverse.val]}
  Minso <- reshape2::melt(Minso, id = "Age", variable.name = "Inso", value.name = "Inso.val")
  print(Minso)
  Mdata <- suppressWarnings(full_join(Mdata, Minso, by = "Age"))
  Mdata <- Mdata %>%
    tidyr::pivot_longer(cols = c(MAAT, MAP), names_to = "Param", values_to = "Value") %>%
    dplyr::mutate(SD = ifelse(Param == "MAAT", SD_MAAT, SD_MAP)) %>%
    dplyr::select(Age, Group, Subgroup, Param, Value, SD, Inso, Inso.val, N)
  Mdata <- Mdata[!is.na(Mdata$Value),]
  
  #### Select particular curves ####
  if(is.null(List.keep) == F){
    Row.to.keep <- c()
    for (i in 1:length(List.keep)) {
      Ti <- List.keep[[i]]
      Ri <- which(Mdata$Param == Ti[1] & Mdata$Group == Ti[2] & Mdata$Inso == Ti[3])
      Row.to.keep <- append(Row.to.keep, Ri)
    }
    Mdata <- Mdata[Row.to.keep,]
  }
  
  #### Remove ages ####
  if(!is.null(Remove.ages)){
    Mdata <- Mdata[!Mdata$Age %in% Remove.ages,]
  }
  
  #### Rescale ####
  if(Rescale.reconstruction == T){print("zup")
    Mdata <- Mdata %>%
      group_by(Group, Subgroup, Param) %>%
      mutate(Value = Scaling.fun(Value)) %>%
      ungroup()}
  
  #### CI ####
  if(Add.CI == T){
    Mdata$se <- Mdata$SD / sqrt(Mdata$N)  # Standard error
    Mdata$t_crit <- qt((1 + CI.perc/100)/2, pmax(Mdata$N - 1, 1))  # avoid df=0
    Mdata$CI_lower <- Mdata$Value - Mdata$t_crit * Mdata$se
    Mdata$CI_upper <- Mdata$Value + Mdata$t_crit * Mdata$se
    CI <- geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), color = NA, size = 0.15, linetype = "dashed", alpha = CI.alpha, show.legend = F)
  }
  else{
    CI <- NULL
  }
  
  #### Facet by groups ####
  if(is.null(Facet.order) == F){
    Mdata$Group <- factor(Mdata$Group, levels = Facet.order)
  }
  #### Graphical settings ####
  My_glob_theme <- theme(panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid = element_blank(), legend.key = element_blank(),
                         strip.background = element_blank(), plot.margin = unit(c(0,0,0,0), "pt"))
  
  if(is.null(Add.vlines) == F){My_vlines <- geom_vline(xintercept = Add.vlines, col = "grey30", lty = 2, alpha = 0.7)}
  else{My_vlines <- NULL}
  
  if(Plot.2.merge == T){Tit.x <- element_blank(); Lab.x <- element_blank()}
  else{Tit.x <- element_text(); Lab.x <- element_text()}
  if(is.null(Limits)){Limits <- c(round(min(Mdata$Age, na.rm = T), digits = -1), round(max(Mdata$Age, na.rm = T), digits = -1))}
  Strip.x <- element_text(size = Title.size, face = "bold")
  if(Show.strip == T){Strip.y <- element_text(size = Title.size, face = "bold")}
  else{Strip.y <- element_blank()}
  Lab.x <- element_text(angle = Time.lab.rot, hjust = 1, vjust = 1)
  
  #### Plot ####
  p <- ggplot(Mdata, aes(x = Age, y = Value, color = Param, fill = Param))+ 
    My_vlines +
    geom_line(show.legend = F)+
    facet_grid(Group ~ Param)+ geom_point(show.legend = F)+ CI+
    geom_line(inherit.aes = F, aes(x = Age, y = Inso.val, color = Inso), linetype = "dashed", linewidth = Line.inso.width, alpha = 1) +
    guides(color = guide_legend(ncol = Leg.ncol), fill = guide_legend(ncol = NULL))+
    scale_x_continuous(name = Title.age, breaks = round(seq(Limits[1], Limits[2], by = Select.interv)), expand = c(0.02, 0.02))+
    Col_scale + Fill_scale + ylab("z-scores")+
    theme(axis.title.x = Tit.x, axis.text.x = Lab.x, strip.text.x = Strip.x, strip.text.y = Strip.y, legend.position = Leg.pos)+
    My_glob_theme
  
  #### Return ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){ggsave(p, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm", limitsize = F)}
    else{ggsave(Save.plot)}}
  if(is.null(Save.plot.Rds) == F){saveRDS(p, Save.plot.Rds)}
  
  if(is.null(Save.data.Rds) == F){saveRDS(Res.save, Save.data.Rds)}
  if(Return.plot == T){return(p)}
  else{return(Mdata)}
  
}


Prep.2.multiproxy.plot <- function(L, MD, Param.clim = "MAAT", Keep.CI = F){
  L <- Filter(Negate(is.null), L)
  if(Keep.CI == F){
    L <- lapply(L, function(x) x <- x[names(x) %in% c("Age", Param.clim)])
    L <- Filter(function(df) !all(is.na(df[Param.clim])), L)
  }
  else{
    L <- lapply(L, function(x) x <- x[names(x) %in% c("Age", names(x)[grep(Param.clim, names(x))])])
    L <- Filter(function(df) !all(is.na(df[grep(Param.clim, names(df))])), L)
  }
  names(L) <- paste(MD$Site_full[match(names(L), MD$ID)], MD$Proxy_classes[match(names(L), MD$ID)], sep = ".")
  return(L)
}


Table.Latex.settings <- function(Table = NULL, Group.order = NULL, Round.coord = 3, Keep.columns = NULL, Rename.columns = NULL, Save.path = NULL, Abbrev.OK = F,
                                 Print.tab.head = F,
                                 Save.csv = F, Save.Latex = T, Names.2.rows = F, Blanck.cells = T, LateX.caption = NULL, Latex.lab = NULL) {
  
  #### Subfunctions ####
  add_ID <- function(data, Group.order) {
    # Generate prefixes A, B, C, ...
    prefixes <- LETTERS[seq_along(Group.order)]
    names(prefixes) <- Group.order
    data$ID <- with(data, paste0(
      prefixes[Area],
      ave(seq_along(Area), Area, FUN = seq_along)
    ))
    
    return(data)
  }
  
  blank_duplicates <- function(df) {as.data.frame(lapply(df, function(col) {c(col[1], ifelse(col[-1] == head(col, -1), "-", col[-1]))}), stringsAsFactors = FALSE)}
  
  #### Cleaning tables ####
  if(any(grepl("\\|TRA$", Table$ID))){Table <- Table[grepl("\\|TRA$", Table$ID),]}
  
  if(is.null(Group.order) == T){Group.order <- unique(Table$Area)}
  else{Table$Area <- factor(Table$Area, levels = Group.order)}
  Table <- Table[order(Table$Site_full),]
  Table <- Table[order(Table$Area),]
  
  Table <- add_ID(Table, Group.order = Group.order)
  Table$Area <- as.character(Table$Area)
  
  
  if(is.null(Round.coord) == F & is.null(Table["Latitude"]) == F){Table$Latitude <- as.character(round(Table$Latitude, digits = Round.coord))}
  if(is.null(Round.coord) == F & is.null(Table["Longitude"]) == F){Table$Longitude <- as.character(round(Table$Longitude, digits = Round.coord))}
  if(is.null(Keep.columns) == F){Table <- Table[Keep.columns]}
  if(is.null(Rename.columns) == F){names(Table) <- Rename.columns}
  if(Blanck.cells == T){Table[c(1:6)] <- blank_duplicates(Table[c(1:6)])}
  Table$Region[Table$Region == "-"] <- "" 
  
  Table$References[is.na(Table$References)] <- ""
  Table$References[Table$References == "My cores"] <- ""
  
  Table$References <- paste("\\citet{", gsub("  ", " ", Table$References), "}", sep = "")
  Table$`Proxy type` <- gsub("_", "-", Table$`Proxy type`)
  Table$`Proxy type` <- gsub("n-alkane", "\\\\textit{n}-alkane", Table$`Proxy type`)
  Table$`Proxy type` <- gsub("Paleodust", "MS, GS, SC", Table$`Proxy type`)
  Table$`Proxy type` <- gsub("\\+", "\\, ", Table$`Proxy type`)
  Table$`Proxy type` <- gsub("d18O", "\\$\\\\mathrm\\{\\\\delta\\^\\{18\\}O\\}\\$", Table$`Proxy type`)
  Table$`Proxy type` <- gsub("d13C", "\\$\\\\mathrm\\{\\\\delta\\^\\{13\\}C\\}\\$", Table$`Proxy type`)
  
  Table$Database[Table$Database == "No"] <- "Orig. pub."
  Table$Database[Table$Database == "LegacyPollen2.0, Neotoma"] <- "LegacyPollen2.0"
  Table$Database[grepl("EMBSeCBIO|Neotoma|Legacy|ISOMED|PANGAEA|Sisal", Table$Database)] <- paste("\\texttt{", Table$Database[grepl("EMBSeCBIO|Neotoma|Legacy|ISOMED|PANGAEA|Sisal", Table$Database)], "}", sep = "")
  Table$Database[Table$Database == "LLDB"] <- "\\citet{Street-Perrott_1989}"
  Table$Database[Table$Database == "LL_URSS_DB"] <- "\\citet{Tarasov_1996}"
  Table$Database[Table$Database == "Lin_2024"] <- "\\citet{Lin_2024}"
  Table$Database[Table$Database == "ALDS"] <- "\\citet{Li_2024b}"
  
  Table$`Climate proxy`[Table$`Climate proxy` == "Both"] <- "T, M"
  Table$`Climate proxy`[Table$`Climate proxy` == "MAAT"] <- "T"
  Table$`Climate proxy`[Table$`Climate proxy` == "MAP"] <- "M"
  
  if(Abbrev.OK == T){
    Table$`Archive type`[Table$`Archive type` == "Alluvial/Fluvial"] <- "Alluv./Fluv."
  }
  
  #### Save and export ####
  if(is.null(Save.path) == F){
    dir.create(file.path(dirname(Save.path), ""), showWarnings = F, recursive = T)
    if(Save.csv == T){
      Save.path.csv <- paste(Save.path, "csv", sep = ".")
      write.table(Table, file = Save.path.csv, row.names=T, col.names=NA, sep=",", dec = ".")
    }
    
    if(Save.Latex == T){
      library(xtable)
      Save.path.Latex <- paste(Save.path, "tex", sep = ".")
      Position.line <- which(Table$Region != "")
      Position.line <- c(-1,Position.line, nrow(Table)+1)
      
      if(Names.2.rows == T){
        L1 <- gsub(" .*", "", names(Table))
        L2 <- gsub(".* ", "", names(Table))
        Table <- rbind(L2, Table)
        names(Table) <- L1
        Table[1, which(Table[1,] == names(Table))] <- ""
        names(Table) <- paste("\\textbf{", names(Table), "}", sep = "")
        Position.line[Position.line == 0] <- 1
      }
      
      Tlatex <- xtable(Table, caption = LateX.caption, type = "latex", label = Latex.lab)
      print.xtable(Tlatex, file = Save.path.Latex, booktabs = T, include.rownames = F, comment = F,
                   caption.placement = "top", sanitize.text.function = function(x){x}, 
                   tabular.environment = "longtable", floating = F,
                   scale.box = '\textwidth',
                   hline.after = Position.line)
    }}
  
  if(Print.tab.head == T){print(head(Table))}
  return(Table)
}

Seasonal.insol <- function(Import.path = NULL, Save.Rds = NULL, Save.plot = NULL, Limits = NULL, Y.lim = NULL, Show.each.month = F,
                           Manual.y.line = NULL, Show.Facet = T, Keep.only = NULL, CI.lines = T,
                           Imported.raw.data = T, H = 800, W = 500, Remove.annual = F, Show.CI = T, Anomaly = F, Date.Anomaly = -100) {
  #### Settings ####
  if(is.null(Import.path) == T){warning("**** No data imported ! ****"); return(NULL)}
  Months <- c("Jan", "Fev", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  names(Months) <- seq(1,12)
  Season <- c(rep("Winter", 2), rep("Spring", 3), rep("Summer", 3), rep("Autumn", 3), "Winter")
  names(Season) <- seq(1,12)
  
  
  #### Import months data ####
  Iso.path <- list.files(path = Import.path, pattern = ".csv", full.names = T, all.files = T, recursive = F)
  Iso.path <- Iso.path[grepl(paste(Months, collapse = "|"), Iso.path)]
  
  if(Imported.raw.data == T){
    List.inso <- lapply(Iso.path, function(df){data.frame(read.table(file = df, header = F, sep = "", strip.white = T))})
    List.inso <- lapply(List.inso, function(x){names(x) <- c("Age", "eccentricity", "climatic_precession", "obliquity", "insolation"); return(x)})
    List.inso <- lapply(List.inso, function(x){x$Age <- -x$Age*1000; return(x)})
  }
  else{List.inso <- lapply(Iso.path, function(df){data.frame(read.csv(file = df, sep = "\t", dec=".", header = T))})}
  
  names(List.inso) <- gsub(".*_|\\.csv", "", Iso.path)
  List.inso <- List.inso[Months]
  Age <- List.inso[[1]]["Age"]
  
  List.inso <- Filter(Negate(is.null), List.inso)
  
  
  if(Anomaly == T){
    X.ano <- which(Age[1] == Date.Anomaly)
    List.inso <- lapply(List.inso, function(y){cbind(y[1], apply(y[2:ncol(y)], 2, function(x) x <- x - x[X.ano]))})
  }
  
  
  #### Calculate seasonal values ####
  for (i in 1:length(List.inso)) {List.inso[[i]] <- List.inso[[i]][names(List.inso[[i]]) == "insolation"]; names(List.inso[[i]]) <- i}
  List.inso <- do.call(cbind, List.inso)
  List.inso <- cbind(Age, List.inso)
  
  List.inso.m <- reshape2::melt(List.inso, "Age")
  List.inso.m$Season <- Season[List.inso.m$variable]
  List.inso.m$Months <- Months[List.inso.m$variable]
  
  Inso.annual <- List.inso.m
  Inso.sea.m <- aggregate(value ~ Age + Season, List.inso.m[c("Age", "Season", "value")], FUN = mean)
  Inso.sea.sd <- aggregate(value ~ Age + Season, List.inso.m[c("Age", "Season", "value")], FUN = sd)
  Inso.sea.m$MaxI <- Inso.sea.m$value + 1.86*Inso.sea.sd$value/3
  Inso.sea.m$MinI <- Inso.sea.m$value - 1.86*Inso.sea.sd$value/3
  
  Inso.annual.m <- aggregate(value ~ Age, Inso.annual[c("Age", "value")], FUN = mean)
  Inso.annual.sd <- aggregate(value ~ Age, Inso.annual[c("Age", "value")], FUN = sd)
  Inso.annual.m$MaxI <- Inso.annual.m$value + 1.86*Inso.annual.sd$value/3
  Inso.annual.m$MinI <- Inso.annual.m$value - 1.86*Inso.annual.sd$value/3
  Inso.annual.m$Season <- "Annual"
  Inso.annual.m <- Inso.annual.m[names(Inso.sea.m)]
  if(Remove.annual == F){
    Inso.plot <- rbind(Inso.sea.m, Inso.annual.m)
    Inso.plot$Season <- factor(Inso.plot$Season, levels = c("Annual", "Winter", "Spring", "Summer", "Autumn"))
  }
  else{
    Inso.plot <- Inso.sea.m
    Inso.plot$Season <- factor(Inso.plot$Season, levels = c("Winter", "Spring", "Summer", "Autumn"))
  }
  names(Inso.annual.m)[names(Inso.annual.m) == "value"] <- "Annual"
  
  #### Show each month ####
  if(Show.each.month == T){
    Add.months <- Inso.annual
    Add.months$MaxI <- NA
    Add.months$MinI <- NA
    Add.months <- Add.months[c(1,4,3,6,7,5)]
    Inso.plot$Months <- Inso.plot$Season
    Inso.plot <- rbind(Inso.plot, Add.months)
    Mean.curve <- geom_line(aes(y = value, x = Age, color = Months))
  }
  else{
    Inso.plot$Months <- NA
    Mean.curve <- geom_line(aes(y = value, x = Age))}
  
  #### Graph param ####
  Inso_color <- c(Summer = "darkred", Winter = "royalblue", Spring = "darkorange", Autumn = "#64A25F",
                  Jun = "chocolate",
                  
                  
                  Dec = "aquamarine", Jan = "aquamarine4", Fev = "cadetblue3",
                  Mar = "pink", Apr = "orange3", Avr = "orange4"
  )
  if(is.null(Limits) == F){My_lims_x <- xlim(Limits)}
  else{My_lims_x <- NULL}
  
  if(is.null(Y.lim) == F){My_lims_y <- ylim(Y.lim)}
  else{My_lims_y <- NULL}
  
  if(Show.CI == T){
    CI <- geom_ribbon(aes(ymin = MinI, ymax = MaxI, x = Age),  color = NA, alpha = .1)
    if(CI.lines == T){
      Min.CI <- geom_line(aes(y = MinI, x = Age), linetype = "dashed")
      Max.CI <- geom_line(aes(y = MaxI, x = Age), linetype = "dashed")}
    else{Min.CI <- NULL; Max.CI <- NULL}
    
    Mean.line <- "dashed"
  }
  else{CI <- NULL; Min.CI <- NULL; Max.CI <- NULL}
  
  if(is.null(Manual.y.line) == F){Event.clim <- geom_vline(xintercept = Manual.y.line, color = "grey70", linetype = "dashed")}
  else{Event.clim <- NULL}
  if(Show.Facet == T){My_facet <- facet_grid(vars(Season), scale = "free")}
  else{My_facet <- NULL}
  
  if(is.null(Keep.only) == F){Inso.plot <- Inso.plot[Inso.plot$Months %in% Keep.only,]}
  
  if(Anomaly == T){
    Line.zero <- geom_hline(yintercept = 0, color = "grey30", linetype = "dashed")
    Y.tit <- ylab(substitute(Insolation~anomaly~(W/m^2)))
  }
  else{
    Line.zero <- NULL
    Y.tit <- ylab(substitute(Insolation~(W/m^2)))}
  
  #### Plot season ####
  p.inso <- ggplot(data = Inso.plot, aes(fill = Season, color = Season, group = Months))+
    Line.zero + My_lims_x + My_lims_y + My_facet +
    Max.CI + Min.CI +
    Mean.curve +
    scale_color_manual(values = Inso_color)+
    scale_fill_manual(values = Inso_color)+
    CI + Event.clim + Y.tit + xlab("Time (cal. year BP)")+
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = "grey30"))
  
  if(is.null(Save.plot) == F){ggsave(filename = Save.plot, p.inso, width = W*0.026458333, height = H*0.026458333, units = "cm")}
  
  #### Export data ####
  Res.inso <- reshape2::dcast(Inso.sea.m[c("Age", "Season", "value")], Age ~ Season, value.var = "value")
  Res.inso <- cbind(Res.inso, Annual = Inso.annual.m$Annual)
  
  if(is.null(Save.Rds) == F){
    Path.to.create <- gsub("(.*/).*\\.Rds.*","\\1", Save.Rds)
    dir.create(file.path(Path.to.create), showWarnings = F, recursive = T)
    saveRDS(Res.inso, Save.Rds)
    }
  return(Res.inso)
}
