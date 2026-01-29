# This function executes a Dominance analysis on a list of models
# essentially return information about how much variance is explained
# by each predictor that is considered in each of the models in the list
# The analysis is done on a group basis, so it returns metrics about how much 
# variance does every group explains for each model. To facilitate aggregation
# For the analysis to be robust the list of models must contain models 
# that consider ideally minimun a predictor per group, also low numbers of 
# predictors are preferred. 

# after screening and fitting all models, the AIC criterion is used to select the top 10%
# of models and on those the importance is computed. 

#version 29.01.2026

Mdls_dominance_Analysis_OLS_AIC <- function(df, model_list, dataset_name, cfg) {
  
  # START TIMER
  start_time <- Sys.time()
  print(paste(">>> Processing:", dataset_name, "| N =", nrow(df)))
  print(paste("    Started at:", start_time))
  
  # =========================================================
  # 1. DATA TRANSFORMATION LAYER
  # =========================================================
  
  # A. Target Transformation
  y_col <- cfg$response_var
  if(cfg$response_transform == "log") {
    df$target_y <- log(df[[y_col]])
  } else {
    df$target_y <- df[[y_col]]
  }
  
  # B. Predictor Transformation (Log & Logit)
  log_vars <- cfg$Positive_to_Log
  logit_vars <- cfg$Proportions_to_Logit
  
  # Apply Log (with safety)
  for(col in log_vars) {
    if(col %in% names(df)) {
      vals <- df[[col]]
      if(any(vals <= 0, na.rm=TRUE)) {
        vals <- ifelse(vals <= 0, min(vals[vals>0], na.rm=TRUE)/2, vals)
      }
      df[[paste0("log_", col)]] <- log(vals)
    }
  }
  
  # Apply Logit (with safety)
  for(col in logit_vars) {
    if(col %in% names(df)) {
      vals <- df[[col]]
      if(max(vals, na.rm=TRUE) > 1) vals <- vals / 100
      epsilon <- 1e-4
      vals <- pmax(epsilon, pmin(1 - epsilon, vals))
      df[[paste0("logit_", col)]] <- log(vals / (1 - vals))
    }
  }
  
  df <- as.data.frame(df)
  
  
  # =========================================================
  # 2. PHASE 1: SCREENING LOOP (OLS + AIC)  
  # =========================================================
  
  print(paste("    >>> PHASE 1: Screening", length(model_list), "models using OLS + AIC..."))
  
  df_screening_results <- data.frame()
  n_failed <- 0  # Track failures
  
  # Loop through ALL models to get metrics
  for (m_name in names(model_list)) {
    tryCatch({
      raw_preds <- model_list[[m_name]]
      final_preds <- c()
      
      # Name Mapping Logic
      for(p in raw_preds) {
        if(p %in% log_vars && paste0("log_", p) %in% names(df)) {
          final_preds <- c(final_preds, paste0("log_", p))
        } else if(p %in% logit_vars && paste0("logit_", p) %in% names(df)) {
          final_preds <- c(final_preds, paste0("logit_", p))
        } else {
          final_preds <- c(final_preds, p)
        }
      }
      
      # Fit OLS regression
      model_data <- df[, c("target_y", final_preds), drop = FALSE]
      model_data <- na.omit(model_data)
      
      # Check if we have enough data after NA removal
      if(nrow(model_data) < length(final_preds) + 2) {
        warning(paste("Model", m_name, "skipped: insufficient data after NA removal"))
        n_failed <- n_failed + 1
        next
      }
      
      # Scale predictors
      model_data[, final_preds] <- scale(model_data[, final_preds])
      
      # Fit OLS
      formula_str <- paste("target_y ~", paste(final_preds, collapse = " + "))
      ols_model <- lm(as.formula(formula_str), data = model_data)
      
      res_row <- data.frame(
        Dataset     = dataset_name,
        Model_ID    = m_name,
        AIC         = AIC(ols_model),
        BIC         = BIC(ols_model),
        R2_adj      = summary(ols_model)$adj.r.squared,
        RMSE        = sqrt(mean(residuals(ols_model)^2)),
        LogLik      = as.numeric(logLik(ols_model)),
        N_obs       = nrow(model_data),  # Track sample size
        N_preds     = length(final_preds)
      )
      
      df_screening_results <- rbind(df_screening_results, res_row)
      
    }, error = function(e) {
      warning(paste("Model", m_name, "failed:", e$message))
      n_failed <- n_failed + 1
    })
  }
  
  print(paste("    >>> Screening complete.", nrow(df_screening_results), 
              "models succeeded,", n_failed, "failed."))
  
  # =========================================================
  # 3. SELECTION: PICK THE WINNERS  
  # =========================================================
  
  n_select <- ifelse(is.null(cfg$n_top_models), 50, cfg$n_top_models)
  
  df_winners <- df_screening_results %>%
    arrange(AIC) %>%  # Lower AIC = better
    slice_head(n = n_select)
  
  print(paste("    >>> SELECTION: Selected top", nrow(df_winners), "models based on AIC"))
  
  # =========================================================
  #  4. PHASE 2: DOMINANCE ANALYSIS ON WINNERS
  # =========================================================
  print("    >>> PHASE 2: Running Dominance Analysis on Winners...")
  
  df_dominance_results <- data.frame()
  
  for(i in 1:nrow(df_winners)) {
    
    m_name <- df_winners$Model_ID[i]
    raw_preds <- model_list[[m_name]]
    
    # Apply transformations
    final_preds <- c()
    for(p in raw_preds) {
      if(p %in% log_vars) {
        final_preds <- c(final_preds, paste0("log_", p))
      } else if(p %in% logit_vars) {
        final_preds <- c(final_preds, paste0("logit_", p))
      } else {
        final_preds <- c(final_preds, p)
      }
    }
    
    # Prepare data for this specific model
    model_data <- df[, c("target_y", final_preds), drop = FALSE]
    model_data <- na.omit(model_data)
    
    # Scale predictors (critical for fair dominance comparison)
    model_data[, final_preds] <- scale(model_data[, final_preds])
    
    # Fit OLS model for dominance analysis
    formula_str <- paste("target_y ~", paste(final_preds, collapse = " + "))
    ols_model <- lm(as.formula(formula_str), data = model_data)
    
    # Create groups for THIS specific model AND mapping for singles
    model_groups <- list()
    var_to_class_map <- list()  # Track which class each variable belongs to
    
    for(class_name in names(cfg$pool)) {
      class_vars_raw <- cfg$pool[[class_name]]
      
      # Find which transformed predictors belong to this class
      class_vars_transformed <- c()
      for(var in class_vars_raw) {
        transformed_var <- NULL
        
        if(paste0("log_", var) %in% final_preds) {
          transformed_var <- paste0("log_", var)
        } else if(paste0("logit_", var) %in% final_preds) {
          transformed_var <- paste0("logit_", var)
        } else if(var %in% final_preds) {
          transformed_var <- var
        }
        
        if(!is.null(transformed_var)) {
          class_vars_transformed <- c(class_vars_transformed, transformed_var)
          var_to_class_map[[transformed_var]] <- class_name  #Store mapping
        }
      }
      
      # Only include as group if it has 2+ predictors
      if(length(class_vars_transformed) >= 2) {
        model_groups[[class_name]] <- class_vars_transformed
      }
    }
    
    # Run dominance analysis
    tryCatch({
      lmg_result <- relaimpo::calc.relimp(
        ols_model, 
        type = "lmg",
        rela = TRUE,
        groups = if(length(model_groups) > 0) model_groups else NULL,
        groupnames = if(length(model_groups) > 0) names(model_groups) else NULL
      )
      
      # Extract results
      result_names <- names(lmg_result$lmg)
      result_values <- as.numeric(lmg_result$lmg)
      
      # Process each result
      for(j in seq_along(result_names)) {
        result_name <- result_names[j]
        result_value <- result_values[j]
        
        # Determine if this is a group or individual predictor
        if(result_name %in% names(model_groups)) {
          # It's a group
          class_name <- result_name
          n_in_class <- length(model_groups[[result_name]])
        } else {
          # It's an individual predictor - map back to its class
          class_name <- var_to_class_map[[result_name]]
          if(is.null(class_name)) {
            # Fallback: couldn't map (shouldn't happen)
            class_name <- result_name
            n_in_class <- 0
          } else {
            n_in_class <- 1
          }
        }
        
        # Store result
        group_importance <- data.frame(
          Model_ID = m_name,
          Dataset = dataset_name,
          Class = class_name,
          LMG_Importance = result_value * 100,
          Model_R2 = summary(ols_model)$r.squared,
          N_Predictors_Total = length(final_preds),
          N_Predictors_InClass = n_in_class,
          Is_Group = n_in_class >= 2,  # NEW: Flag for grouped vs individual
          stringsAsFactors = FALSE
        )
        
        df_dominance_results <- rbind(df_dominance_results, group_importance)
      }
      
    }, error = function(e) {
      warning(paste("Dominance analysis failed for model", m_name, ":", e$message))
    })
  }
  
  # =========================================================
  # 6. FINAL REPORTING
  # =========================================================
  
  end_time <- Sys.time()
  duration <- round(difftime(end_time, start_time, units = "mins"), 2)
  print(paste("    Finished at:", end_time))
  print(paste("    Total Execution Time:", duration, "minutes"))
  
  # Return results
  return(list(
    models = df_screening_results,
    winners = df_winners,
    dominance = df_dominance_results,
    dominance_summary = df_dominance_results %>%
      group_by(Dataset, Class) %>%
      summarise(
        Mean_Importance = mean(LMG_Importance),
        SD_Importance = sd(LMG_Importance),
        SE_Importance = sd(LMG_Importance) / sqrt(n()),  # Standard Error
        CI_lower = Mean_Importance - qt(0.975, n() - 1) * SE_Importance,  # 95% CI
        CI_upper = Mean_Importance + qt(0.975, n() - 1) * SE_Importance,
        
        Median_Importance = median(LMG_Importance),
        Min_Importance = min(LMG_Importance),
        Max_Importance = max(LMG_Importance),
   
        N_Models = n(),
        N_AsGroup = sum(Is_Group),
        N_AsSingle = sum(!Is_Group),
        .groups = 'drop'
      )
  ))
}

