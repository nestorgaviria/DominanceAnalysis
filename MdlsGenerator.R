#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
#*******************************************************************************
#*****************   Function to generate the models ensemble  *****************
#************   this generates models using between 5-6 predictor  *************
#***********************   checks for duplicates  *************************
#*#***********************   ⬇︎⬇⬇︎⬇⬇︎⬇⬇︎⬇⬇︎⬇︎    ****************************** 
#≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈≈
#version 29.01.2026
# IMPROVED VERSION: Prevents duplicates

generate_models_ensemble <- function(pool, n_random_models, seed) {
  set.seed(seed)
  model_list <- list()
  
  # Track unique model signatures
  used_signatures <- character(0)
  
  # Helper function to create signature
  make_signature <- function(vars) {
    paste(sort(vars), collapse = "|")
  }
  
  # Helper function to add model (only if unique)
  add_model <- function(name, vars) {
    sig <- make_signature(vars)
    if(!(sig %in% used_signatures)) {
      model_list[[name]] <<- vars
      used_signatures <<- c(used_signatures, sig)
      global_counts[vars] <<- global_counts[vars] + 1
      return(TRUE)
    }
    return(FALSE)
  }
  
  # 1. SETUP
  all_preds <- unlist(pool, use.names = FALSE)
  pred_class_map <- data.frame(
    predictor = all_preds,
    class = rep(names(pool), lengths(pool)),
    stringsAsFactors = FALSE
  )
  
  # Global Usage Tracker
  global_counts <- setNames(rep(0, length(all_preds)), all_preds)
  class_names <- names(pool)

  print(paste(" Generating", n_random_models, "Interaction Models..."))
  
  # =========================================================
  # STEP 1: GENERATE BALANCED RANDOM MODELS (Interactions)
  # =========================================================
  
  successful_models <- 0
  max_attempts <- n_random_models * 5  # Allow more attempts
  attempt <- 0
  
  while(successful_models < n_random_models && attempt < max_attempts) {
    attempt <- attempt + 1
    
    # A. Determine Size (Range 3 to 9)
    target_size <- sample(5:6, 1, prob = c(0.5, 0.5))
    
    current_model <- c()
    current_classes <- c()
    
    
    #Force at least 1 predictor from each class first
    class_names <- names(pool)
    for(cls in class_names) {
      available_vars <- pool[[cls]]
      
      # Exclude already selected
      candidates <- setdiff(available_vars, current_model)
      
      if(length(candidates) == 0) {
        # This shouldn't happen, but fallback to any var from class
        candidates <- available_vars
      }
      
      # Weight by usage (favor underused)
      cand_usage <- global_counts[candidates]
      weights <- 1 / (1 + cand_usage)^2
      
      selected <- sample(candidates, 1, prob = weights)
      current_model <- c(current_model, selected)
      current_classes <- c(current_classes, cls)
      global_counts[selected] <- global_counts[selected] + 1
    }
    
    
    # B. Fill remaining slots (if target_size = 6, we need 1 more)
    while(length(current_model) < target_size) {
      
      candidates <- setdiff(all_preds, current_model)
      
      # Class Constraints (Max 2 per class)
      if(length(current_classes) > 0) {
        cls_counts <- table(current_classes)
        full_classes <- names(cls_counts)[cls_counts >= 2]
        cand_classes <- pred_class_map$class[match(candidates, pred_class_map$predictor)]
        candidates <- candidates[!(cand_classes %in% full_classes)]
      }
      
      if(length(candidates) == 0) break
      
      # Weights
      cand_usage <- global_counts[candidates]
      weights <- 1 / (1 + cand_usage)^2 
      
      selected <- sample(candidates, 1, prob = weights)
      current_model <- c(current_model, selected)
      
      selected_cls <- pred_class_map$class[match(selected, pred_class_map$predictor)]
      current_classes <- c(current_classes, selected_cls)
      global_counts[selected] <- global_counts[selected] + 1
    }
    
    # Try to add model (only if unique)
    if(add_model(paste0("Mdl_", successful_models + 1), current_model)) {
      successful_models <- successful_models + 1
    }
  }
  # 
  if(successful_models < n_random_models) {
    warning(paste("Only generated", successful_models, "unique models out of",
                  n_random_models, "requested"))
  }
  
  # =========================================================
  # DIAGNOSTICS
  # =========================================================
  usage_stats <- as.data.frame(global_counts)
  colnames(usage_stats) <- "Times_Selected"
  
  print("--- Balance Check (Top 5 & Bottom 5) ---")
  print(head(usage_stats[order(usage_stats$Times_Selected), , drop=FALSE], 5))
  print(tail(usage_stats[order(usage_stats$Times_Selected), , drop=FALSE], 5))
  
  print(paste("\n>>> Total unique models generated:", length(model_list)))
  
  return(model_list)
}