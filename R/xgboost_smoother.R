## Extract S for a single tree



create_S_from_tree =  function(current_tree_indices_train, current_tree_indices_test,lambda) {
  
  S_train = outer(current_tree_indices_train, current_tree_indices_train, FUN = "==")
  S_train = S_train/(rowSums(S_train)+lambda)
  
  
  
  S_test = outer(current_tree_indices_test, current_tree_indices_train, FUN = "==")
  S_test = S_test/(rowSums(S_test)+lambda)
  
  
  return(list(S_train = S_train, S_test = S_test))
}


## Extract S for a single boosted tree

create_S_from_single_boosted_tree = function(current_tree_indices_train,current_tree_indices_test,S_gb_prev,lambda){
  
  S  = create_S_from_tree(current_tree_indices_train,current_tree_indices_test,lambda)
  
  
  
  if (is.null(S_gb_prev)){
    
    # first tree: just normal tree
    return (list(S_train = S$S_train, S_test = S$S_test))
  }
  
  n_train = length(current_tree_indices_train)
  n_test = length(current_tree_indices_test)
  
  all_nodes = unique(current_tree_indices_train)
  n_nodes = length(all_nodes)
  
  
  node_corrections = matrix(0, nrow = n_nodes, ncol = n_train)
  S_train_correction = matrix(0, nrow = n_train, ncol = n_train)
  S_test_correction = matrix(0, nrow = n_test, ncol = n_train)
  
  
    
   
    
  for (i in 1:length(all_nodes)) {
      
      n = all_nodes[i]
      
      # Create correction matrix
      leaf_id_train = current_tree_indices_train == n
      node_corrections[i, ] = colSums(S_gb_prev[leaf_id_train, , drop = FALSE]) /(sum(leaf_id_train)+lambda)
      
      
      S_train_correction[leaf_id_train, ] = matrix(rep(node_corrections[i, ], sum(leaf_id_train)), nrow = sum(leaf_id_train), byrow = TRUE)
      
      leaf_id_test = current_tree_indices_test == n
      S_test_correction[leaf_id_test, ] = matrix(rep(node_corrections[i, ], sum(leaf_id_test)), nrow = sum(leaf_id_test), byrow = TRUE)
      
  }
  
  
    
  
  #cat("Dimensions of S_train:", dim(S$S_train), "\n")
  #cat("Dimensions of S_train_correction:", dim(S_train_correction), "\n")
  
  
  S_train = S$S_train - S_train_correction
  S_test = S$S_test - S_test_correction
  
 
  
  
  return(list(S_train = S_train, S_test = S_test))
  
}

# put everything together 

create_S_from_gbtregressor = function(model,leaf_indices_train,leaf_indices_test,base_score,output_dir, save_output, compact = FALSE){
  
  if (compact==FALSE){
    
    if (save_output){
      # Check if the directory exists, and create it if it doesn't
      if (!file.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
    }
    
    n_train = nrow(leaf_indices_train)
    n_test = nrow(leaf_indices_test)
    
    lr = model$params$eta
    
    #if (base_score == 0){
      
    S_train_curr = matrix(0, nrow = n_train, ncol = n_train)
    S_test_curr = matrix(0, nrow = n_test, ncol = n_train)
      
    #}
    
    #else if (base_score == "mean"){
     # print("yes")
     # S_train_curr = matrix(1/n_train, nrow = n_train, ncol = n_train)
     # S_test_curr = matrix(1/n_train, nrow = n_test, ncol = n_train)
      #return(list(S_train = S_train_curr ,S_test = S_test_curr))
    #}
      
    
    
  
  
    
    lambda = model$params$lambda
    
    number_of_trees = ncol(leaf_indices_train)
    
    for (col in 1:number_of_trees) {
      
      #cat(sprintf("\rProcessing tree: %d/%d", col, number_of_trees))
      
      #cat("\rProcessing tree:", col)
      #flush.console()
      
      current_tree_indices_train = leaf_indices_train[, col]
      current_tree_indices_test = leaf_indices_test[, col]
      
      S = create_S_from_single_boosted_tree(current_tree_indices_train,current_tree_indices_test,if(col == 1) NULL else S_train_curr , lambda)
      S_train_curr = S_train_curr + lr * S$S_train
      S_test_curr = S_test_curr + lr * S$S_test
      
      if (save_output) {
        # Save the current matrix to an RDS file
        output_filename_train = paste(output_dir, "/S_curr_iteration_train_", col, ".rds", sep = "")
        output_filename_test = paste(output_dir, "/S_curr_iteration_test_", col, ".rds", sep = "")
        
        saveRDS(S_train_curr, file = output_filename_train)
        saveRDS(S_test_curr, file = output_filename_test)
        
        #save the current tree matrix 
        
        output_filename_train = paste(output_dir, "/S_iteration_train_", col, ".rds", sep = "")
        output_filename_test = paste(output_dir, "/S_iteration_test_", col, ".rds", sep = "")
        
        saveRDS(S_train_curr, file = output_filename_train)
        saveRDS(S_test_curr, file = output_filename_test)
        
        
      }
    }
    
    
    
    return (list(S_train = S_train_curr ,S_test = S_test_curr))
  } 
  
}


get_xgboost_weights = function(model,dtrain,dtest,base_score = 0, save_output = FALSE, output_dir = NULL){
  
  model_params = model$params
  # Collect mismatched parameters
  errors = c()
  if (!is.null(model_params$alpha) && model_params$alpha != 0) {
    errors = c(errors, "alpha must be 0")
  }
  if (!is.null(model_params$subsample) && model_params$subsample != 1) {
    errors = c(errors, "subsample must be 1")
  }
  if (!is.null(model_params$max_delta_step) && model_params$max_delta_step != 0) {
    errors = c(errors, "max_delta_step must be 0")
  }
  if (!is.null(model_params$base_score) && model_params$base_score != 0) {
    errors = c(errors, "base_score must be 0")
  }
  
  if (length(errors) > 0) {
    stop(paste("The smoother is available only for specific parameter values:", paste(errors, collapse = ", ")))
  }
  
  
  leaf_indices_train = predict(model, dtrain, predleaf = TRUE)
  leaf_indices_test = predict(model, dtest, predleaf = TRUE)
  
  smoothers = create_S_from_gbtregressor(model,leaf_indices_train,leaf_indices_test,base_score = base_score,save_output = save_output,output_dir = output_dir)
  
  return(list(S_train = smoothers$S_train, S_test = smoothers$S_test))
  
}
  
  

# Metrics

compute_metrics_from_S = function(S, y) {
  # compute predictions
  y_pred = S %*% y
  
  # compute MSE
  mse = mean((y - y_pred)^2)
  
  # compute accuracy
  acc = mean((y_pred > 0) == (y > 0))
  
  # compute trace metric
  eff_p_tr = sum(diag(S))
  
  # compute l2-norm
  l2_norm = mean(sqrt(rowSums(S^2)))
  
  # compute squared l2-norm
  l2_norm_sq = mean(rowSums(S^2))
  
  return(list(mse = mse, acc = acc, eff_p_tr = eff_p_tr, l2_norm = l2_norm, l2_norm_sq = l2_norm_sq))
}


compute_metrics_from_xgboost = function(model,data,target){
  
  pred = predict(model, data)
  
  # Compute RMSE
  mse = mean((pred - target)^2)
  
  # compute accuracy
  acc = mean((pred > 0) == (target > 0))
  
  
  
  
  return(list(mse = mse,acc = acc))
  
  
  
}

reconstruct_predictions = function(bst, dtest) {
  # Get leaf indices from the model
  leaf_indices = predict(bst, dtest, predleaf = TRUE)
  
  # Extract booster model and tree data
  booster = xgb.Booster.complete(bst)
  trees = xgb.model.dt.tree(model = booster)
  
  # Initialize a matrix to store the predicted values for each sample
  predicted_values = matrix(0, nrow = nrow(leaf_indices), ncol = ncol(leaf_indices))
  reconstructed_predictions = numeric(nrow(leaf_indices))
  
  # Loop through each row (sample) in the leaf indices matrix
  for (i in 1:nrow(leaf_indices)) {
    sample_prediction = 0
    # Loop through each column (tree) in the leaf indices matrix
    for (j in 1:ncol(leaf_indices)) {
      # Get the leaf node index for the current sample and tree
      leaf_node = leaf_indices[i, j]
      
      # Find the corresponding row in the trees data.table
      leaf_row = trees[Tree == (j - 1) & Node == leaf_node]
      
      sample_prediction = sample_prediction + leaf_row$Quality
      
      # Extract the predicted value (Quality) and store it in the predicted_values matrix
      predicted_values[i, j] = leaf_row$Quality
    }
    
    reconstructed_predictions[i] = sample_prediction
  }
  
  return(reconstructed_predictions)
}




create_custom_objective <- function(subsample_rate, n_samples) {
  sampled_indices_per_tree <- list()
  
  custom_objective <- function(preds, dtrain) {
    # Get the true labels
    labels <- getinfo(dtrain, "label")
    
    # Initialize gradients and Hessians
    preds <- 1 / (1 + exp(-preds))  # Sigmoid function
    grad <- preds - labels
    hess <- preds * (1 - preds)
    
    # Subsampling
    sampled_indices <- sample(n_samples, size = floor(subsample_rate * n_samples), replace = FALSE)
    
    # Set gradients and Hessians of non-sampled data points to zero
    mask <- numeric(n_samples)
    mask[sampled_indices] <- 1
    grad <- grad * mask
    hess <- hess * mask
    
    # Record the sampled indices
    sampled_indices_per_tree[[length(sampled_indices_per_tree) + 1]] <<- sampled_indices
    
    return(list(grad = grad, hess = hess))
  }
  
  # Function to retrieve the sampled indices
  get_indices <- function() {
    sampled_indices_per_tree
  }
  
  list(
    objective = custom_objective,
    get_indices = get_indices
  )
}
