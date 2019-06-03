library('tidyverse')
cross_fold_splits = function(partitions, print_splits = FALSE){
  
  # Visual helper function
  hline = function(width){
    return(str_c(rep_along(1:width, '_'), collapse = ''))
  }
  
  # Check partitions
  if( length(partitions) < 3 ){
    stop("At least 3 partitions are required")
  }
  
  # Set containing for final output
  splits = list()
  
  # Calculate the total number of splits (model yeilded)
  n_models = length(partitions) * (length(partitions) - 1)
  
  # Print split header to STDOUT
  if( print_splits ){
    cat(hline(80), "\n",
        "Based on (", paste(partitions,collapse = ", "),
        "), the following ",n_models," splits were created:\n",
        hline(80), "\n", sep='')
  }
  
  # Set model counter
  model_i = 0
  
  # Create splits
  for( test_partition in partitions ){
    
    # Remove test partition from all partitions to yield cross fold partitions
    cross_fold_partitions = setdiff(partitions, test_partition)
    for( validation_partition in cross_fold_partitions ){
      
      # Increment model counter
      model_i  = model_i + 1
      
      # Set id for current model
      zero = ''
      if( nchar(model_i) == 1 ){ zero = '0' }
      model_id = paste0("model_", zero, model_i)
      
      # Remove validation partition from cross fold partitions to yield
      # training partitions
      training_partitions = setdiff(cross_fold_partitions, validation_partition)
      
      # Initialise new model container
      splits[[model_id]] = list()
      
      # Set test partition
      splits[[model_id]][["test_partition"]]       = test_partition
      
      # Set validation partition
      splits[[model_id]][["validation_partition"]] = validation_partition
      
      # Set training partitions
      splits[[model_id]][["training_partitions"]]  = training_partitions
      
      # Print current composition to STDOUT
      if( print_splits ){
        cat("Test partition = ", test_partition,
            ", Validation partition = ", validation_partition,
            ", Training partitions = ", paste0(training_partitions,collapse = ', '),
            "\n",sep='')
      }
    }
  }
  # Print footer to STDOUT
  if( print_splits ){
    cat(hline(80), "\n", sep = '')
  }
  return(splits)  
}
data_splits = cross_fold_splits(partitions = 1:5, print_splits = TRUE)
