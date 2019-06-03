# Clear workspace
# ------------------------------------------------------------------------------
rm(list=ls())
setwd("/home/matt/Documents/Courses/Algorithms_in_Bioinformatics/Project/")

# Load libraries
# ------------------------------------------------------------------------------
library('tidyverse')
library('ggseqlogo')
library('keras')

# Load functions
# ------------------------------------------------------------------------------
source(file = 'R/99_functions.R')
source(file = 'R/roc_curve.R')
source(file = 'R/auc.R')


# Set seeds for random number generation to ensure reproducibility
# ------------------------------------------------------------------------------
set.seed(1)
use_session_with_seed(1)

# Load data
# ------------------------------------------------------------------------------
pep_dat = read_table2(file = 'data/alpha_data_X_balanced_6_short.csv') 




# Let's see if we can match the QR codes of the 15-mer peptides to whether they
# bind to HLA-A*02:01 or not...

# Neural network training
# After exploring the training data we are now ready to train some neural
# networks to predict peptide binding to DRB1-03:01. The code below implements
# training of a feedforward neural network (FNN). Your task will be to add a
# convolutional neural network (CNN) and train this in parallel to the FNN.
# Then you should compare the FNNs and CNNs performance on the task of
# predicting peptide binding to MHCII molecule DRB1-03:01.
#
# But first please go through the code and make sure you understand what is
# going on.
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

build_CNN = function(k_pep, n_enc, n_hid, f_act, l_rate, f_loss){
  
  #   k_pep    Length of peptide
  #   n_end    Number of encoding values per amino acid in peptide
  #   n_hid    Number of hidden neurons
  #   f_act    Activation function
  #   l_rate   Learning rate
  #   f_loss   Loss function
  
  # Instantiate new keras sequential model
  mdl = keras_model_sequential()
  
  # Set layers and associated hyper parameters
  mdl %>% 
    layer_conv_1d(input_shape = c(k_pep, n_enc), filters = n_hid * 2, padding = "same",
                  kernel_initializer = "glorot_normal", kernel_size = 9, activation = f_act) %>%
    layer_global_max_pooling_1d() %>% 
    layer_dense(units = n_hid, activation = f_act) %>%
    layer_dense(units = 1, activation = f_act)
  
  # Compile model with chosen hyperparameters
  mdl %>% compile(loss = f_loss, optimizer = optimizer_adam(lr = l_rate), metrics = f_loss)
  
  # Done, return
  return(mdl)
}

pep_predict = function(peptide, model){
  pep_enc  = peptide %>% pep_encode_bl50 #%>% array_reshape(c(dim(.)[1], dim(.)[2] * dim(.)[3]))
  pep_pred = model %>% predict(pep_enc) %>% as.vector
  return(pep_pred)
}

# Set data variables
# ------------------------------------------------------------------------------

# Peptide length
k_pep = 129 #alpha
#k_pep = 138 #beta

# Number of encoding values per amino acid in peptide
n_enc = 23

#Set the hyperparameters
model_activation    = "sigmoid"
model_learning_rate = 0.001
model_optimizer     = optimizer_adam(lr = model_learning_rate)
model_loss          = 'mean_squared_error'
model_metrics       = 'mean_squared_error'
model_epochs        = 1  # add more
model_batch_size    = 20
model_n_hidden      = 2  # try with this

#Building the model (pre-defined in the function build_CNN)
CNN_model = build_CNN(k_pep, n_enc, n_hid = model_n_hidden, f_act = model_activation, l_rate = model_learning_rate, f_loss = model_loss)


#We will use partition 6 as independent test set and then do 5-fold cross validation on partitions 1-5
pep_dat_test_set = pep_dat %>% filter(partition == 6)
pep_dat = pep_dat %>% filter(partition != 6)

partitions  = pep_dat %>% count(partition) %>% pull(partition) # 1 2 3 4 5
data_splits = cross_fold_splits(partitions = partitions, print_splits = TRUE)

#Train the model
meta_data = list()
model_dir = 'models/'
for( mdl_i in names(data_splits) ){
  # Set partitions
  test_partition       = data_splits[[mdl_i]]$test_partition
  validation_partition = data_splits[[mdl_i]]$validation_partition
  training_partitions  = data_splits[[mdl_i]]$training_partitions
  
  # Set model file
  mdl_file = str_c(model_dir, mdl_i, ".hdf5")
  
  # Set callbacks used for early stopping
  callbacks_list = list( callback_early_stopping(monitor = "mean_squared_error", patience = 10),
                         callback_model_checkpoint(filepath = mdl_file, monitor = "val_loss",
                                                   save_best_only = TRUE) )
  # Extract and encode current training data
  pep_dat_training = pep_dat %>% filter(partition %in% training_partitions)
  x_train          = pep_dat_training %>% pull(peptide) %>% pep_encode_bl50 #%>% array_reshape(c(dim(.)[1], dim(.)[2] * dim(.)[3]))
  y_train          = pep_dat_training %>% pull(target) %>% array
  
  # Extract and encode current validation data
  pep_dat_validation = pep_dat %>% filter(partition %in% validation_partition)
  x_validation       = pep_dat_validation %>% pull(peptide) %>% pep_encode_bl50 #%>% array_reshape(c(dim(.)[1], dim(.)[2] * dim(.)[3]))
  y_validation       = pep_dat_validation %>% pull(target) %>% array
  
  # Do training
  history = CNN_model %>% fit(x = x_train, y = y_train, batch_size = model_batch_size, epochs = model_epochs,
                          callbacks = callbacks_list, validation_data = list(x_validation, y_validation))
  # Save meta data
  meta_data[[mdl_i]] = list(); meta_data[[mdl_i]]$mdl_file = mdl_file ; meta_data[[mdl_i]]$history = history
}


#TEST
for( mdl_id in names(meta_data) ){
  
  # Load current model  
  mdl = load_model_hdf5(filepath = meta_data[[mdl_id]]$mdl_file)
  meta_data[[mdl_id]]$model = mdl
  
  # Do prediction. Note how prediction are ONLY made for data points from the
  # test partition of the current model!
  pred_id = str_c('pred_', mdl_id)
  pep_dat = pep_dat %>% 
    mutate(!!pred_id := ifelse( partition == data_splits[[mdl_id]]$test_partition,
                                pep_predict(peptide = peptide, model = mdl), NA))
  
}
# Calculate mean predictions
pep_dat = pep_dat %>% mutate(pred_model_mean = select(., contains('pred_')) %>% rowMeans(na.rm = TRUE))

# Plot training
plot_training(history)

# Plot performance
plot_performance(CNN_model, X = x_validation, y = y_validation)

auc = auc(pep_dat$pred_model_mean,pep_dat$target)

roc = roc_curve(pep_dat$pred_model_mean,pep_dat$target)
