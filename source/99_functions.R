# ------------------------------------------------------------------------------
# 36625 Algorithms in Bioinformatics - DTU Bioinformatics
# June 2018 - Leon Eyrich Jessen - jessen@bioinformatics.dtu.dk
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Various helper functions for running the exercises...
# ------------------------------------------------------------------------------



load_BLOSUM50 = function(){
  blosum_mat = read.table(file = 'data/BLOSUM50.txt')
  return(blosum_mat[1:23,1:23])
}



pep_check = function(pep){
  
  # Check if 'pep' is a character vector
  if( !( is.vector(pep) & is.character(pep) ) ){
    stop("'pep' has to be a character vector with one or more peptides!")
  }
  # Check if all sequences are of equal length
  if( pep %>% nchar %>% unique %>% length > 1 ){
    stop("All peptides must have equal number of positions")
  }
  # Check if 'pep' contain non-allowed amino acid residue characters
  if( pep %>% str_c(collapse='') %>% str_detect(pattern = "[^ARNDCQEGHILKMFPSTWYVX-]")){
    stop("Non standard amino acid residue character found.
         Only 'ARNDCQEGHILKMFPSTWYVX-' allowed")
  }
  return(TRUE)
}



pep_split = function(pep){
  
  # Check input
  pep_check(pep = pep)
  
  # Convert to matrix
  # do.call applies a function to the list returned from args
  # so rbind to form matrix each of the elements in the list returned
  # by strsplit
  return( do.call(what = rbind, args = strsplit(x = pep, split = '')) )
}



pep_encode_bl50 = function(pep){
  
  # Check input vector
  pep_check(pep = pep)
  
  # Set encoding matrix
  BLOSUM50 = load_BLOSUM50()
  enc_mat  = BLOSUM50 / 5
  enc_mat  = as.matrix(enc_mat)

  # Then we convert the vector of peptides to a matrix
  # with dimensions 'm x n' = 'n_peps x length_peps'
  p_mat = pep %>% pep_split
  
  # Assign meaningful variable names to dimensions
  n_peps = length(pep)     # i'th row
  l_peps = nchar(pep[1])   # j'th column
  l_enc  = ncol(enc_mat)   # k'th slice
  
  # Finally we define our output tensor as a 3d array
  # with dimensions n_peps x l_peps x l_enc (l_enc = 21)
  o_tensor = array(data = NA, dim = c(n_peps, l_peps, l_enc))
  for( i in 1:n_peps ){
    pep_i_residues = p_mat[i,]
    pep_img        = enc_mat[pep_i_residues,]
    o_tensor[i,,]  = pep_img
  }
  return(o_tensor)
}



plot_training = function(history){
  p = tibble(epoch = seq(1,history$params$epochs),
         training = history$metrics$mean_squared_error,
         validation = history$metrics$val_mean_squared_error) %>%
    gather(Set, MSE, -epoch) %>%
    ggplot(aes(x = epoch, y = MSE, colour = Set)) +
    geom_point() +
    geom_line() +
    theme_bw()
  return(p)
}



plot_combi_training = function(histories){
  p = tibble(epoch = seq(1, histories$FNN_history$params$epochs),
             FNN_training   = histories$FNN_history$metrics$mean_squared_error,
             FNN_validation = histories$FNN_history$metrics$val_mean_squared_error,
             CNN_training   = histories$CNN_history$metrics$mean_squared_error,
             CNN_validation = histories$CNN_history$metrics$val_mean_squared_error) %>% 
    gather(Set, MSE, -epoch) %>%
    ggplot(aes(x = epoch, y = MSE, colour = Set)) +
    geom_line() +
    theme_bw()
  return(p)
}



plot_performance = function(model, X, y){
  y_true = y %>% as.vector
  y_pred = model %>% predict(X) %>% as.vector
  pcc = cor(y_pred, y_true)
  mse = sum((y_true - y_pred) ** 2) / length(y_true)
  p = tibble(y_pred = y_pred, y_true = y_true) %>% 
    ggplot(aes(x = y_pred, y = y_true)) + 
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, colour = "black") +
    ggtitle(label = "", subtitle = str_c("Pearson's Correlation Coefficient = ", round(pcc,4),
                                         ", Mean Squared Error = ", round(mse,4))) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_bw()
  return(p)
}



plot_combi_performance = function(FNN_mdl, CNN_mdl, X, y){
  y_true = y %>% as.vector
  FNN_y_pred = FNN_mdl %>% predict(X) %>% as.vector
  CNN_y_pred = CNN_mdl %>% predict(X) %>% as.vector
  FNN_pcc = cor(FNN_y_pred, y_true, method = "pearson")
  CNN_pcc = cor(CNN_y_pred, y_true, method = "pearson")
  FNN_mse = sum((y_true - FNN_y_pred) ** 2) / length(y_true)
  CNN_mse = sum((y_true - CNN_y_pred) ** 2) / length(y_true)
  p = tibble(FNN_y_pred = FNN_y_pred, CNN_y_pred = CNN_y_pred, y_true = y_true) %>%
    gather(model, y_pred, -y_true) %>% 
    ggplot(aes(x = y_pred, y = y_true, colour = model)) + 
    geom_point() +
    geom_smooth(method = "lm") +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed', colour = "black") +
    ggtitle(label = "",
            subtitle = str_c("FNN: Pearson's Correlation Coefficient = ", round(FNN_pcc,4),
                             ", Mean Squared Error = ", round(FNN_mse,4), "\n",
                             "CNN: Pearson's Correlation Coefficient = ", round(CNN_pcc,4),
                             ", Mean Squared Error = ", round(CNN_mse,4))) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    theme_bw()
  return(p)
}



plot_pep_img = function(pep){
  
  # Check input vector
  pep_check(pep = pep)
  if( length(pep) > 1 ){
    stop("Only one peptide as input to plot_pep_img()")
  }
  
  # Set encoding matrix
  BLOSUM50 = load_BLOSUM50()
  enc_mat  = BLOSUM50 / 5
  enc_mat  = as.matrix(enc_mat)
  
  # Then we convert the vector of peptides to a matrix
  # with dimensions 'm x n' = 'n_peps x length_peps'
  p_mat = pep %>% pep_split
  
  # and use the residues to index the encoding matrix
  PSSM = enc_mat[p_mat[1,],]
  
  # Set names of peptide positions
  residues = paste(rownames(PSSM), seq(1, nrow(PSSM)), sep = '_')
  rownames(PSSM) = residues
  
  # Plot
  p = PSSM %>% as.data.frame %>% rownames_to_column %>% as_tibble %>%
    rename(peptide = rowname) %>% 
    mutate(peptide = factor(peptide, levels = rev(residues))) %>% 
    gather(encoding_value, value, -peptide) %>% 
    ggplot(aes(x = encoding_value, y = peptide, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    theme_bw()
  
  # Done, return
  return(p) 
  
}
