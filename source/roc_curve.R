roc_curve = function(y_pred, y_true, threshold = NULL){
  print(y_true)
  if( !is.null(threshold) ){
    y_true = (y_true > threshold) + 0
  }
  if( !all(y_true %in% c(0,1)) ){
    stop("Only binary class assignment (0,1) in y_true is accepted!\n",
         "  - Did you swap y_true and y_pred in your function call?\n",
         "  - Do you need to set a threshold for defining neg/pos y_true values?")
  }
  if( all(y_true == 0) ){
    stop("Cannot calculate ROC-curve values when all y_true are negatives!")
  }
  if( all(y_true == 1) ){
    stop("Cannot calculate ROC-curve values when all y_true are positives!")
  }
  y_true = y_true[order(y_pred, decreasing=TRUE)]
  ROC = data.frame(FPR = cumsum( !y_true ) / sum( !y_true ),
                   TPR = cumsum( y_true ) / sum( y_true ))
  return( ROC )
}
