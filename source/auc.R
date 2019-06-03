auc = function(y_pred, y_true){
  roc = roc_curve(y_true = y_true, y_pred = y_pred)
  x  = roc$FPR
  y  = roc$TPR
  dx = c(diff(x), 0)
  dy = c(diff(y), 0)
  return( sum(y * dx) + sum(dy * dx)/2 )
}
