print.multinomAR <- function(out){
  
  #preliminaries
  cat("\nCall:\n")
  print(out$call)
  cat("\n")
  cat("Predicted total frequencies\n")
  print(out$mTot_new)
  cat("\n")
  cat("Predicted increase of total frequencies\n")
  print(out$mDiff_new)
  cat("\n")
  
}