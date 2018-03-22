# AUXILIARY FUNCTIONS
# 1. conversion : input type conversion


# 1. aux.conversion -------------------------------------------------------
#' @keywords internal
#' @noRd
aux.conversion <- function(x){
  if (is.character(x)){
    x = as.numeric(as.factor(unlist(strsplit(x,split=""))))
  } else if (is.factor(x)){
    x = as.numeric(x)
  } else {
    x = as.numeric(as.factor(x))
  }
  return(round(x))
}






