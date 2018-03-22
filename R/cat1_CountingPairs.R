## CAT1 : Counting Pairs
#   01. single01_chisq    : Chi-Squared Coefficient
#   02. single02_rand     : Rand Index
#   03. single03_adjrand  : Adjusted Rand Index
#   04. single04_fmi      : Fowlkes-Mallows Index
#   05. single05_mirkin   : Mirkin Metric
#   06. single06_jaccard  : Jaccard Index
#   07. single07_pd       : Partition Difference
#   16. single16_wallace1 : Wallace Criterion Type 1
#   17. single17_wallace2 : Wallace Criterion Type 2
#   18. single18_overlap  : Overlap Coefficient
#   19. single19_sdc      : Sorensen-Dice Coefficient
#   20. single20_smc      : Simple Matching Coefficient
#   22. single22_tversky  : Tversky Index


# 01. single01_chisq ------------------------------------------------------
#' @keywords internal
#' @noRd
single01_chisq <- function(confmat, scx, scy, n){
  # 1. preliminary
  nx = length(scx)
  ny = length(scy)

  # 2. main computation
  output = 0
  for (i in 1:nx){
    for (j in 1:ny){
      m_ij = confmat[i,j]
      e_ij = scx[i]*scy[j]/n
      output = output+ ((m_ij-e_ij)^2)/e_ij
    }
  }
  return(output)
}

# 02. single02_rand -------------------------------------------------------
#' @keywords internal
#' @noRd
single02_rand <- function(pairmat, n){
  n11 = pairmat[2,2]
  n00 = pairmat[1,1]
  output = (2*(n00+n11))/(n*(n-1))
  return(output)
}

# 03. single03_adjrand ----------------------------------------------------
#' @keywords internal
#' @noRd
single03_adjrand <- function(confmat,scx,scy,n){
  # 1-1. preprocessing
  nk = length(scx)
  nl = length(scy)
  # 1-2. computing basic elements
  t1 = 0
  for (i in 1:nk){
    tx = scx[i]
    t1 = t1+ (tx*(tx-1))/2
  }
  t2 = 0
  for (j in 1:nl){
    ty = scy[j]
    t2 = t2+ (ty*(ty-1))/2
  }
  t3 = (2*t1*t2)/(n*(n-1))
  summ = 0
  for (i in 1:nk){
    for (j in 1:nl){
      tgt = confmat[i,j]
      summ = summ+(tgt*(tgt-1))/2
    }
  }
  # 1-3. gathering up
  output = (summ-t3)/(((t1+t2)/2)-t3)
  return(output)
}

# 04. single04_fmi --------------------------------------------------------
#' @keywords internal
#' @noRd
single04_fmi <- function(pairmat){
  # 4-1. separate
  n11 = pairmat[1,1]
  n01 = pairmat[2,1]
  n10 = pairmat[1,2]
  # 4-2. compute
  output = n11/sqrt((n11+n10)*(n11+n01))
  return(output)
}

# 05. single05_mirkin -----------------------------------------------------
#' @keywords internal
#' @noRd
single05_mirkin <- function(confmat, scx, scy){
  # 1. preliminary
  nk = length(scx)
  nl = length(scy)
  # 2. main iteration
  output = sum((scx^2)) + sum((scy^2))
  for (i in 1:nk){
    for (j in 1:nl){
      tgt = confmat[i,j]
      output = output - (2*(tgt^2))
    }
  }
  return(output)
}

# 06. single06_jaccard ----------------------------------------------------
#' @keywords internal
#' @noRd
single06_jaccard <- function(pairmat){
  # 1. separate
  n11 = pairmat[1,1]
  n01 = pairmat[2,1]
  n10 = pairmat[1,2]
  # 2. compute
  output = n11/(n11+n10+n01)
  return(output)
}


# 07. single07_pd ---------------------------------------------------------
#' @keywords internal
#' @noRd
single07_pd <- function(pairmat){
  output = pairmat[2,2]
  return(output)
}

# 16. single16_wallace1 ---------------------------------------------------
#' @keywords internal
#' @noRd
single16_wallace1 <- function(pairmat, scx){
  n11   = pairmat[1,1]
  denom = sum((scx*(scx-1))/2)
  return(n11/denom)
}

# 17. single17_wallace2 ---------------------------------------------------
#' @keywords internal
#' @noRd
single17_wallace2 <- function(pairmat, scy){
  n11   = pairmat[1,1]
  denom = sum((scy*(scy-1))/2)
  return(n11/denom)
}

# 18. single18_overlap ----------------------------------------------------
#' @keywords internal
#' @noRd
single18_overlap <- function(pairmat){
  x = pairmat[1,1]+pairmat[1,2]
  y = pairmat[1,1]+pairmat[2,1]

  output = pairmat[1,1]/min(x,y)
  return(output)
}

# 19. single19_sdc --------------------------------------------------------
#' @keywords internal
#' @noRd
single19_sdc <- function(pairmat){
  TP = pairmat[1,1]
  FP = pairmat[1,2]
  FN = pairmat[2,1]

  output = (2*TP)/((2*TP)+FN+FP)
  return(output)
}

# 20. single20_smc --------------------------------------------------------
#' @keywords internal
#' @noRd
single20_smc <- function(pairmat){
  output = (sum(diag(pairmat)))/(sum(pairmat))
  return(output)
}

# 22. single22_tversky ----------------------------------------------------
#  Tanimoto Coefficient is special case of Tversky Index with (alpha,beta)=(1,1)
#  alpha=0.5=beta : Dice's Coefficient / SDC
#  alpha=1=beta   : Tanimoto Coefficient
#  sym=FALSE : original tversky
#  sym=TRUE  : a variant introduced on Wikipedia
#' @keywords internal
#' @noRd
single22_tversky <- function(pairmat,alpha,beta,sym){
  TP = pairmat[1,1]
  FP = pairmat[1,2]
  FN = pairmat[2,1]

  if (!sym){
    output = TP/(TP+(alpha*FP)+(beta*FN))
  } else {
    a = min(FP,FN)
    b = max(FP,FN)
    output = TP/(TP+(beta*(alpha*a+(1-alpha)*b)))
  }
  return(output)
}
