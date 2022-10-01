## CAT2 : Set Overlaps/Matching
#   08. single08_f        : F-Measure
#   09. single09_mhm      : Meila-Heckerman Measure
#   10. single10_mmm      : Maximum-Match Measure
#   11. single11_vdm      : Van Dongen Measure


# 08. single08_f ----------------------------------------------------------
#' @keywords internal
#' @noRd
single08_f <- function(scx,scy,n){
  # # 1. preliminary
  # kk = length(scx)
  # ll = length(scy)
  #
  # # 2. computation
  # output = 0
  # for (i in 1:kk){
  #   tx = scx[i]
  #   vecvaly = (2*tx*scy)/(tx+scy)
  #   output = output + (tx*max(vecvaly))
  # }
  # output = output/n
  # return(output)

  # 1. preliminary
  nx = length(scx)
  ny = length(scy)

  # 2. matrix entries
  fij = array(0, c(nx, ny))
  for (i in 1:nx){
    ci = scx[i]
    for (j in 1:ny){
      cj = scy[j]
      fij[i,j] = 2*ci*cj/(ci+cj)
    }
  }

  # 3. compute
  output = 0
  for (i in 1:nx){
    # output = output + (max(as.vector(fij[i,]))*scx[i]/n)
    output = output + (max(as.vector(fij[i,]))*scx[i]/base::sum(scx))
  }
  return(output)
}

# 09. single09_mhm --------------------------------------------------------
#' @keywords internal
#' @noRd
single09_mhm <- function(confmat, n){
  # 1. get size
  kk = dim(confmat)[1]
  # 2. compute
  output = 0
  for (i in 1:kk){
    output = output+max(confmat[i,])
  }
  output = output/n
  return(output)
}

# 10. single10_mmm --------------------------------------------------------
#' @keywords internal
#' @noRd
single10_mmm <- function(confmat,n,nk,nl){
  # 1. preprocessing
  minsize = min(nk,nl)
  # 2. iteration
  output = 0
  for (i in 1:minsize){
    output = output + max(confmat[i,])
  }
  # 3. return
  output = output/n
  return(output)
}

# 11. single11_vdm --------------------------------------------------------
#' @keywords internal
#' @noRd
single11_vdm <- function(confmat,n){
  # 1. preliminary
  nk = nrow(confmat)
  nl = ncol(confmat)

  # 2. iteration
  output = 2*n
  for (i in 1:nk){
    output = output-max(confmat[i,])
  }
  for (j in 1:nl){
    output = output-max(confmat[,j])
  }
  return(output)
}
