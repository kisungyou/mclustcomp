## CAT3 : Information Theory
#   12. single12_mi       : Mutual Information
#   13. single13_nmi1     : NMI by Strehl & Ghosh
#   14. single14_nmi2     : NMI by Fred & Jain
#   15. single15_vi       : Variation of Information



# 12. single12_mi ---------------------------------------------------------
#' @keywords internal
#' @noRd
single12_mi <- function(Ixy){
  return(Ixy)
}

# 13. single13_nmi1 -------------------------------------------------------
#' @keywords internal
#' @noRd
single13_nmi1 <- function(Ixy,Hx,Hy,threps){
  altthr = c(threps,1e-10)
  denom  = sqrt(Hx*Hy)
  if (denom<threps){
    denom = max(altthr)
  }
  output = Ixy/denom

  if (output<1e-5){output = 0}
  if (output>1-(1e-5)){output = 1}

  return(output)
}

# 14. single14_nmi2 -------------------------------------------------------
#' @keywords internal
#' @noRd
single14_nmi2 <- function(Ixy,Hx,Hy,threps){
  altthr = c(threps,1e-10)
  denom  = Hx+Hy
  if (denom<threps){
    denom = max(altthr)
  }
  output = 2*Ixy/denom

  if (output<1e-5){output = 0}
  if (output>1-(1e-5)){output = 1}

  return(output)
}

# 15. single15_vi ---------------------------------------------------------
#' @keywords internal
#' @noRd
single15_vi <- function(Ixy,Hx,Hy,threps){
  # 1. prep
  altthr = c(threps,1e-7)
  if (Hx<threps){
    valx = max(altthr)
  } else {
    valx = Hx
  }
  if (Hy<threps){
    valy = max(altthr)
  } else {
    valy = Hy
  }

  # 2. compute
  output = valx+valy-2*Ixy
  if (output<1e-5){output = 0}
  if (output>1-(1e-5)){output = 1}
  return(output)
}
