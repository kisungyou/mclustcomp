## CAT3 : Information Theory
#   12. single12_mi       : Mutual Information
#   13. single13_nmi1     : NMI by Strehl & Ghosh
#   14. single14_nmi2     : NMI by Fred & Jain
#   15. single15_vi       : Variation of Information
#   23. single23_jent     : Joint Entropy
#   24. single24_nmi3     : NMI by Danon
#   25. single25_nvi      : Normalized Variation of Information



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
  correct = min(threps)
  altthr = c(threps,1e-10)
  denom  = sqrt(Hx*Hy)
  if (denom<threps){
    denom = min(altthr)
  }
  output = Ixy/denom

  if (output<correct){output = 0}
  if (output>(1-correct)){output = 1}

  return(output)
}

# 14. single14_nmi2 -------------------------------------------------------
#' @keywords internal
#' @noRd
single14_nmi2 <- function(Ixy,Hx,Hy,threps){
  correct = min(threps)
  altthr = c(threps,1e-10)
  denom  = Hx+Hy
  if (denom<threps){
    denom = min(altthr)
  }
  output = 2*Ixy/denom

  if (output<correct){output = 0}
  if (output>(1-correct)){output = 1}

  return(output)
}

# 15. single15_vi ---------------------------------------------------------
#' @keywords internal
#' @noRd
single15_vi <- function(Ixy,Hx,Hy,threps){
  correct = min(threps)
  # 1. prep
  altthr = c(threps,1e-10)
  if (Hx<threps){
    valx = min(altthr)
  } else {
    valx = Hx
  }
  if (Hy<threps){
    valy = min(altthr)
  } else {
    valy = Hy
  }

  # 2. compute
  output = valx+valy-2*Ixy
  return(output)
}

# 23. single23_jent -------------------------------------------------------
# now from main script, 0-corrected matrix is given
#' @keywords internal
#' @noRd
single23_jent <- function(Pxy){
  # Compute the Value
  output = -sum(Pxy*log2(Pxy))
  return(output)
}

# 24. single24_nmi3 -------------------------------------------------------
#' @keywords internal
#' @noRd
single24_nmi3 <- function(Hx,Hy,Pxy){
  # 1. compute joint entropy from 'single23_jent'
  Hxy = -sum(Pxy*log2(Pxy))
  # 2. compute
  output = (Hx+Hy-Hxy)/((Hx+Hy)/2)
  return(output)
}

# 25. single25_nvi --------------------------------------------------------
#' @keywords internal
#' @noRd
single25_nvi <- function(Hx,Hy,Ixy,threps){
  # 1. prep
  altthr = c(threps,1e-10)
  if (Hx<threps){    Hx = min(altthr)  }
  if (Hy<threps){    Hy = min(altthr)  }

  # 2. compute
  output = (((Hx-Ixy)/Hx) + ((Hy-Ixy)/Hy))/2
  return(output)
}
