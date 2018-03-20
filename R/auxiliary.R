# 01. Adjusted Rand Index -------------------------------------------------
score01_adjrand <- function(confmat,scx,scy,n){
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

# 02. Chi-Squared ---------------------------------------------------------
score02_chisq <- function(confmat,scx,scy,n){
  # 2-1. preliminary
  nx = length(scx)
  ny = length(scy)

  # 2-2. main iteration
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

# 03. F-measure -----------------------------------------------------------
score03_f <- function(scx,scy,n){
  # 3-1. preliminary
  kk = length(scx)
  ll = length(scy)

  # 3-2. computation
  output = 0
  for (i in 1:kk){
    tx = scx[i]
    vecvaly = (2*tx*scy)/(tx+scy)
    output = output + (tx*max(vecvaly))
  }
  output = output/n
  return(output)
}

# 04. FMI -----------------------------------------------------------------
score04_fmi <- function(pairmat){
  # 4-1. separate
  n11 = pairmat[1,1]
  n01 = pairmat[2,1]
  n10 = pairmat[1,2]
  # 4-2. compute
  output = n11/sqrt((n11+n10)*(n11+n01))
  return(output)
}

# 05. Jaccard -------------------------------------------------------------
score05_jaccard <- function(pairmat){
  # 5-1. separate
  n11 = pairmat[1,1]
  n01 = pairmat[2,1]
  n10 = pairmat[1,2]
  # 5-2. compute
  output = n11/(n11+n10+n01)
  return(output)
}

# 06. MHM -----------------------------------------------------------------
score06_mhm <- function(confmat,n){
  # 6-1. get size
  kk = dim(confmat)[1]
  # 6-2. compute
  output = 0
  for (i in 1:kk){
    output = output+max(confmat[i,])
  }
  output = output/n
  return(output)
}

# 07. Mirkin --------------------------------------------------------------
score07_mirkin <- function(confmat,scx,scy){
  # 7-1. preliminary
  nk = length(scx)
  nl = length(scy)
  # 7-2. main iteration
  output = sum((scx^2)) + sum((scy^2))
  for (i in 1:nk){
    for (j in 1:nl){
      tgt = confmat[i,j]
      output = output - (2*(tgt^2))
    }
  }
  return(output)
}

# 08. MMM -----------------------------------------------------------------
score08_mmm <- function(confmat,n,nk,nl){
  # 8-1. preprocessing
  minsize = min(nk,nl)
  # 8-2. iteration
  output = 0
  for (i in 1:minsize){
    output = output + max(confmat[i,])
  }
  # 8-3. return
  output = output/n
  return(output)
}

# 09. nmi1 by Strehl & Ghosh ----------------------------------------------
score09_nmi1 <- function(Ixy,Hx,Hy,threps){
  altthr = c(threps,1e-7)
  denom  = sqrt(Hx*Hy)
  if (denom<threps){
    denom = max(altthr)
  }
  output = Ixy/denom

  if (output<1e-5){output = 0}
  if (output>1-(1e-5)){output = 1}

  return(output)
}

# 10. nmi2 by Fred & Jain -------------------------------------------------
score10_nmi2 <- function(Ixy,Hx,Hy,threps){
  altthr = c(threps,1e-7)
  denom  = Hx+Hy
  if (denom<threps){
    denom = max(altthr)
  }
  output = 2*Ixy/denom

  if (output<1e-5){output = 0}
  if (output>1-(1e-5)){output = 1}

  return(output)
}

# 11. overlap -------------------------------------------------------------
score11_overlap <- function(pairmat){
  x = pairmat[1,1]+pairmat[1,2]
  y = pairmat[1,1]+pairmat[2,1]

  output = pairmat[1,1]/min(x,y)
  return(output)
}

# 12. PD ------------------------------------------------------------------
score12_pd <- function(pairmat){
  output = pairmat[2,2]
  return(output)
}

# 13. Rand Index ----------------------------------------------------------
score13_rand <- function(pairmat,n){
  output = 2*(pairmat[2,2]+pairmat[1,1])/((n-1)*n)
  return(output)
}

# 14. Sorensen-Dice Coefficient -------------------------------------------
score14_sdc <- function(pairmat){
  TP = pairmat[1,1]
  FP = pairmat[1,2]
  FN = pairmat[2,1]

  output = (2*TP)/(2*TP+FN+FP)
  return(output)
}

# 15. Simple Matching Coefficient -----------------------------------------
score15_smc <- function(pairmat){
  output = (sum(diag(pairmat)))/(sum(pairmat))
  return(output)
}


# 16. tanimoto : special case of tversky with (alpha,beta)=(1,1) ----------
#  use score17_tversky with some special paramter choice

# 17. tversky -------------------------------------------------------------
#  double rcpp_tversky(NumericMatrix pairmat, const double alpha, const double beta, const int sym){
#  alpha=0.5=beta : Dice's Coefficient / SDC
#  alpha=1=beta   : Tanimoto Coefficient
#  sym=FALSE : original tversky
#  sym=TRUE  : a variant introduced on Wikipedia
score17_tversky <- function(pairmat,alpha,beta,sym){
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


# 18. VDM -----------------------------------------------------------------
score18_vdm <- function(confmat,n){
  # 18-1. preliminary
  nk = nrow(confmat)
  nl = ncol(confmat)

  # 18-2. iteration
  output = 2*n
  for (i in 1:nk){
    output = output-max(confmat[i,])
  }
  for (j in 1:nl){
    output = output-max(confmat[,j])
  }
  return(output)
}


# 19. VI ------------------------------------------------------------------
score19_vi <- function(Ixy,Hx,Hy,threps){
  # 19-1. prep
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

  # 19-2. compute
  output = valx+valy-2*Ixy
  if (output<1e-5){output = 0}
  if (output>1-(1e-5)){output = 1}
  return(output)
}




