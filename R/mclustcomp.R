#' Measures for Comparing Clusterings
#'
#' Given two partitions \eqn{C_1} and \eqn{C_2}, it returns community comparison scores
#' corresponding with a set of designated methods. Note that two label vectors should be
#' of same length having either numeric or factor type. Currently we have 3 classes of methods
#' depending on methodological philosophy behind each. See below for the taxonomy.
#'
#' @section Category 1. Counting Pairs:
#' \tabular{cl}{
#' TYPE \tab FULL NAME \cr
#' \code{'adjrand'}  \tab \href{https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index}{Adjusted Rand index}.\cr
#' \code{'chisq'}    \tab \href{https://en.wikipedia.org/wiki/Chi-squared_test}{Chi-Squared Coefficient}.\cr
#' \code{'fmi'}      \tab \href{https://en.wikipedia.org/wiki/Fowlkes-Mallows_index}{Fowlkes-Mallows index}.\cr
#' \code{'jaccard'}  \tab \href{https://en.wikipedia.org/wiki/Jaccard_index}{Jaccard index}.\cr
#' \code{'mirkin'}   \tab Mirkin Metric, or Equivalence Mismatch Distance. \cr
#' \code{'overlap'}  \tab \href{https://en.wikipedia.org/wiki/Overlap_coefficient}{Overlap Coefficient}, or Szymkiewicz-Simpson coefficient.\cr
#' \code{'pd'}       \tab Partition Difference.\cr
#' \code{'rand'}     \tab \href{https://en.wikipedia.org/wiki/Rand_index}{Rand Index}.\cr
#' \code{'sdc'}      \tab \href{https://en.wikipedia.org/wiki/Sorensen-Dice_coefficient}{Sørensen–Dice Coefficient}.\cr
#' \code{'smc'}      \tab \href{https://en.wikipedia.org/wiki/Simple_matching_coefficient}{Simple Matching Coefficient}.\cr
#' \code{'tanimoto'} \tab \href{https://en.wikipedia.org/wiki/Jaccard_index#Tanimoto_similarity_and_distance}{Tanimoto index}.\cr
#' \code{'tversky'}  \tab \href{https://en.wikipedia.org/wiki/Tversky_index}{Tversky index}. Tanimoto Coefficient and Dice's coefficient are special cases with (alpha,beta) = (1,1) and (0.5,0.5), respectively.\cr
#' \code{'wallace1'} \tab Wallace Criterion Type 1.\cr
#' \code{'wallace2'} \tab Wallace Criterion Type 2.
#' }
#'
#' @section Category 2. Set Overlaps/Matching:
#' \tabular{cl}{
#' TYPE \tab FULL NAME \cr
#' \code{'f'}   \tab F-Measure.
#' \code{'mhm'} \tab Meila-Heckerman Measure.
#' \code{'mmm'} \tab Maximum-Match Measure.
#' \code{'vdm'} \tab Van Dongen Measure.
#' }
#'
#' @section Category 3. Information Theory:
#' \tabular{cl}{
#' TYPE \tab FULL NAME \cr
#' \code{'mi'}   \tab Mutual Information. \cr
#' \code{'nmi1'} \tab \href{https://en.wikipedia.org/wiki/Mutual_information#Normalized_variants}{Normalized Mutual Information} by Strehl and Ghosh. \cr
#' \code{'nmi2'} \tab \href{https://en.wikipedia.org/wiki/Mutual_information#Normalized_variants}{Normalized Mutual Information} by Fred and Jain. \cr
#' \code{'vi'}   \tab \href{https://en.wikipedia.org/wiki/Variation_of_information}{Variation of Information}.
#' }
#'
#' @param x,y vectors of clustering labels
#' @param types \code{"all"} for returning scores for every available measure.
#' Either a single score name or a vector of score names can be supplied. See the section
#' for the list of the methods for details.
#' @param tversky.param a list of parameters for Tversky index; \code{alpha} and \code{beta} for
#' weight parameters, and \code{sym}, a logical where \code{FALSE} stands for original method, \code{TRUE}
#' for a revised variant to symmetrize the score. Default (alpha,beta)=(1,1).
#'
#' @return a data frame with columns \code{types} and corresponding \code{scores}.
#'
#' @examples
#' ## example 1. compare two identical clusterings
#' x = sample(1:5,20,replace=TRUE) # label from 1 to 5, 10 elements
#' y = x                           # set two labels x and y equal
#' mclustcomp(x,y)                 # show all results
#'
#' ## example 2. selection of a few methods
#' z = sample(1:4,20,replace=TRUE)           # generate a non-trivial clustering
#' cmethods = c("jaccard","tanimoto","rand") # select 3 methods
#' mclustcomp(x,z,types=cmethods)            # test with the selected scores
#'
#' ## example 3. tversky.param
#' tparam = list()                           # create an empty list
#' tparam$alpha = 2
#' tparam$beta  = 3
#' tparam$sym   = TRUE
#' mclustcomp(x,z,types="tversky")           # default set as Tanimoto case.
#' mclustcomp(x,z,types="tversky",tversky.param=tparam)
#'
#'
#' @references
#' \insertRef{strehl_cluster_2003}{mclustcomp}
#'
#' \insertRef{meila_comparing_2007}{mclustcomp}
#'
#' \insertRef{goos_comparing_2003}{mclustcomp}
#'
#' \insertRef{wagner_comparing_2007}{mclustcomp}
#'
#' \insertRef{albatineh_similarity_2006}{mclustcomp}
#'
#' \insertRef{mirkin_eleven_2001}{mclustcomp}
#'
#' \insertRef{rand_objective_1971}{mclustcomp}
#'
#' \insertRef{kuncheva_using_2004}{mclustcomp}
#'
#' \insertRef{fowlkes_method_1983}{mclustcomp}
#'
#' \insertRef{dongen_performance_2000}{mclustcomp}
#'
#' \insertRef{jaccard_distribution_1912}{mclustcomp}
#'
#' \insertRef{li_combining_2010}{mclustcomp}
#'
#' \insertRef{larsen_fast_1999}{mclustcomp}
#'
#' \insertRef{meila_experimental_2001}{mclustcomp}
#'
#' \insertRef{cover_elements_2006}{mclustcomp}
#'
#' \insertRef{ana_robust_2003}{mclustcomp}
#'
#' \insertRef{wallace_comment_1983}{mclustcomp}
#'
#' \insertRef{simpson_mammals_1943}{mclustcomp}
#'
#' \insertRef{dice_measures_1945}{mclustcomp}
#'
#' \insertRef{segaran_programming_2007}{mclustcomp}
#'
#' \insertRef{tversky_features_1977}{mclustcomp}
#'
#' @export
mclustcomp <- function(x,y,types="all",tversky.param=list()){
  #------------------------------------------------------------------------
  ## PREPROCESSING
  #   1. size argument
  if ((!is.vector(x))||(!is.vector(y))){
    stop("* mclustcomp : input 'x' and 'y' should both be a vector of class labels.")
  }
  n = length(x)
  if (length(y)!=n){
    stop("* mclustcomp : two vectors should be of same size.")
  }

  #   2. type conversion and unique vector
  x  = aux.conversion(x)
  y  = aux.conversion(y)
  ux = unique(x)
  uy = unique(y)

  if (length(ux)==1){    warning("* mclustcomp : 'x' is a trivial clustering.")  }
  if (length(uy)==1){    warning("* mclustcomp : 'y' is a trivial clustering.")  }
  if (length(ux)==n){    warning("* mclustcomp : 'x' is the singleton clustering.")  }
  if (length(uy)==n){    warning("* mclustcomp : 'y' is the singleton clustering.")  }

  #   3. tversky parameter
  listdot = as.list(environment())
  if ("tversky.param" %in% names(listdot)){
    tversky.param = listdot$tversky.param
  } else {
    tversky.param = list()
  }
  if (!("alpha" %in% names(tversky.param))){tversky.param$alpha = 1}
  if (!("beta" %in% names(tversky.param))){tversky.param$beta = 1}
  if (!("sym" %in% names(tversky.param))){tversky.param$sym = FALSE}
  if (tversky.param$alpha < 0){stop("* mclustcomp : tversky.param$alpha should be >= 0.")}
  if (tversky.param$beta < 0){stop("* mclustcomp : tversky.param$beta should be >= 0.")}
  if (!is.logical(tversky.param$sym)){stop("* mclustcomp : tversky.param$sym should
                                           be a logical variable; FALSE for original Tversky index, TRUE for a variant.")}


  #------------------------------------------------------------------------
  ## PRELIMINARY COMPUTATIONS
  ## Prelim1 : CONFUSION MATRIX of size(length(ux),length(uy))
  confmat = get.confusion(x,y,ux,uy)
  ## Prelim2 : size of each cluster
  scx = get.commsize(x,ux)
  scy = get.commsize(y,uy)
  ## Prelim3 : comembership matrix of (2,2)
  pairmat = get.pair(x,y)

  ## Prelim4 : probability-related stuffs for Mutual Information
  threps = min(1e-10,10*(.Machine$double.eps))
  probs  = get.probs(confmat,scx,scy,n,threps)

  ## Control : type.out
  ## Case 1  : Single Argument
  ##  {"all" or single name}
  ## Case 2  : a vector of names; c("f","rand")
  type_allnames = c("adjrand","chisq","f","fmi","jaccard","mhm","mirkin","mmm",
                    "mi","nmi1","nmi2","overlap","pd","rand","sdc","smc","tanimoto",
                    "tversky","vdm","vi","wallace1","wallace2")
  type_out   = unique(types)
  if ("all" %in% type_out){
    type_test = sort(type_allnames)
  } else {
    type_test = sort(type_out)  # this type test is the one we should generate again
  }

  #------------------------------------------------------------------------
  ## MAIN COMPUTATION
  type_score = rep(0,length(type_test))
  for (i in 1:length(type_test)){
    type_score[i] = mclustsingle(n,x,y,ux,uy,scx,scy,confmat,pairmat,probs,threps,type_test[i],tversky.param)
  }

  #------------------------------------------------------------------------
  ## RETURN RESULTS
  result = data.frame(types=type_test,scores=type_score)
  return(result)
}



# COMPUTE :: single measure branching -------------------------------------
## Original Implementation of 19 methods
mclustsingle <- function(n,x,y,ux,uy,scx,scy,confmat,pairmat,probs,threps,type,tversky.param){
  # Missing parameters for score08_mmm
  nk = length(scx)
  nl = length(scy)
  # Sepearting probs for NMI and VIs
  Ixy = probs$Ixy
  Hx  = probs$Hx
  Hy  = probs$Hy
  # Tversky parameter
  t.alpha = tversky.param$alpha
  t.beta  = tversky.param$beta
  t.sym   = tversky.param$sym

  switch(type,
         "chisq"    = {output = single01_chisq(confmat,scx,scy,n)},
         "rand"     = {output = single02_rand(pairmat,n)},
         "adjrand"  = {output = single03_adjrand(confmat,scx,scy,n)},
         "fmi"      = {output = single04_fmi(pairmat)},
         "mirkin"   = {output = single05_mirkin(confmat,scx,scy)},
         "jaccard"  = {output = single06_jaccard(pairmat)},
         "pd"       = {output = single07_pd(pairmat)},
         "f"        = {output = single08_f(scx,scy,n)},
         "mhm"      = {output = single09_mhm(confmat,n)},
         "mmm"      = {output = single10_mmm(confmat,n,nk,nl)},
         "vdm"      = {output = single11_vdm(confmat,n)},
         "mi"       = {output = single12_mi(Ixy)},
         "nmi1"     = {output = single13_nmi1(Ixy,Hx,Hy,threps)},
         "nmi2"     = {output = single14_nmi2(Ixy,Hx,Hy,threps)},
         "vi"       = {output = single15_vi(Ixy,Hx,Hy,threps)},
         "wallace1" = {output = single16_wallace1(pairmat,scx)},
         "wallace2" = {output = single17_wallace2(pairmat,scy)},
         "overlap"  = {output = single18_overlap(pairmat)},
         "sdc"      = {output = single19_sdc(pairmat)},
         "smc"      = {output = single20_smc(pairmat)},
         "tanimoto" = {output = single22_tversky(pairmat,1,1,FALSE)},
         "tversky"  = {output = single22_tversky(pairmat,t.alpha,t.beta,t.sym)}
         )

  # return output
  return(output)
}
