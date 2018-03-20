#' Measures for Comparing Clusterings
#'
#' Given two partitions \eqn{C_1} and \eqn{C_2}, it returns community comparison scores
#' corresponding with a set of designated methods. Note that two label vectors should be
#' of same length having either numeric or factor type.
#'
#'
#' @name mclustcomp
#' @aliases mclustcomp
#' @section list of the methods:
#' \describe{
#'   \item{\code{'adjrand'}}{\href{https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index}{"Adjusted Rand index"}. See [1985.Hubert].}
#'   \item{\code{'chisq'}}{\href{https://en.wikipedia.org/wiki/Chi-squared_test}{"Chi-squared coefficient"}. See [2001.Mirkin].}
#'   \item{\code{'f'}}{"F-Meausre". See [1999.Larsen].}
#'   \item{\code{'fmi'}}{\href{https://en.wikipedia.org/wiki/Fowlkes-Mallows_index}{"Fowlkes-Mallows index"}. See [1983.Fowlkes].}
#'   \item{\code{'jaccard'}}{\href{https://en.wikipedia.org/wiki/Jaccard_index}{"Jaccard index"}. See [1901.Jaccard].}
#'   \item{\code{'mhm'}}{"Meila-Heckerman Measure". See [2001.Meila].}
#'   \item{\code{'mirkin'}}{"Mirkin Metric", also known as "Equivalence Mismatch Distance". See [2000.van Dongen].}
#'   \item{\code{'mmm'}}{"Maximum-Match Measure". See [2006.Meila].}
#'   \item{\code{'nmi1'}}{\href{https://en.wikipedia.org/wiki/Mutual_information#Normalized_variants}{"Normalized Mutual Information"} by Strehl and Ghosh. See [2003.Strehl].}
#'   \item{\code{'nmi2'}}{\href{https://en.wikipedia.org/wiki/Mutual_information#Normalized_variants}{"Normalized Mutual Information"} by Fred and Jain. See [2003.Fred].}
#'   \item{\code{'overlap'}}{\href{https://en.wikipedia.org/wiki/Overlap_coefficient}{"Overlap coefficient"}. Also called as "Szymkiewicz-Simpson coefficient". See [1934.Szymkiewicz].}
#'   \item{\code{'pd'}}{"Partition Difference". See [2004.Li].}
#'   \item{\code{'rand'}}{\href{https://en.wikipedia.org/wiki/Rand_index}{"Rand index"}. See [1971.Rand].}
#'   \item{\code{'sdc'}}{\href{https://en.wikipedia.org/wiki/Sorensen-Dice_coefficient}{"Sørensen–Dice coefficient"}. Also known as "Sørensen index", "F1 score" or "Dice's coefficient". See [1948.Sørensen].}
#'   \item{\code{'smc'}}{\href{https://en.wikipedia.org/wiki/Simple_matching_coefficient}{"Simple Matching Coefficient"}. See [2007.Segaran].}
#'   \item{\code{'tanimoto'}}{\href{https://en.wikipedia.org/wiki/Jaccard_index#Tanimoto_similarity_and_distance}{"Tanimoto index"}. See [1958.Tanimoto].}
#'   \item{\code{'tversky'}}{\href{https://en.wikipedia.org/wiki/Tversky_index}{"Tversky index"}. Tanimoto coefficient (\code{'tanimoto'}) and Dice's coefficient (\code{'sdc'}) are special cases of Tversky index
#'   when (alpha,beta) = (1,1) and (0.5,0.5), respectively. See [1977.Tversky].}
#'   \item{\code{'vdm'}}{"van Dongen measure". See [2000.van Dongen].}
#'   \item{\code{'vi'}}{\href{https://en.wikipedia.org/wiki/Variation_of_information}{"Variation of Information"}. See [2003.Kraskov].}
#' }
#' @references [1901.Jaccard] Jaccard, P. (1901) \emph{Étude comparative de la distribution florale dans une portion des Alpes et des Jura.} Bulletin de la Société Vaudoise des Sciences Naturelles, 37:547-579.
#' @references [1934.Szymkiewicz] Szymkiewicz, D. (1934) \emph{Une contribution statistique a la geographie floristique.} Acta Societatis Botanicorum Poloniae, Vol.34(3):249-265.
#' @references [1948.Sørensen] Sørensen, T. (1948) \emph{A method of establishing groups of equal amplitude in plant sociology based on similarity of species and its application
#' to analyses of the vegetation on Danish commons.} Kongelige Danske Videnskabernes Selskab, Vol.5(4):1-34.
#' @references [1958.Tanimoto] Tanimoto, T. (1958) \emph{An Elementary Mathematical theory of Classification and Prediction.} Internal IBM Technical Report.
#' @references [1971.Rand] Rand, W.M. (1971) \emph{Objective criteria for the evaluation of clustering methods.} Journal of the American Statistical Association, Vol.66(336):846-850.
#' @references [1977.Tversky] Tversky, A. (1977) \emph{Features of Similarity.} Psychological Reviews, Vol.84(4):327-352.
#' @references [1983.Fowlkes] Fowlkes, E. B. and Mallows, C. L. (1983) \emph{A Method for Comparing Two Hierarchical Clusterings.} Journal of the American
#' Statistical Association, Vol.78(383):553-569.
#' @references [1985.Hubert] Hubert, L. and Arabie, P. (1985) \emph{Comparing partitions}. Journal of Classification, Vol.2(1):193-218.
#' @references [1999.Larsen] Larsen, B. and Aone, C. (1999) \emph{Fast and effective text mining using linear-time document clustering.}
#' Proceedings of the fifth ACM SIGKDD international conference on Knowledge discovery and data mining, 16-22.
#' @references [2000.van Dongen] van Dongen, S. (2000) \emph{Performance Criteria for Graph Clustering and Markov Cluster Experiment}. Centrum voor Wiskunde en Informatica,
#' Technical Report INS-R0012.
#' @references [2001.Meila] Meila, M. and Heckerman, D. (2001) \emph{An Experimental Comparison of Model-Based Clustering Methods}. Machine Learning, Vol.42(1-2):9-29.
#' @references [2001.Mirkin] Mirkin, B. (2001) \emph{Eleven Ways to Look at the Chi-Squared Coefficient for Contingency Tables}. The American Statistician, Vol.55(2):111-120.
#' @references [2003.Fred] Fred, A. L.N. and Jain, A. K. (2003) \emph{Robust Data Clustering.} Proceedings of IEEE Computer Society Conference on Computer Vision
#' and Pattern Recognition, CVPR, (3):128-136.
#' @references [2003.Kraskov] Kraskov, A., Stögbauer, H., Andrzejak, R.G., and Grassberger, P. (2003) \emph{Hierarchical Clustering Based on Mutual Information.} arXiv:q-bio/0311039.
#' @references [2003.Strehl] Strehl, A. and Ghosh, J. (2003) \emph{Cluster ensembles - a knowledge reuse framework for combining multiple partitions.} The Journal of Machine Learning Research, Vol.3:583-617.
#' @references [2004.Li] Li, T., Ogihara, M., and Ma, S. (2004) \emph{On combining multiple clusterings.} Proceedings of the thirteenth ACM international conference on Information and knowledge management, 294-303.
#' @references [2006.Meila] Meila, M. (2006) \emph{Comparing clusterings-an information based distance.} Journal of Multivariate Analysis, Vol.98(5):873-895.
#' @references [2007.Segaran] Segaran, T. (2007) \emph{Programming Collective Intelligence.} O'Reilly Media, ISBN-10:0596529325.
#'
#' @param x,y vectors of clustering labels
#' @param type.out \code{"all"} for returning scores for every available measure.
#' Either a single score name or a vector of score names can be supplied. See the section
#' for the list of the methods for details.
#' @param tversky.param a list of parameters for Tversky index; \code{alpha} and \code{beta} for
#' weight parameters, and \code{sym}, a logical where \code{FALSE} stands for original method, \code{TRUE}
#' for a revised variant to symmetrize the score. Default (alpha,beta)=(1,1).
#'
#' @return a data frame where each element is a score in accordance with a name of the method as a key.
#'
#' @examples
#' ## example 1. compare two identical clusterings
#' x = sample(1:5,10,replace=TRUE) # label from 1 to 5, 10 elements
#' y = x                           # set two labels x and y equal
#' mclustcomp(x,y)                 # show all results
#'
#' ## example 2. selection of a few methods
#' z = sample(1:4,10,replace=TRUE)           # generate a non-trivial clustering
#' cmethods = c("jaccard","tanimoto","rand") # select 3 methods
#' mclustcomp(x,z,type.out=cmethods)         # test with the selected scores
#'
#' ## example 3. tversky.param
#' tparam = list()                           # create an empty list
#' tparam$alpha = 2
#' tparam$beta  = 3
#' tparam$sym   = TRUE
#' mclustcomp(x,z,type.out="tversky")        # default set as Tanimoto case.
#' mclustcomp(x,z,type.out="tversky",tversky.param=tparam)
#'
#' @export
mclustcomp <- function(x,y,type.out="all",tversky.param=list()){
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


  ## Prelim1 : CONFUSION MATRIX of size(length(ux),length(uy))
  confmat = get.confusion(x,y,ux,uy)
  ## Prelim2 : size of each cluster
  scx = get.commsize(x,ux)
  scy = get.commsize(y,uy)
  ## Prelim3 : comembership matrix of (2,2)
  pairmat = get.pair(x,y)

  ## Prelim4 : probability-related stuffs for Mutual Information
  threps = min(1e-10,.Machine$double.eps)
  probs  = get.probs(confmat,scx,scy,n,threps)

  ## Control : type.out
  ## Case 1  : Single Argument
  ##  {"all" or single name}
  ## Case 2  : a vector of names; c("f","rand")
  type.allnames = c("adjrand","chisq","f","fmi","jaccard","mhm","mirkin","mmm","nmi1","nmi2","overlap","pd","rand","sdc","smc","tanimoto","tversky","vdm","vi")
  type.out   = unique(type.out)
  if ("all" %in% type.out){
    type.test = type.allnames
  } else {
    type.test = type.out  # this type test is the one we should generate again
  }

  ## Main Computation
  ## ("adjrand","chisq","f","fmi","jaccard","mhm","mirkin","mmm","nmi1","nmi2","overlap","pd","rand","sdc","smc","tanimoto","tversky","vdm","vi")
  type.score = array(0,c(1,length(type.test)))
  for (i in 1:length(type.test)){
    type.score[i] = mclustsingle(n,x,y,ux,uy,scx,scy,confmat,pairmat,probs,threps,type.test[i],tversky.param)
  }

  ## Return results
  result = data.frame(type.score)
  colnames(result) = type.test
  orderres = order(names(result)) # reordering by names
  result = result[orderres]

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
         "adjrand"  = {output = score01_adjrand(confmat,scx,scy,n)},
         "chisq"    = {output = score02_chisq(confmat,scx,scy,n)},
         "f"        = {output = score03_f(scx,scy,n)},
         "fmi"      = {output = score04_fmi(pairmat)},
         "jaccard"  = {output = score05_jaccard(pairmat)},
         "mhm"      = {output = score06_mhm(confmat,n)},
         "mirkin"   = {output = score07_mirkin(confmat,scx,scy)},
         "mmm"      = {output = score08_mmm(confmat,n,nk,nl)},
         "nmi1"     = {output = score09_nmi1(Ixy,Hx,Hy,threps)},
         "nmi2"     = {output = score10_nmi2(Ixy,Hx,Hy,threps)},
         "overlap"  = {output = score11_overlap(pairmat)},
         "pd"       = {output = score12_pd(pairmat)},
         "rand"     = {output = score13_rand(pairmat,n)},
         "sdc"      = {output = score14_sdc(pairmat)},
         "smc"      = {output = score15_smc(pairmat)},
         "tanimoto" = {output = score17_tversky(pairmat,1,1,FALSE)},
         "tversky"  = {output = score17_tversky(pairmat,t.alpha,t.beta,t.sym)},
         "vdm"      = {output = score18_vdm(confmat,n)},
         "vi"       = {output = score19_vi(Ixy,Hx,Hy,threps)}
         )
  # return output
  return(output)
}
