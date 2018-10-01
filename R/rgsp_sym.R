#' @name    rgsp_sym
#' @aliases rgsp_sym
#' @title Repetitive Group Sampling Plan Based on Cpk under Symmetric Case
#' @description Calculates Sample Number and Average Sample Number for Repetitive Group Sampling Plan based on Cpk under symmetric case as given in Aslam et al. (2013)
#'
#' @param .p1    Acceptable Quality Level (AQL) Probability
#' @param .p2    Limiting Quality Level (LQL) Probability
#' @param .alpha Producer's alpha-risk
#' @param .beta  Consumer's beta-risk
#' @param .nums  Number of samples with replacement at each iteration
#' @param .rep   Number of iterations
#'
#'
#'
#' @return Sample Number and Average Sample Number
#'
#' @author
#' \enumerate{
#'          \item Muhammad Yaseen (\email{myaseen208@@gmail.com})
#'          \item Muhammad Aslam (\email{aslam_ravian@@hotmail.com})
#'          \item Sami Ullah (\email{samiullahuos@@gmail.com})
#'          \item Muhammad Kashif (\email{mkashif@@uaf.edu.pk})
#'          }
#'
#' @references
#' Aslam, M., Wu, C., Jun, C., Azam, M. and Itay, N. (2013).
#'  Developing a variables repetitive group sampling plan based on process capability index Cpk with unknown mean and variance.
#'  \emph{Journal of Statistical Computation and Simulation}.
#'  \strong{83}(8):1507-1517. (\href{https://www.tandfonline.com/doi/abs/10.1080/00949655.2012.663374}{https://www.tandfonline.com/doi/abs/10.1080/00949655.2012.663374})
#'
#' @importFrom dplyr filter select
#' @importFrom magrittr %>%
#' @importFrom stats pnorm qnorm runif
#' @importFrom tibble tibble
#'
#' @export
#'
#' @examples
#'
#' rgsp_sym(
#'     .p1     = 0.0010
#'   , .p2     = 0.0020
#'   , .alpha  = 0.0500
#'   , .beta   = 0.1000
#'   , .nums   = 10000
#'   , .rep    = 10 # 1000
#' )
#'
#'
if(getRversion() >= "2.15.1"){
  utils::globalVariables(
    c(
        ".p1"
      , ".p2"
      , ".alpha"
      , ".beta"
      , ".nums"
      , ".rep"
      , "ASN"
      , "ASN2"
      , "LP1"
      , "LP2"
      , "ka"
      , "kr"
      , "n"
      , "p1"
      , "p2"
      , "r2"
    )
  )
}

rgsp_sym <- function(.p1, .p2, .alpha, .beta, .nums, .rep){
  UseMethod("rgsp_sym")
}

#' @export
#' @rdname rgsp_sym

rgsp_sym.default <- function(.p1, .p2, .alpha, .beta, .nums, .rep){

  zpu1 <- qnorm(p = .p1/2, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  zpl1 <- qnorm(p = .p1/2, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  C1   <- zpu1/3
  C2   <- zpl1/3

  zp1  <- qnorm(p = .p1/2, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  zp2  <- qnorm(p = .p2/2, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)

  Out0 <- tibble::tibble()
  for (i in 1:.rep){
    Out1  <- tibble::tibble(
        p1     = .p1
      , p2     = .p2
      , ka     = runif(n = .nums, min =  0, max = min(C1, C2))
      , kr     = runif(n = .nums, min =  0, max = ka)
      , n      = sample(x = 2:500, size = .nums, replace = TRUE)

      , LP1    = (2*pnorm(q = (zp1-3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1) - 1)/
                 (2*pnorm(q = (zp1-3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1) -
                  2*pnorm(q = (zp1-3*kr)*sqrt(n/(1+(9*(kr^2))/2)), mean = 0, sd = 1) + 1)

      , LP2    = (2*pnorm(q = (zp2-3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1) - 1)/
                 (2*pnorm(q = (zp2-3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1) -
                  2*pnorm(q = (zp2-3*kr)*sqrt(n/(1+(9*(kr^2))/2)), mean = 0, sd = 1) + 1)

      , r2      = 1/
                    (2*pnorm(q = (zp2-3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1) -
                     2*pnorm(q = (zp2-3*kr)*sqrt(n/(1+(9*(kr^2))/2)), mean = 0, sd = 1) + 1)

      , ASN2    = n/(2*pnorm(q = (zp2-3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1)-
                     2*pnorm(q = (zp2-3*kr)*sqrt(n/(1+(9*(kr^2))/2)), mean = 0, sd = 1) + 1)
    ) %>%

      dplyr::filter(LP1 >= (1-.alpha), LP2 <= .beta, r2 <= 3) %>%
      dplyr::filter(ASN2 == min(ASN2)) %>%
      dplyr::select(p1, p2, n, ka,  kr, ASN = ASN2)

    Out0 <- rbind(Out0, Out1)
    }
  Out <-
    Out0 %>%
    dplyr::filter(ASN == min(ASN))
  return(Out)
}

