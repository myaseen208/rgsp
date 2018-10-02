#' @name    rgsp_asym2
#' @aliases rgsp_asym2
#' @title Repetitive Group Sampling Plan Based on Cpk under asymmetric Case 2
#' @description Calculates Sample Number and Average Sample Number for Repetitive Group Sampling Plan based on Cpk under asymmetric case 2 as given in Aslam et al. (2013)
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
#' rgsp_asym2(
#'     .p1     = 0.001
#'   , .p2     = 0.003
#'   , .alpha  = 0.050
#'   , .beta   = 0.100
#'   , .nums   = 10000
#'   , .rep    = 10 # 1000
#' )
#'
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
      , "D1"
      , "D2"
      , "N1"
      , "N2"
    )
  )
}

rgsp_asym2 <- function(.p1, .p2, .alpha, .beta, .nums, .rep){
  UseMethod("rgsp_asym2")
}

#' @export
#' @rdname rgsp_asym2

rgsp_asym2.default <- function(.p1, .p2, .alpha, .beta, .nums, .rep){

  zpl1  <- qnorm(p =   .p1/3, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  zpu1  <- qnorm(p = 2*.p1/3, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  C1    <- zpl1/3
  C2    <- zpu1/3

  zp11  <- qnorm(p =   .p1/3, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  zp12  <- qnorm(p = 2*.p1/3, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  zp21  <- qnorm(p =   .p2/3, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)
  zp22  <- qnorm(p = 2*.p2/3, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)

  Out0 <- tibble::tibble()
  for (i in 1:.rep){
    Out1  <- tibble::tibble(
        p1     = .p1
      , p2     = .p2
      , ka     = runif(n = .nums, min =  0, max = min(C1,C2))
      , kr     = runif(n = .nums, min =  0, max = ka)
      , n      = sample(x = 2:500, size = .nums, replace = TRUE)
      , N1     = pnorm(q =  (zp11 - 3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1) -
                 pnorm(q = -(zp12 - 3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1)
      , D1     = N1+ 1 - pnorm(q =  (zp11 - 3*kr)*sqrt(n/(1+(9*(kr^2))/2)), mean = 0, sd = 1) +
                         pnorm(q = -(zp12 - 3*kr)*sqrt(n/(1+(9*(kr^2))/2)), mean = 0, sd = 1)
      , LP1    = N1/D1

      , N2     = pnorm(q =  (zp21 - 3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1) -
                 pnorm(q = -(zp22 - 3*ka)*sqrt(n/(1+(9*(ka^2))/2)), mean = 0, sd = 1)

      , D2     = N2+ 1 - pnorm(q =  (zp21 - 3*kr)*sqrt(n/(1+(9*(kr^2))/2)), mean = 0, sd = 1) +
                         pnorm(q = -(zp22 - 3*kr)*sqrt(n/(1+(9*(kr^2))/2)), mean = 0, sd = 1)
      , LP2    = N2/D2

      , ASN2    = n/D2
        ) %>%

      dplyr::filter(LP1 >= (1-.alpha), LP2 <= .beta, ka <= 1) %>%
      dplyr::filter(ASN2 == min(ASN2)) %>%
      dplyr::select(p1, p2, n, ka,  kr, ASN = ASN2)
    Out0 <- rbind(Out0, Out1)
    }
  Out <-
  Out0 %>%
    dplyr::filter(ASN == min(ASN))
  return(Out)
}

