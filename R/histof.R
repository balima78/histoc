#' test for ABO compatibility
#'
#' @description ABO compatibility test between donor and candidate
#' @param cABO A character from 'A', 'B', 'AB', 'O'
#' @param dABO A character from 'A', 'B', 'AB', 'O'
#' @param iso A logical value T/F
#' @return A logical value T/F
#' @examples
#' abo(cABO = 'A', dABO = 'A', iso = TRUE)
#' @export
abo <- function(cABO, dABO, iso){

  if(iso == TRUE){
    value <- cABO == dABO
  } else {
    value <- ifelse(dABO == 'O', TRUE,
                    ifelse(dABO == 'A' & cABO %in% c('A','AB'),TRUE,
                           ifelse(dABO == 'B' & cABO %in% c('B','AB'),TRUE,
                                  ifelse(dABO == 'AB' & cABO == 'AB', TRUE, FALSE)
                                  )
                           )
                    )
  }

  return(value)
}

#' number of HLA mismatchs
#'
#' @description Computes the number of HLA mismatchs between donor and candidate
#' @param dA donor's HLA-A typing
#' @param dB donor's HLA-B typing
#' @param dDR donor's HLA-DR typing
#' @param cA candidate's HLA-A typing
#' @param cB candidate's HLA-B typing
#' @param cDR candidate's HLA-DR typing
#' @return mmA number of HLA-A mismatchs between \code{dA} and \code{cA};
#' mmB number of HLA-B mismatchs between \code{dB} and \code{cB};
#' mmDR number of HLA-DR mismatchs between \code{dA}DRand \code{cDR};
#' and mmHLA as the sum of mmA + mmB + mmDR
#' @examples
#' mmHLA(dA = c('01','02'), dB = c('05','07'), dDR = c('01','04'), cA = c('01','02'), cB = c('03','15'), cDR = c('04','07'))
#' @export
mmHLA <- function(dA, dB, dDR, cA, cB, cDR){

  # mismatchs HLA-A
  mmA <- ifelse(dA[1] %in% cA & dA[2] %in% cA, 0,
                ifelse(dA[1] %in% cA | dA[2] %in% cA, 1, 2) )
  # mismatchs HLA-B
  mmB <- ifelse(dB[1] %in% cB & dB[2] %in% cB, 0,
                ifelse(dB[1] %in% cB | dB[2] %in% cB, 1, 2) )
  # mismatchs HLA-DR
  mmDR <- ifelse(dDR[1] %in% cDR & dDR[2] %in% cDR, 0,
                ifelse(dDR[1] %in% cDR | dDR[2] %in% cDR, 1, 2) )

  mmHLA <- mmA + mmB + mmDR

  list(mmHLA = mmHLA,
       mmA = mmA,
       mmB = mmB,
       mmDR = mmDR)
}

#' virtual crossmatch (XM)
#'
#' @description returns candidates' virtual crossmatch againts donor's HLA typing
#' @param dA donor's HLA-A typing
#' @param dB donor's HLA-B typing
#' @param dDR donor's HLA-DR typing
#' @param df.abs data frame with candidates' antibodies
#' @return A dataframe with candidates' ID and XM result POS/NEG
#' @examples
#' xmatch(dA = c('01','02'), dB = c('05','07'), dDR = c('01','04'), df.abs = abs)
#' @export
xmatch <- function(dA, dB, dDR, df.abs){

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package \"dplyr\" must be installed to use this function.",
      call. = FALSE
    )
  }

  dhla <- c(paste0('A',dA[1]),
            paste0('A',dA[2]),
            paste0('B',dB[1]),
            paste0('B',dB[2]),
            paste0('DR',dDR[1]),
            paste0('DR',dDR[2]))

  df.abs %>%
    dplyr::mutate(res = abs %in% dhla) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(xm = ifelse(sum(res)>0, "POS","NEG"))
}

#' Hiperimunized classification
#'
#' @description returns candidates' hiperimunized classification according to a cutoff value
#' @param cPRA candidate's cPRA value
#' @param cutoff A value to compare candidate's cPRA
#' @return A logical value T/F when cPRA >= cutoff
#' @examples
#' hiper(cPRA = 99, cutoff = 85)
#' @export
hiper <- function(cPRA, cutoff){

  value <- cPRA >= cutoff
  return(value)
}

#' Senior Program classification
#'
#' @description Returns 1 when candidates' belongs to Senior Program. Prioritization for older patients for older donors, followed for young patients for younger donors, and for last the remaining patients.
#' @param dage A numeric value with donor's age
#' @param cage A numeric value with candidate's age
#' @return The value 1 for a candidates older than 65 with also a donor older than 65
#' @examples
#' sp(dage = 66, cage = 70)
#' @export
sp <- function(dage, cage){

  value <- ifelse(dage >= 65 & cage >= 65, 1,
                  ifelse(dage < 65 & cage < 65, 2,3))
  return(value)
}
