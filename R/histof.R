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
abo <- function(cABO = 'A', dABO = 'A', iso = TRUE){

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
#' mmHLA(dA = c('1','2'), dB = c('5','7'), dDR = c('1','4'), cA = c('1','2'), cB = c('03','15'), cDR = c('04','07'))
#' @export
mmHLA <- function(dA = c('1','2'), dB = c('5','7'), dDR = c('1','4'),
                  cA = c('1','2'), cB = c('3','15'), cDR = c('4','7')){

  mmA = NULL
  mmB = NULL
  mmDR = NULL

  # verify function parameters
  if(!is.character(dA)){stop("donor's HLA-A typing is not valid!\n")}
  if(!is.character(dB)){stop("donor's HLA-B typing is not valid!\n")}
  if(!is.character(dDR)){stop("donor's HLA-DR typing is not valid!\n")}
  if(!is.character(cA)){stop("candidate's HLA-A typing is not valid!\n")}
  if(!is.character(cB)){stop("candidate's HLA-B typing is not valid!\n")}
  if(!is.character(cDR)){stop("candidate's HLA-DR typing is not valid!\n")}

  # compute missmatches
  mmA<-if_else((dA[1] %in% cA & dA[2] %in% cA) | (dA[1] %in% cA & (is.na(dA[2]) | dA[2] == "")), 0,
               if_else(dA[1] %in% cA | dA[2] %in% cA, 1,
                       if_else(!dA[1] %in% cA & (is.na(dA[2]) | dA[2] == ""), 1,
                               if_else(dA[1] == dA[2], 1,2))))

  mmB<-if_else((dB[1] %in% cB & dB[2] %in% cB) | (dB[1] %in% cB & (is.na(dB[2]) | dB[2] == "")), 0,
               if_else(dB[1] %in% cB | dB[2] %in% cB, 1,
                       if_else(!dB[1] %in% cB & (is.na(dB[2]) | dB[2] == ""), 1,
                               if_else(dB[1] == dB[2], 1,2))))

  mmDR<-if_else((dDR[1] %in% cDR & dDR[2] %in% cDR) | (dDR[1] %in% cDR & (is.na(dDR[2]) | dDR[2] == "")), 0,
                if_else(dDR[1] %in% cDR | dDR[2] %in% cDR, 1,
                        if_else(!dDR[1] %in% cDR & (is.na(dDR[2]) | dDR[2] == ""), 1,
                                if_else(dDR[1] == dDR[2],1,2))))

  # resume results
  mmHLA = mmA + mmB + mmDR
  mm = c(mmA,mmB,mmDR,mmHLA)
  names(mm) <- c("mmA","mmB","mmDR","mmHLA")

  return(mm)
}

#' virtual crossmatch (XM)
#'
#' @description returns candidates' virtual crossmatch againts donor's HLA typing
#' @param dA donor's HLA-A typing
#' @param dB donor's HLA-B typing
#' @param dDR donor's HLA-DR typing
#' @param df.abs data frame with candidates' antibodies
#' @return A dataframe with candidates' ID and xm result POS/NEG
#' @examples
#' xmatch(dA = c('1','2'), dB = c('5','7'), dDR = c('1','4'), df.abs = cabs)
#' @export
xmatch <- function(dA = c('1','2'),
                   dB = c('5','7'),
                   dDR = c('1','4'),
                   df.abs = cabs){

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
hiper <- function(cPRA, cutoff = 85){

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

#' TRANSPLANTSCORE (Tx Score)
#'
#' @description Returns the estimated 5-year event (mortality or graft failure combined outcome) probability as described by Molnar, el al (2017).
#' @param ageR A numeric value with recipient's age
#' @param race A character value with recipient's race from the options: 'White', 'Black', 'Hispanic', 'Other'
#' @param causeESRD A numeric value with recipient's cause of End-Stage Renal Disease, with options: 'Other', 'Diabetes', 'Hypertension', 'Glomerulonephritis', 'Cystic disease'
#' @param timeD A numeric value with recipient's time on dialysis (months)
#' @param diabetesR A logical value with recipient's diabetic status
#' @param coronary A logical value with recipient's coronary artery disease status
#' @param albumin A numeric value with recipient's albumin (g/dL)
#' @param hemoglobin A numeric value with recipient's hemoglobin (g/dL)
#' @param ageD A numeric value with donor's age
#' @param diabetesD A logical value with donor's diabetic status, with options: 'Absence', 'Presence', 'Unknown'
#' @param ECD A logical value regarding Extended Criteria Donor
#' @param mmHLA_A A numeric value (0, 1, 2) with the number of HLA-A mismatchs
#' @param mmHLA_B A numeric value (0, 1, 2) with the number of HLA-B mismatchs
#' @param mmHLA_DR A numeric value (0, 1, 2) with the number of HLA-DR mismatchs
#' @return 5 year probability for combined outcome of mortality or graft failure
#' @examples
#' txscore(ageR = 20, race = "White", causeESRD = "Other", timeD = 12, diabetesR = F, coronary = F, albumin = 1.5, hemoglobin = 10, ageD = 30, diabetesD = "Absence", ECD = F, mmHLA_A = 0, mmHLA_B = 0, mmHLA_DR = 0)
#' @source \url{https://balima.shinyapps.io/scoreTx/}
#' @export
txscore <- function(ageR = 20
                    , race = "White"
                    #, insurance = 0
                    , causeESRD = "Other"
                    , timeD = 12 #
                    , diabetesR = F
                    , coronary = F
                    , albumin = 1.5
                    , hemoglobin = 10
                    , ageD = 30
                    , diabetesD = "Absence"
                    , ECD = F
                    #, mmHLA = "0"
                    , mmHLA_A = 0
                    , mmHLA_B = 0
                    , mmHLA_DR = 0
){

  mmHLA_ <- as.numeric(mmHLA_A) + as.numeric(mmHLA_B) + as.numeric(mmHLA_DR)
  mmHLA <- ifelse(mmHLA_ == 0 , '0',
                  ifelse(mmHLA_ < 4, '1-3', '4-6'))

  ageR <- ifelse(ageR < 35 , 0.0993,
                 ifelse(ageR <50 , -0.0784,
                        ifelse(ageR < 65, 0, 0.1881)))
  race <- ifelse(race == "White", 0,
                 ifelse(race == "Black", 0.1609,
                        ifelse(race == "Hispanic", -0.2554, -0.4475)))
  causeESRD <- ifelse(causeESRD == "Diabetes", 0,
                      ifelse(causeESRD == "Hypertension", 0.1541,
                             ifelse(causeESRD == "Glomerulonephritis", 0.1447,
                                    ifelse(causeESRD == "Cystic Disease", -0.1870, 0.3209))))
  timeD <- ifelse(timeD < 12, 0,
                  ifelse(timeD < 36, -0.2618,
                         ifelse(timeD < 61, -0.3747, -0.1432)))
  diabetesR <- ifelse(diabetesR == T, 0.3021, 0)
  coronary <- ifelse(coronary == T, 0.2617, 0)
  albumin <- (albumin - 4)*(-0.2644)
  hemoglobin <- (hemoglobin - 12.3)*(-0.0451)
  ageD <- (ageD - 39)*0.0059
  diabetesD <- ifelse(diabetesD == "Absence", 0,
                      ifelse(diabetesD == "Presence", 0.4596, -0.3308))
  ECD <- ifelse(ECD == T, 0.2082, 0)
  mmHLA <- ifelse(mmHLA == "0" , 0,
                  ifelse(mmHLA == "1-3", 0.3241, 0.3115))

  LP <- ageR + race + causeESRD + timeD + diabetesR + coronary + albumin + hemoglobin + ageD + diabetesD + ECD + mmHLA

  gamma <- 0.916

  PS = gamma * LP

  prob5y <- round((1-0.752292^exp(PS))*100,2)

  list(LP = LP
       , gamma = gamma
       , PS = PS
       , prob5y = prob5y)

}
