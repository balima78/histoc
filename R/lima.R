#' Candidates' Color Priority
#'
#' @description Classification of candidates according to waiting list time on dialysis' quartiles and two cPRA cutoff values.
#' @param data A data frame with information for candidates' waiting list.
#' @param q2 A numerical value for the median of candidates' waiting list.
#' @param q3 A numerical value for the third quartile of candidates' waiting list.
#' @param cPRA1 A numerical value (0-100) for the lower cPRA cutoff.
#' @param cPRA2 A numerical value (0-100) for the higher cPRA cutoff. cPRA2 must be greater than cPRA1.
#' @return A data frame with a new column 'cp' (color priority)
#' @examples
#' cp(data = candidates, q2 = 100, q3 = 60, cPRA1 = 50, cPRA2 = 85)
#' @export
cp <- function(data = candidates, q2 = 60, q3 = 100, cPRA1 = 50, cPRA2 = 85){

  # verify function parameters
  if(cPRA2 < cPRA1){
    stop("higher cPRA cutoff value (cPRA2) must be greater than lower cPRA cutoff (cPRA1)!\n")
  } else if (q2 >= q3){
    stop("median time on dialysis quartiles must be lower than third quartile: q2 < q3!\n")
    }

  data <- data %>%
    dplyr::mutate(cp = ifelse(urgent == 1, 1,
                            ifelse(cPRA >= cPRA2 | dialysis >= q3, 2,
                                   ifelse(cPRA >= cPRA1 | dialysis >= q2, 3, 4)
                                   )
                       )
           )

  data$cp <- factor(data$cp, levels = 1:4, labels = c('Red', 'Orange', 'Yellow', 'Green'))

  return(data)
  }


#' Candidates' ordering according to Lima's algorithm
#'
#' @description Ordering of waitlisted candidates for a given donor and according to to Lima's algorithm.
#' @param iso A logical value for isogroupal compatibility.
#' @param dABO A character value with ABO blood group.
#' @param dA donor's HLA-A typing.
#' @param dB donor's HLA-B typing.
#' @param dDR donor's HLA-DR typing.
#' @param dage A numeric value with donor's age.
#' @param data A data frame containing demographics and medical information for a group of waitlisted transplant candidates with color priority classification.
#' @param df.abs A data frame with candidates' antibodies.
#' @param n A positive integer to slice the first candidates.
#' @return An ordered data frame with a column 'cp' (color priority), 'sp', 'hi' and 'mmHLA'.
#' @examples
#' lima1(iso = TRUE, dABO = "A", dA = c("1","2"), dB = c("15","44"), dDR = c("1","4"), dage = 65,  data = cp(candidates),  df.abs = abs, n = dim(data)[1])
#' @export
lima1 <- function(iso = TRUE
                  , dABO = "O"
                  , dA = c("1","2"), dB = c("15","44"), dDR = c("1","4")
                  , dage = 60
                  , df.abs = cabs
                  , data = candidates
                  , n = 2
                  , ...){

  n <- max(1, n)

    merge(cp(data = data, ...),
          xmatch(dA = dA, dB = dB, dDR = dDR, df.abs = df.abs),
          all.x=TRUE) %>%
    # cp(data = data) %>%
    # left_join(
    #   xmatch(dA = dA, dB = dB, dDR = dDR, df.abs = df.abs)
    # ) %>%
    rowwise() %>%
    mutate(donor_age = dage,
           SP = sp(cage = age, dage = dage),
           HI = hiper(cPRA = cPRA),
           compBlood = abo(iso = iso, dABO = dABO, cABO = bg),
           mmA = mmHLA(dA = dA, dB = dB, dDR = dDR,
                       cA = c(A1,A2), cB = c(B1,B2), cDR = c(DR1,DR2))[["mmA"]],
           mmB = mmHLA(dA = dA, dB = dB, dDR = dDR,
                       cA = c(A1,A2), cB = c(B1,B2), cDR = c(DR1,DR2))[["mmB"]],
           mmDR = mmHLA(dA = dA, dB = dB, dDR = dDR,
                        cA = c(A1,A2), cB = c(B1,B2), cDR = c(DR1,DR2))[["mmDR"]],
           mmHLA = mmA + mmB + mmDR) %>%
    ungroup() %>%
    # only candidates ABO and HLA compatibles and those in the same group age with the donor are selectec
    filter(compBlood == TRUE & (xm == FALSE | is.na(xm)) & SP <3) %>%
    # order by Lima's algorithm
    arrange(desc(SP), cp, mmHLA, desc(dialysis)) %>%
    # keep only the first n
    slice(1:n) %>%
    select(ID, bg,
           A1, A2, B1, B2, DR1, DR2,
           mmA, mmB, mmDR, mmHLA,
           age, donor_age, dialysis, cPRA, HI,
           cp, SP)
}
