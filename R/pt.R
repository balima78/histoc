#' Points for age differences
#'
#' @description Punctuation given for age difference between candidates and donors
#' @param dage A numeric value with donor's age.
#' @param cage A numeric value with candidate's age.
#' @param pts A numerical value for the points to age difference
#' @return A numerical value for pre-defined points
#' @examples
#' pts_age(dage = 60, cage = 40, pts = 4)
#' @export
pts_age <- function(dage = 60
                    , cage = 40
                    , pts = 4){
  # verify function parameters
  if(!is.numeric(dage) | dage < 18 | dage > 99){
    stop("donor's age is not valid!\n")}
  if(!is.numeric(cage) | cage < 18 | cage > 99){
    stop("candidate's age is not valid!\n")}
  if(!is.numeric(pts) | pts < 1 | pts > 20){
    stop("age points are not valid!\n")}

  pts<-ifelse((dage > 60 & cage < 55) | (dage < 40 & cage > 55), 0, pts)

  return(pts)

}

#' Points for cPRA sensitization
#'
#' @description Punctuation given for sensitized patients according to cPRA value
#' @param cPRA A numerical value (0-100) for cPRA
#' @param pts.80 A numerical value for the points to a cPRA >= 80
#' @param pts.50 A numerical value for the points to a cPRA >= 50
#' @return A numerical value for pre-defined points
#' @examples
#' pts_PRA(cPRA = 0, pts.80 = 8, pts.50 = 4)
#' @export
pts_PRA <- function(cPRA = 0
                    , pts.80 = 8
                    , pts.50 = 4){
  # verify function parameters
  if(!is.numeric(cPRA) | cPRA < 0 | cPRA > 100){
    stop("PRA value is not valid!\n")}
  if(!is.numeric(pts.80) | pts.80 < 0 | pts.80 > 100){
    stop("attributed points for a PRA >= 80% is not valid!\n")}
  if(!is.numeric(pts.50) | pts.50 < 0 | pts.50 > 100){
    stop("attributed points for a PRA >= 50% is not valid!\n")}

  pts<-if_else(cPRA >= 80, pts.80,
               if_else(cPRA >= 50, pts.50, 0))

  return(pts)

}

#' Points for HLA mismatches
#'
#' @description Punctuation given according to HLA mismatchs (mm) for item A) to E) from PT's algorithm
#' @param itemA  Points for HLA fullmatch (no mm for HLA-A, B and DR)
#' @param itemB  Points without mm for HLA-B and DR
#' @param itemC  Points with 1 mm for HLA-B and DR
#' @param itemD  Points with 1 mm for HLA-B and 1 mm for DR
#' @param itemE  Points for remaing possibilities
#' @return A numerical value for pre-defined points
#' @examples
#' pts_HLA(itemA = 12, itemB = 8, itemC = 4, itemD = 2, itemE = 1,
# dA = c('1','2'), dB = c('5','7'), dDR = c('1','4'),
# cA = c('1','2'), cB = c('3','15'), cDR = c('4','7'))
#' @export
pts_HLA <- function(itemA = 12
                    , itemB = 8
                    , itemC = 4
                    , itemD = 2
                    , itemE = 1
                    , mm.A = 0
                    , mm.B = 0
                    , mm.DR = 0){
  # verify function parameters
  if(!is.numeric(itemA) | itemA < 0 | itemA > 99){
    stop("points for 0 mmHLA (full match) is not valid!\n")}
  if(!is.numeric(itemB) | itemB < 0 | itemB > 99){
    stop("points for 0 mmB and mmDR is not valid!\n")}
  if(!is.numeric(itemC) | itemC < 0 | itemC > 99){
    stop("points for 1 mmB or mmDR is not valid!\n")}
  if(!is.numeric(itemD) | itemD < 0 | itemD > 99){
    stop("points for 1 mmB and 1 mmDR is not valid!\n")}
  if(!is.numeric(itemE) | itemE < 0 | itemE > 99){
    stop("points for more than 2 mmB and mmDR is not valid!\n")}

  mm <- list('mmA' = mm.A,
             'mmB' = mm.B,
             'mmDR' = mm.DR,
             'mmHLA' = mm.A + mm.B + mm.DR)

  pts<-if_else(mm[["mmHLA"]] == 0, itemA,
               if_else(mm[["mmB"]]+mm[["mmDR"]] == 0, itemB,
                       if_else(mm[["mmB"]]+mm[["mmDR"]] == 1, itemC,
                               if_else(mm[["mmB"]] == 1 & mm[["mmDR"]] == 1, itemD,
                                       itemE))))
  return(pts)
}


#' Matching punctuation' according to 2007 PT's algorithm
#'
#' @description Ordering of waitlisted candidates for a given donor and according to PT's algorithm.
#' @param iso A logical value for isogroupal compatibility.
#' @param dABO A character value with ABO blood group.
#' @param dA donor's HLA-A typing.
#' @param dB donor's HLA-B typing.
#' @param dDR donor's HLA-DR typing.
#' @param dage A numeric value with donor's age.
#' @param data A data frame containing demographics and medical information for a group of waitlisted transplant candidates with color priority classification.
#' @param df.abs A data frame with candidates' antibodies.
#' @param pts.80 A numerical value for the points to a cPRA >= 80
#' @param pts.50 A numerical value for the points to a cPRA >= 50
#' @param pts.dial
#' @param pts.age A numerical value for the points to age difference
#' @param n A positive integer to slice the first candidates.
#' @return An ordered data frame with a column 'cp' (color priority), 'sp', 'hi' and 'mmHLA'.
#' @examples
#' pt1(iso = TRUE, dABO = "A", dA = c("1","2"), dB = c("15","44"), dDR = c("1","4"),
# dage = 65,  data = cp(candidates),  df.abs = abs, n = 2)
#' @export
pt1 <- function(iso = TRUE
                , dABO = "O"
                , dA = c("1","2"), dB = c("15","44"), dDR = c("1","4")
                , dage = 65
                , df.abs = cabs
                , data = candidates
                , pts.80 = 8
                , pts.50 = 4
                , pts.dial = 0.1
                , pts.age = 4
                , n = 2
                , ...){

  n <- max(1, n)

  merge(data,
        xmatch(dA = dA, dB = dB, dDR = dDR, df.abs = df.abs),
        all.x=TRUE) %>%
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
           mmHLA = mmA + mmB + mmDR,
           ptsHLA = pts_HLA(mm.A = mmA, mm.B = mmB, mm.DR = mmDR, ...),
           ptsPRA = pts_PRA(cPRA = cPRA, pts.80 = pts.80, pts.50 = pts.50),
           ptsage = pts_age(dage = dage, cage = age, pts = pts.age),
           ptsdial = pts.dial * dialysis,
           ptsPT = ptsHLA + ptsPRA + ptsage + ptsdial
           ) %>%
    ungroup() %>%
    # only candidates ABO and HLA compatibles and those in the same group age with the donor are selectec
    filter(compBlood == TRUE & (xm == FALSE | is.na(xm))) %>%
    # order by Lima's algorithm
    arrange(HI, desc(ptsPT)) %>%
    # keep only the first n
    slice(1:n) %>%
    select(ID, bg,
           A1, A2, B1, B2, DR1, DR2,
           mmA, mmB, mmDR, mmHLA,
           age, donor_age, dialysis, cPRA, HI,
           ptsPT, SP,
           ptsHLA, ptsPRA, ptsage, ptsdial)

}
