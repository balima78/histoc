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
    mutate(cp = ifelse(urgent == 1, 1,
                            ifelse(cPRA >= cPRA2 | dialysis >= q3, 2,
                                   ifelse(cPRA >= cPRA1 | dialysis >= q2, 3, 4)
                                   )
                       )
           )

  data$cp <- factor(data$cp, levels = 1:4, labels = c('Red', 'Orange', 'Yellow', 'Green'))

  return(data)
  }
