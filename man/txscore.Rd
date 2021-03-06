% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/histof.R
\name{txscore}
\alias{txscore}
\title{TRANSPLANTSCORE (Tx Score)}
\source{
\url{https://balima.shinyapps.io/scoreTx/}
}
\usage{
txscore(
  ageR = 20,
  race = "White",
  causeESRD = "Other",
  timeD = 12,
  diabetesR = F,
  coronary = F,
  albumin = 1.5,
  hemoglobin = 10,
  ageD = 30,
  diabetesD = "Absence",
  ECD = F,
  mmHLA_A = 0,
  mmHLA_B = 0,
  mmHLA_DR = 0
)
}
\arguments{
\item{ageR}{A numeric value with recipient's age}

\item{race}{A character value with recipient's race from the options: 'White', 'Black', 'Hispanic', 'Other'}

\item{causeESRD}{A numeric value with recipient's cause of End-Stage Renal Disease, with options: 'Other', 'Diabetes', 'Hypertension', 'Glomerulonephritis', 'Cystic disease'}

\item{timeD}{A numeric value with recipient's time on dialysis (months)}

\item{diabetesR}{A logical value with recipient's diabetic status}

\item{coronary}{A logical value with recipient's coronary artery disease status}

\item{albumin}{A numeric value with recipient's albumin (g/dL)}

\item{hemoglobin}{A numeric value with recipient's hemoglobin (g/dL)}

\item{ageD}{A numeric value with donor's age}

\item{diabetesD}{A logical value with donor's diabetic status, with options: 'Absence', 'Presence', 'Unknown'}

\item{ECD}{A logical value regarding Extended Criteria Donor}

\item{mmHLA_A}{A numeric value (0, 1, 2) with the number of HLA-A mismatchs}

\item{mmHLA_B}{A numeric value (0, 1, 2) with the number of HLA-B mismatchs}

\item{mmHLA_DR}{A numeric value (0, 1, 2) with the number of HLA-DR mismatchs}
}
\value{
5 year probability for combined outcome of mortality or graft failure
}
\description{
Returns the estimated 5-year event (mortality or graft failure combined outcome) probability as described by Molnar, el al (2017).
}
\examples{
txscore(ageR = 20, race = "White", causeESRD = "Other", timeD = 12, diabetesR = F, coronary = F, albumin = 1.5, hemoglobin = 10, ageD = 30, diabetesD = "Absence", ECD = F, mmHLA_A = 0, mmHLA_B = 0, mmHLA_DR = 0)
}
