% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pt.R
\name{pt1}
\alias{pt1}
\title{Matching punctuation' according to 2007 PT's algorithm}
\usage{
pt1(
  iso = TRUE,
  dABO = "O",
  dA = c("1", "2"),
  dB = c("15", "44"),
  dDR = c("1", "4"),
  dage = 65,
  df.abs = cabs,
  data = candidates,
  pts.80 = 8,
  pts.50 = 4,
  pts.dial = 0.1,
  pts.age = 4,
  n = 2,
  ...
)
}
\arguments{
\item{iso}{A logical value for isogroupal compatibility.}

\item{dABO}{A character value with ABO blood group.}

\item{dA}{donor's HLA-A typing.}

\item{dB}{donor's HLA-B typing.}

\item{dDR}{donor's HLA-DR typing.}

\item{dage}{A numeric value with donor's age.}

\item{df.abs}{A data frame with candidates' antibodies.}

\item{data}{A data frame containing demographics and medical information for a group of waitlisted transplant candidates with color priority classification.}

\item{pts.80}{A numerical value for the points to a cPRA >= 80}

\item{pts.50}{A numerical value for the points to a cPRA >= 50}

\item{pts.age}{A numerical value for the points to age difference}

\item{n}{A positive integer to slice the first candidates.}
}
\value{
An ordered data frame with a column 'cp' (color priority), 'sp', 'hi' and 'mmHLA'.
}
\description{
Ordering of waitlisted candidates for a given donor and according to PT's algorithm.
}
\examples{
pt1(iso = TRUE, dABO = "A", dA = c("1","2"), dB = c("15","44"), dDR = c("1","4"),
}
