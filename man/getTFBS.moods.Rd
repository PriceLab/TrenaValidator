% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrenaValidator.R
\name{getTFBS.moods,TrenaValidator-method}
\alias{getTFBS.moods,TrenaValidator-method}
\alias{getTFBS.moods}
\title{return the TFBS}
\usage{
\S4method{getTFBS.moods}{TrenaValidator}(
  obj,
  tbl.regions,
  match.threshold,
  conservation.threshold,
  motifs
)
}
\arguments{
\item{obj}{An instance of the TrenaValidator class}

\item{tbl.regions}{a data.frame with chrom, start, end columns}

\item{match.threshold}{numeric, 1e-4 for example}

\item{conservation.threshold}{numeric, between zero and one}

\item{motifs}{a list of motifs, the result of a MotifDb query}
}
\value{
a data.frame
}
\description{
return the TFBS
}
