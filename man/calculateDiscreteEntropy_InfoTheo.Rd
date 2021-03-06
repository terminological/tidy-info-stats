% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyDiscreteEntropy.R
\name{calculateDiscreteEntropy_InfoTheo}
\alias{calculateDiscreteEntropy_InfoTheo}
\title{calculate entropy of an optionally discrete value (X) using a infotheo library}
\usage{
calculateDiscreteEntropy_InfoTheo(
  df,
  groupVars,
  infoTheoMethod = "mm",
  collect = FALSE,
  ...
)
}
\arguments{
\item{df}{- may be grouped, in which case the grouping is interpreted as different types of discrete variable}

\item{groupVars}{- the column of the discrete value (X)}

\item{infoTheoMethod}{- the method of the entropy estimator ("mm","emp","shrink","sg")}

\item{collect}{- if false (the default) this will fail on a dbplyr table as this is not supported in SQL}
}
\value{
a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
}
\description{
calculate entropy of an optionally discrete value (X) using a infotheo library
}
