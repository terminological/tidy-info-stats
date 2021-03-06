% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyDiscreteEntropy.R
\name{calculateDiscreteEntropy_Grassberger}
\alias{calculateDiscreteEntropy_Grassberger}
\title{calculate entropy of an optionally discrete value (X) using a histogram approach using the following method}
\usage{
calculateDiscreteEntropy_Grassberger(df, groupVars, countVar = NULL, ...)
}
\arguments{
\item{df}{- may be grouped, in which case the grouping is interpreted as different types of discrete variable}

\item{groupVars}{- the columns of the discrete value quoted by the vars() function (e.g. ggplot facet_wrap)}
}
\value{
a dataframe containing the disctinct values of the groups of df, and for each group an entropy value (H). If df was not grouped this will be a single entry
}
\description{
P. Grassberger, “Entropy Estimates from Insufficient Samplings,” arXiv [physics.data-an], 29-Jul-2003 [Online]. Available: http://arxiv.org/abs/physics/0307138
}
\details{
but with a digamma based function (rather than harmonics) detailed in eqns 31 & 35.
For our purposes we fix l=0 to give the form in eqn 27. The error in this method is supposedly better for undersampled cases (where number of bins similar to number of samples)
}
