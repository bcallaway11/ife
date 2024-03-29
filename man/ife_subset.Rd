% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ife_subset.R, R/staggered_ife_subset.R
\name{ife_subset}
\alias{ife_subset}
\title{ife_subset}
\usage{
ife_subset(data, g, tp, nife, anticipation = 0, ...)

ife_subset(data, g, tp, nife, anticipation = 0, ...)
}
\arguments{
\item{data}{the full dataset}

\item{g}{the current group}

\item{tp}{the current time period}

\item{nife}{the number of interactive fixed effects}

\item{anticipation}{number of periods before participating
treatment that treatment can have an effect on outcomes}

\item{...}{extra arguments to get the subset correct}
}
\value{
list that contains correct subset of data, \code{n1}
number of observations
in this subset, and \code{disidx} a vector of the correct ids for this
subset.

list that contains correct subset of data, \code{n1}
number of observations
in this subset, and \code{disidx} a vector of the correct ids for this
subset.
}
\description{
A function for obtaining the correct subset in a staggered treatment
adoption setup with panel data.  This is the subset of not yet treated
observations in a particular time period and in the previous
\code{nife+1} periods where \code{nife} is the number of interactive fixed effects.

A function for obtaining the correct subset in a staggered treatment
adoption setup with panel data.  This is the subset of not yet treated
observations in a particular time period and in the previous
\code{nife+1} periods where \code{nife} is the number of interactive fixed effects.
}
