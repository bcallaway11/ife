% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup_cov_pre_test.R
\name{setup_cov_pre_test}
\alias{setup_cov_pre_test}
\title{Setup data for covariate pre-test}
\usage{
setup_cov_pre_test(
  yname,
  gname,
  tname,
  idname,
  data,
  nife = 1,
  anticipation = 0,
  ...
)
}
\arguments{
\item{yname}{Name of outcome in \code{data}}

\item{gname}{Name of group in \code{data}}

\item{tname}{Name of time period in \code{data}}

\item{idname}{Name of id in \code{data}}

\item{data}{balanced panel data}

\item{nife}{the number of interactive fixed effects to include in the model}

\item{anticipation}{Number of periods that treatment is anticipated.  Default
is 0.  This is in ``periods''; e.g., code will work in time periods are
equally spaced but 2 years apart.  In this case, to allow for treatment
anticipation of 2 year (<=> 1 period), set \code{anticipation = 1}.}

\item{...}{extra arguments}
}
\value{
\code{pte_params} object
}
\description{
Sets up the data for conducting covariate pre-test
}
