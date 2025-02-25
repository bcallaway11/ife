% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/staggered_ife_attgt.R
\name{staggered_ife_attgt}
\alias{staggered_ife_attgt}
\title{staggered_ife_attgt}
\usage{
staggered_ife_attgt(
  gt_data,
  nife = 1,
  xformla = ~1,
  anticipation = 0,
  ret_ife_regs = FALSE,
  ...
)
}
\arguments{
\item{gt_data}{data frame that is local to the specific groups
and times for which we'll be computing a treatment efect
estimate.}

\item{nife}{The number of interactive fixed effects.  Default is 1.}

\item{xformla}{Formula for which covariates to include in the model.  Default is ~1.}

\item{anticipation}{Number of periods that treatment is anticipated.  Default
is 0.  This is in ``periods''; e.g., code will work in time periods are
equally spaced but 2 years apart.  In this case, to allow for treatment
anticipation of 2 year (<=> 1 period), set \code{anticipation = 1}.}

\item{ret_ife_regs}{Whether or not to return the first stage ife regressions; default is FALSE.}

\item{...}{extra arguments}
}
\value{
\code{ptetools::attgt_if} object
}
\description{
Computes estimates of group-time average treatment
effects in an interactive treatment effects model for untreated
potential outcomes by exploiting staggered treatment adoption as
in Callaway and Tsyawo (2023).  This function is based on the local-
differencing approach (similar to the approach proposed in
Callaway and Karami).  See \code{staggered_ife_attgt2} for the main approach
discussed in the paper where all pre-treatment periods are used
to estimate the interactive fixed effects model.
}
