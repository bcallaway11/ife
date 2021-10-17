#' @title Individual Panel Data about Job Displacement 
#'
#' @description The data used in Callaway and Karami (2021)
#'
#' @format A data frame with 18264 rows and 10 columns
#' \describe{
#'   \item{id}{individual identifier}
#'   \item{displaced}{whether an individual is ever displaced from their job in any time period in the sample}
#'   \item{EDUC}{education category}
#'   \item{first.displaced}{for displaced workers, the year that they were displaced; for non-displaced workers, it is equal to 0}
#'   \item{race}{race}
#'   \item{gender}{gender}
#'   \item{afqt}{armed forces qualification test normalized score}
#'   \item{year}{the year of the observation}
#'   \item{earn}{yearly earnings}
#'   \item{learn}{log of yearly earnings}
#' }
#' @source Callaway and Karami (2021)
"job_displacement_data"
