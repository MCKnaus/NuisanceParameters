#' Pennsylvania Re-employment Bonus Experiment
#'
#' This dataset contains individual-level data from the Pennsylvania Re-employment Bonus Experiment,
#' originally conducted in the 1980s by the U.S. Department of Labor (Bilias, 2000). The experiment aimed
#' to evaluate the incentive effects of cash bonuses on early reemployment among unemployment insurance (UI) claimants.
#'
#' UI claimants were randomly assigned either to a control group or to one of six treatment groups, which differed in
#' bonus amount, qualification period, and workshop offer:
#'
#' \tabular{llll}{
#'   \strong{Group} \tab \strong{Bonus Amount} \tab \strong{Qualification Period} \tab \strong{Workshop Offer} \cr
#'   Control \tab 0 \tab 0 \tab No \cr
#'   Treatment 1 \tab Low (3 x WBA) \tab Short (6 Weeks) \tab Yes \cr
#'   Treatment 2 \tab Low (3 x WBA) \tab Long (12 Weeks) \tab Yes \cr
#'   Treatment 3 \tab High (6 x WBA) \tab Short (6 Weeks) \tab Yes \cr
#'   Treatment 4 \tab High (6 x WBA) \tab Long (12 Weeks) \tab Yes \cr
#'   Treatment 5 \tab Initially High (6 x WBA) \tab Long (12 Weeks) \tab Yes \cr
#'   Treatment 6 \tab High (6 x WBA) \tab Long (12 Weeks) \tab No
#' }
#'
#' Claimants in the treatment groups were offered a cash bonus if they secured employment within a specified
#' qualification period and retained the job for a minimum duration. The control group was subject to standard UI rules.
#'
#' @format A data frame with \eqn{N} rows and \eqn{K} variables. Covariates include demographic characteristics,
#' employment history, location, and occupation:
#' \describe{
#'   \item{abdt}{Chronological time of enrollment in the experiment}
#'   \item{tg}{Treatment group (bonus amount - qualification period)}
#'   \item{inuidur1}{Duration of the first unemployment spell (weeks)}
#'   \item{female}{1 if female, 0 if male}
#'   \item{black}{1 if Black, 0 otherwise}
#'   \item{hispanic}{1 if Hispanic, 0 otherwise}
#'   \item{othrace}{1 if non-white, non-Black, non-Hispanic, 0 otherwise}
#'   \item{dep}{Number of dependents}
#'   \item{q1}{Dummy indicating enrollment in quarter 1}
#'   \item{q2}{Dummy indicating enrollment in quarter 2}
#'   \item{q3}{Dummy indicating enrollment in quarter 3}
#'   \item{q4}{Dummy indicating enrollment in quarter 4}
#'   \item{q5}{Dummy indicating enrollment in quarter 5}
#'   \item{q6}{Dummy indicating enrollment in quarter 6}
#'   \item{recall}{1 if claimant expected to be recalled to previous job, 0 otherwise}
#'   \item{agelt35}{1 if age < 35, 0 otherwise}
#'   \item{agegt54}{1 if age > 54, 0 otherwise}
#'   \item{durable}{1 if occupation in durable manufacturing, 0 otherwise}
#'   \item{nondurable}{1 if occupation in non-durable manufacturing, 0 otherwise}
#'   \item{lusd}{1 if filed in low unemployment short-duration area, 0 otherwise}
#'   \item{husd}{1 if filed in high unemployment short-duration area, 0 otherwise}
#'   \item{muld}{1 if filed in moderate unemployment long-duration area, 0 otherwise}
#' }
#'
#' @details
#' The dataset is frequently used in labor economics and program evaluation to study the effects of financial incentives
#' on unemployment duration. Notably, most covariates are categorical or dummy variables; the response variable is continuous.
#'
#' @source
#' Kaggle: \url{https://www.kaggle.com/code/victorchernozhukov/analyzing-rct-reemployment-experiment/notebook}
#' 
#' For further details on the data: Bilias, 2000. Readme file available at 
#' \url{http://qed.econ.queensu.ca/jae/2000-v15.6/bilias/readme.b.txt}
#' 
#' @examples
#' data(pa_reemployment)
#' summary(pa_reemployment)
#' table(pa_reemployment$tg)
"pa_reemployment"

