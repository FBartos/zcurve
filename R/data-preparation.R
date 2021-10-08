#' @title Prepare data for z-curve
#' 
#' @description \code{zcurve_data} is used to prepare data for the 
#' [zcurve()] function. The function transform strings containing 
#' reported test statistics \code{"z", "t", "f", "chi", "p"} into two-sided 
#' p-values. Test statistics reported as inequalities are as considered 
#' to be censored as well as test statistics reported with low accuracy 
#' (i.e., rounded to too few decimals). See details for more information.
#' 
#' @param data a vector strings containing the test statistics.
#' @param rounded an optional argument specifying rounding to be applied. Defaults 
#' to \code{TRUE} which automatically extracts the number of decimals from input. 
#' Another option is \code{FALSE} to treat all input as exact values or a numeric 
#' vector with values specifying precision of the input (with 0 representing exact 
#' values).
#' @param stat_precise an integer specifying the numerical precision of 
#' \code{"z", "t", "f"} statistics treated as exact values.
#' @param p_precise an integer specifying the numerical precision of 
#' p-values treated as exact values.
#'
#' @details By default, the function extract the type of test statistic:
#' \itemize{
#'  \item{\code{"F(df1, df2)=x"}}{F-statistic with df1 and df2 degrees of freedom,}
#'  \item{\code{"chi(df)=x"}}{Chi-square statistic with df degrees of freedom,}
#'  \item{\code{"t(df)=x"}}{for t-statistic with df degrees of freedom,}
#'  \item{\code{"z=x"}}{for z-statistic,}
#'  \item{\code{"p=x"}}{for p-value.}
#' }
#' The input is not case sensitive and automatically removes empty spaces. Furthermore, 
#' inequalities (\code{"<"} and \code{">"}) can be used to denote censoring. I.e., that 
#' the p-value is lower than \code{"x"} or that the test statistic is larger than \code{"x"} 
#' respectively. The automatic de-rounding procedure (if \code{rounded = TRUE}) treats 
#' p-values with less decimal places than specified in \code{p_precise} or test statistics 
#' with less decimal places than specified in \code{stat_precise} as censored on an interval 
#' that could result in a given rounded value. I.e., a \code{"p = 0.03"} input would be 
#' de-rounded as a p-value lower than 0.035 but larger than 0.025. 
#'
#'
#' @return An object of type \code{"zcurve_data"}.
#' @export zcurve_data
#'
#' @examples
#' # Specify a character vector containing the test statistics
#' data <- c("z = 2.1", "t(34) = 2.21", "p < 0.03", "F(2,23) > 10", "p = 0.003")
#'
#' # Obtain the z-curve data object
#' data <- zcurve_data(data)
#' 
#' # inspect the resulting object
#' data
#' @seealso [zcurve()], [print.zcurve_data()], [head.zcurve_data()]
zcurve_data <- function(data, rounded = TRUE, stat_precise = 2, p_precise = 3){
  
  data <- tolower(data) 
  data <- gsub(" ", "", data) 
  
  # deal with chi^2
  data <- gsub("pchisq", "c", data)
  data <- gsub("chi", "c", data)
  
  # extract the values
  stat_type <- substr(data, 1, 1)
  stat_val  <- substr(data, regexpr("[=]|[<]|[>]", data) + 1, nchar(data))
  stat_df1  <- ifelse(stat_type %in% c("t", "f", "c"), substr(data, regexpr("\\(", data) + 1, regexpr("[,]|[\\)]", data) - 1), NA)
  stat_df2  <- ifelse(stat_type == "f",           substr(data, regexpr(",", data) + 1, regexpr("[\\)]", data) - 1), NA)
  censored  <- grepl("<", data) | grepl(">", data)
  digits    <- ifelse(regexpr("\\.", data) == -1, 0, nchar(data) - regexpr("\\.", data))
  
  # check the input
  if(any(!stat_type %in% c("t", "z", "p", "f", "c")))
    stop(paste0("Unknown test statistic: ", paste0("'", unique(stat_type[!stat_type %in% c("t", "z", "p", "f", "c")]),"'", collapse = ", "), "."))
  
  
  # check that all matches are numeric
  stat_val <- tryCatch(
    as.numeric(stat_val),
    warning = function(w) stop(paste0("The following input could not be decoded: ", paste0("'", data[which(is.na(suppressWarnings(as.numeric(stat_val))))], "'", collapse = ", "), "."), call. = FALSE)
  )
  stat_df1 <- tryCatch(
    as.numeric(stat_df1),
    warning = function(w) stop(paste0("The following input could not be decoded: ", paste0("'", data[which(is.na(suppressWarnings(as.numeric(stat_df1))))], "'", collapse = ", "), "."), call. = FALSE)
  )
  stat_df2 <- tryCatch(
    as.numeric(stat_df2),
    warning = function(w) stop(paste0("The following input could not be decoded: ", paste0("'", data[which(is.na(suppressWarnings(as.numeric(stat_df2))))], "'", collapse = ", "), "."), call. = FALSE)
  )
  
  # set rounding (0 = un-rounded due to automatic conversion)
  if(length(rounded) == 1 && !rounded){
    # deal with the values as precise values
    rounded <- rep(0, length(data))
  }else if(length(rounded) == 1 && rounded){
    # specify automatic rounding
    rounded <- rep(FALSE, length(data))
    rounded[stat_type == "p" & digits < p_precise]    <- digits[stat_type == "p" & digits < p_precise]
    rounded[stat_type != "p" & digits < stat_precise] <- digits[stat_type != "p" & digits < stat_precise]
  }else{
    # use user specify rounding
    if(length(rounded) != length(data))
      stop("The rounding indicator does not match the lenght of data input.")
    if(!is.numeric(rounded))
      stop("The rounding indicator is not numeric.")
    if(any(rounded < 0))
      stop("The rounding indicator must be non-negative.")
  }
  
  # prepare empty containers
  p_vals    <- rep(NA, length(data))
  p_vals.lb <- rep(NA, length(data))
  p_vals.ub <- rep(NA, length(data))
  
  # compute and allocate the p-values accordingly
  for(i in seq_along(data)){
    if(rounded[i] == 0 && !censored[i]){
      # precise non-censored values
      p_vals[i] <- tryCatch(
        switch(
          stat_type[i],
          "f" = stats::pf(stat_val[i], df1 = stat_df1[i], df2 = stat_df2[i], lower.tail = FALSE),
          "c" = stats::pchisq(stat_val[i], df = stat_df1[i], lower.tail = FALSE),
          "t" = stats::pt(abs(stat_val[i]), df = stat_df1[i], lower.tail = FALSE) * 2,
          "z" = stats::pnorm(abs(stat_val[i]), lower.tail = FALSE) * 2,
          "p" = stat_val[i]
        ),
        warning = function(w) stop(paste0("The following input could not be decoded: '", data[i], "'."))
      )
    }else if(rounded[i] == 0 && censored[i]){
      # precise censored values
      p_vals.ub[i] <- tryCatch(
        switch(
          stat_type[i],
          "f" = stats::pf(stat_val[i], df1 = stat_df1[i], df2 = stat_df2[i], lower.tail = FALSE),
          "c" = stats::pchisq(stat_val[i], df = stat_df1[i], lower.tail = FALSE),
          "t" = stats::pt(abs(stat_val[i]), df = stat_df1[i], lower.tail = FALSE) * 2,
          "z" = stats::pnorm(abs(stat_val[i]), lower.tail = FALSE) * 2,
          "p" = stat_val[i]
        ),
        warning = function(w) stop(paste0("The following input could not be decoded: '", data[i], "'."))
      )
      p_vals.lb[i] <- 0
    }else if(rounded[i] != 0 && !censored[i]){
      # rounded non-censored values
      p_vals.ub[i] <- tryCatch(
        switch(
          stat_type[i],
          "f" = stats::pf(stat_val[i] - 0.5 * 10^-digits[i] , df1 = stat_df1[i], df2 = stat_df2[i], lower.tail = FALSE),
          "c" = stats::pchisq(stat_val[i] - 0.5 * 10^-digits[i], df = stat_df1[i], lower.tail = FALSE),
          "t" = stats::pt(abs(stat_val[i]) - 0.5 * 10^-digits[i], df = stat_df1[i], lower.tail = FALSE) * 2,
          "z" = stats::pnorm(abs(stat_val[i]) - 0.5 * 10^-digits[i], lower.tail = FALSE) * 2,
          "p" = stat_val[i] + 0.5 * 10^-digits[i]
        ),
        warning = function(w) stop(paste0("The following input could not be decoded: '", data[i], "'."))
      )
      p_vals.lb[i] <- tryCatch(
        switch(
          stat_type[i],
          "f" = stats::pf(stat_val[i] + 0.5 * 10^-digits[i] , df1 = stat_df1[i], df2 = stat_df2[i], lower.tail = FALSE),
          "c" = stats::pchisq(stat_val[i] + 0.5 * 10^-digits[i], df = stat_df1[i], lower.tail = FALSE),
          "t" = stats::pt(abs(stat_val[i]) + 0.5 * 10^-digits[i], df = stat_df1[i], lower.tail = FALSE) * 2,
          "z" = stats::pnorm(abs(stat_val[i]) + 0.5 * 10^-digits[i], lower.tail = FALSE) * 2,
          "p" = stat_val[i] - 0.5 * 10^-digits[i]
        ),
        warning = function(w) stop(paste0("The following input could not be decoded: '", data[i], "'."))
      )
    }else if(rounded[i] != 0 && !censored[i]){
      # rounded censored values
      p_vals.ub[i] <- tryCatch(
        switch(
          stat_type[i],
          "f" = stats::pf(stat_val[i] - 0.5 * 10^-digits[i] , df1 = stat_df1[i], df2 = stat_df2[i], lower.tail = FALSE),
          "c" = stats::pchisq(stat_val[i] - 0.5 * 10^-digits[i], df = stat_df1[i], lower.tail = FALSE),
          "t" = stats::pt(abs(stat_val[i]) - 0.5 * 10^-digits[i], df = stat_df1[i], lower.tail = FALSE) * 2,
          "z" = stats::pnorm(abs(stat_val[i]) - 0.5 * 10^-digits[i], lower.tail = FALSE) * 2,
          "p" = stat_val[i] + 0.5 * 10^-digits[i]
        ),
        warning = function(w) stop(paste0("The following input could not be decoded: '", data[i], "'."))
      )
      p_vals.lb[i] <- 0
    }
  }
  
  output <- list(
    precise  = data.frame(
      "input" = data[!is.na(p_vals)],
      "p"     = p_vals[!is.na(p_vals)]
    ),
    censored = data.frame(
      "input" = data[!is.na(p_vals.lb)],
      "p.lb"  = p_vals.lb[!is.na(p_vals.lb)],
      "p.ub"  = p_vals.ub[!is.na(p_vals.ub)]
    )
  )
  class(output) <- "zcurve_data"
  
  return(output)
}

### methods
#' Prints a z-curve data object
#' @param x z-curve data object
#' @param ... Additional arguments
#' @export  print.zcurve_data
#' @rawNamespace S3method(print, zcurve_data)
#' @seealso [zcurve_data()]
print.zcurve_data <- function(x, ...){
  cat(paste0("Object of class z-curve data with ", nrow(x$precise), " precise and ", nrow(x$censored), " censored p-values.\n\n"))
  cat("Precise p-values:\n")
  print(x$precise, ...)
  cat("\n")
  cat("Censored p-values:\n")
  print(x$censored, ...)
}

#' Prints first few rows of a z-curve data object
#' @param x z-curve data object
#' @param ... Additional arguments
#' @export  head.zcurve_data
#' @rawNamespace S3method(head, zcurve_data)
#' @seealso [zcurve_data()]
head.zcurve_data <- function(x, ...){
  cat(paste0("Object of class z-curve data with ", nrow(x$precise), " precise and ", nrow(x$censored), " censored p-values.\n\n"))
  cat("Precise p-values:\n")
  print(head(x$precise, ...))
  cat("\n")
  cat("Censored p-values:\n")
  print(head(x$censored, ...))
}
