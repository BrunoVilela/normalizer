#' Function to normalize distributions
#'
#'@description The function will apply several transformations,
#'and indicate the transformation that best approximates a normal distribution.
#'
#'@param x a numeric vector without absent data.
#'@param plotit Logical If TRUE histograms of each transformation will be ploted.
#'@param return_all Logical. If TRUE it will return all values.
#'Default is FALSE, and only the selected vector is returned.
#'@param ... Other parameters passed to the \code{hist} function.
#'
#'@details
#'The function will apply eigth different transformations to the data:
#'\enumerate{
#'\item Square root
#'\item Natural logarithm
#'\item Cube root
#'\item Box Cox
#'\item Logarithm in the base 10
#'\item Inverse (1 / x)
#'\item Inverse Natural logarithm (revert the data, apply a natural logarithm
#'and then invert it back)
#'\item Tukey Ladder of Powers (if the vector is between 3 and 5000 long)
#'}
#'The \code{normalizer} function uses the standardized D value calculated by a
#' Kolmogorov-Smirnov test using the function \code{stats::ks.test} to determine
#' the transformation that most approximate a normal distribution.
#'
#'@author Bruno Vilela (email: \email{bvilela@wustl.edu})
#'
#'@return The function will return a list with transformed vector
#'(or the raw data) that is closer to a normal distribution and the D value.
#' If \code{return_all} equals TRUE, than a list with all transformations
#' is returned.
#'
#'@examples
#'# Example 1:
#' library(letsR)
#' x <- as.vector(na.exclude(values(temp)))
#' x <- x / 100
#' x_t <- normalizer(x, plotit = TRUE)
#'
#'# Example 2:
#' library(datasets)
#' data(mtcars)
#' mpg_t <- normalizer(mtcars$mpg, plotit = TRUE, return_all = TRUE)


normalizer <- function(x, plotit = FALSE,
                       return_all = FALSE,
                       ...) {

  if (!is.numeric(x) | !is.vector(x)) {
    stop("x is not a numeric vector")
  }
  if (any(is.na(x)) | any(is.infinite(x)) | any(is.nan(x))) {
    stop("x contains NA, NANs or infinite values")
  }

  if (!is.logical(plotit)) {
    stop("plotit is not logical")
  }
  if (any(x < 0)) {
    x1 <- x + abs(min(x)) + 1
  } else {
    x1 <- x + 1
  }
  trans <- list()
  trans[[1]] <- x
  trans[[2]] <- sqrt(x1)
  trans[[3]] <- log(abs(x1))
  trans[[4]] <- sign(x1) * abs(x1) ^ (1 / 3)# Cube root transformation
  # Boxcox
  Box <- MASS::boxcox(x1 ~ 1, plotit = FALSE)
  Cox <- data.frame(Box$x, Box$y)
  Cox2 <- Cox[with(Cox, order(-Cox$Box.y)),]
  lambda <- Cox2[1, "Box.x"]
  lambda <- ifelse(lambda == 0, 0.000001, lambda)
  trans[[5]] <- (x1 ^ lambda - 1) / lambda
  # Inverse
  trans[[6]] <- 1 / x1
  #
  trans[[7]] <- log10(x1)
  #
  inv <- x1 * -1
  inv2 <- (inv - min(inv)) + 1
  trans[[8]] <- log(inv2) * -1

  if (length(x1) < 5000 & length(x1) > 3) {
    f <- file()
    sink(file = f)
    trans[[9]] <- rcompanion::transformTukey(x1, plotit = FALSE)
    sink()
    close(f)
  }
  trans2 <- lapply(trans, function(x){(x - mean(x)) / sd(x)})
  likelihood <- sapply(trans2, .fitD)
  qual <- which.min(likelihood)

  transformations <- c("raw",
                       "sqrt",
                       "log",
                       "cube root",
                       "boxcox",
                       "inverse(1/x)",
                       "log10",
                       "inverse log",
                       "Tukey")

  if (plotit) {
    par(mfrow = c(3, 3))
    for (i in 1:length(trans)) {
      hist(trans[[i]],
           main = transformations[i],
           xlab = "", ...)
    }
  }
  names(trans) <- transformations[1:length(trans)]
  names(likelihood) <- transformations[1:length(trans)]
  if (return_all) {
    result <- list("Transformed_vectors" = trans, "D" = likelihood)
    return(result)
  } else {
    result <- list(trans[[qual]], likelihood[qual])
    names(result) <- c(transformations[[qual]], "D")
    return(result)
  }
}
.fitD <- function(z) {
  suppressWarnings(stats::ks.test(z, "pnorm")$statistic)
}

# .fitD <- function(z) {
#   MASS::fitdistr(z, "normal")$loglik
# }
#
