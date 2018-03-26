#' Function to normalize distributions
#'
#'@param x a vector.
#'@param plotit Boolean.
#'
#'@examples
#' library(letsR)
#' x <- as.vector(na.exclude(values(temp)))
#' x <- x / 100
#' x_t <- normalizer(x, plotit = TRUE)

normalizer <- function(x, plotit = FALSE) {
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
    trans[[9]] <- rcompanion::transformTukey(x1, plotit = FALSE)
  }

  qual <- which.min(sapply(trans, .fitD))

  transformations <- c("raw",
                       "sqrt",
                       "log",
                       "cube root",
                       "boxcox",
                       "inverse(1/x)",
                       "log10",
                       "inverse log",
                       "Tukey")
  print(transformations[qual])

  if (plotit) {
    par(mfrow = c(3, 3))
    for (i in 1:length(trans)) {
      hist(trans[[i]], 30,
           main = transformations[i],
           xlab = "")
    }
  }


  return(trans[[qual]])
}
.fitD <- function(z) {
  MASS::fitdistr(z, "normal")$loglik
}

