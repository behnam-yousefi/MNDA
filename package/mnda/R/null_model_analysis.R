## Given a value (x) and a set of null values, one can calculate the p.value for the extremeness of x
## using two methods:
## 1) rank-based
## 2) normal distribution-based
## 3) count-based

#' Calculate p.value for x given a set of null values using ranks
#'
#' @param x numeric value
#' @param null.values a numeric vector of null distribution samples
#' @param alternative alterative test including: \code{"two.sided"} [default],
#' \code{"greater"}, and  \code{"less"}
#'
#' @return p.value
#' @export
#'
#' @examples
#' p.val = p_val_rank(1, 1:100)
#'
p_val_rank = function(x, null.values, alternative = "two.sided"){

  y = c(x, null.values)

  if (alternative == "greater"){
    r = Rank(y, decreasing = TRUE)
    p.val = r[1]/length(y)
  }else if (alternative == "less"){
    r = Rank(y, decreasing = FALSE)
    p.val = r[1]/length(y)
  }else if (alternative == "two.sided"){
    r1 = Rank(y, decreasing = FALSE)
    r2 = Rank(y, decreasing = TRUE)
    p.val = min(c(r1[1],r2[1])) / (length(y)/2)
  }

  return(p.val)
}

#' Calculate p.value for x given set of null values using a Gaussian null pdf
#'
#' @param x numeric value
#' @param null.values a numeric vector of null distribution samples
#' @param alternative alterative test including: \code{"two.sided"} [default],
#' \code{"greater"}, and  \code{"less"}
#'
#' @return p.value
#' @export
#'
#' @examples
#' p.val = p_val_norm(1, rnorm(1000,0,1))
#'
p_val_norm = function(x, null.values, alternative = "two.sided"){

  mu = mean(null.values)
  std = stats::sd(null.values)

  if(alternative == "less")
    p.val = stats::pnorm(x, mu, std, lower.tail = TRUE)
  else if(alternative == "greater")
    p.val = stats::pnorm(x, mu, std, lower.tail = FALSE)
  else if(alternative == "two.sided")
    p.val = 2 * stats::pnorm(abs(x), mu, std, lower.tail = FALSE)

  return(p.val)
}

#' Calculate p.value for x given set of null values using counts
#'
#' @param x numeric value
#' @param null.values a numeric vector of null distribution samples
#' @param alternative alterative test including: \code{"two.sided"} [default],
#' \code{"greater"}, and  \code{"less"}
#'
#' @return p.value
#' @export
#'
#' @examples
#' p.val = p_val_count(1, 1:100)
#'
p_val_count = function(x, null.values, alternative = "two.sided"){

  y = c(x, null.values)

  if (alternative == "greater"){
    p.val = sum(null.values>x) / length(y)
  }else if (alternative == "less"){
    p.val = sum(null.values<x) / length(y)
  }else if (alternative == "two.sided"){
    p.val1 = sum(null.values>x) / length(y)
    p.val2 = sum(null.values<x) / length(y)
    p.val = p.val1 + p.val2
  }

  return(p.val)
}





