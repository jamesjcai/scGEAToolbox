%' the weightedZISpearman function calculates weighted rho\*,
%' where rho\* is described in Pimentel et al (2009).
%' This association measure is defined for zero-inflated,
%' non-negative random variables.
%'
%' @title weightedZISpearman
%' @param w weight vector, values should be between 0 and 1
%' @param x x and y are non-negative data vectors
%' @param y x and y are non-negative data vectors
%' @return \code{numeric} weighted rho* association value between x and y
%'
%'  Pimentel, Ronald Silva, "Kendall's Tau and Spearman's Rho for
%'  Zero-Inflated Data" (2009). Dissertations. 721.
%'  https://scholarworks.wmich.edu/dissertations/721

%' @examples
%'
%' x = pmax(0,rnorm(100))
%' y = pmax(0,rnorm(100))
%' w = runif(100)
%' weightedZISpearman(x,y,w)
%'
%' @export

weightedZISpearman <- function(x, y, w = 1) {

  % needs the original values, not the ranks

  if (any(x < 0 | y < 0)) {
    stop("x and/or y values have negative values")
  }
  if (length(x) != length(y)) {
    stop("x and y should have the same length")
  }
  if (length(w) == 1) {
    w <- rep(w, length(x))
  }

  posx = x > 0
  posy = y > 0
  pospos = posx & posy

  p_11 = sum(w * pospos)/sum(w)
  p_00 = sum(w * (!posx & !posy))/sum(w)
  p_01 = sum(w * (!posx & posy))/sum(w)
  p_10 = sum(w * (posx & !posy))/sum(w)
  rho_11 = weightedSpearman(x = x[pospos], y = y[pospos], w = w[pospos])
  rho_star = p_11 * (p_01 + p_11) * (p_10 + p_11) * rho_11 +
    3*(p_00 * p_11 - p_10 * p_01)

  if (is.na(rho_star)) {
    print("Zero inflated Spearman correlation is undefined,
          returning Spearman correlation")
    rho = weightedSpearman(x = x, y = y, w = w)
    return(rho)
  }

  return(rho_star)
}

