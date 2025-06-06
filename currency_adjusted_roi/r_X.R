r_X <- function(r_Y, e_XY, e_XY_prime) {
  R_XY <- (e_XY_prime / e_XY) - 1
  r_X <- (r_Y - R_XY) / (1 + R_XY)
  r_X
}
