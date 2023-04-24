library(purrr)
library(bvpSolve)

# Lenght of the tube (cm):
L <- 6
# Velocities (cm/h):
v <- c(0.001, 0.0417, 0.2043, 0.2807, 0.3857, 0.4521, 0.509, 0.5315, 0.536, 0.5405)
# The corresponding colors:
colors <- c("purple","blue","blue4","cadetblue3","cadetblue","green4","green2",
            "yellow3","orange3","red")
# Diffusion constant (cm2/h):
D <- 0.2
# Maximal growth rate (/h):
r <- 0.42
# Initial guesses of of Phi = F/F_in:
phi0 <- c(0, 0, 0.1, 0.1, 0.1, 0.1, 0.15, 0.15, 0.88, 1.0073009)
# Number of discreteization points
N <- 70

phi_values <- function(v_val, phi0_val) {
  kappa <- 0.1 / (1 / v_val)
  lambda <- r * D / v_val^2
  sigma <- L * v_val / D
  s <- seq(0, sigma, le = N) # independent variable
  
  func <- function(s, Y, pars){
    with (as.list(Y), {
      dphi <- f
      df <- f + (lambda * phi * (1 - phi)) / (kappa + phi)
      return(list(c(dphi, df)))
    })
  }
  
  bound <- function(j, Y, pars) {
    with (as.list(Y), {
      if (j == 1) return (phi - f - 1)
      if (j == 2) return (f)
    })
  }
  
  sguess <- seq(0, sigma, le = N)
  yguess <- matrix((rep(c(phi0_val, phi0_val - 1), N)), ncol = N,
                   dimnames = list(c("phi", "f")))
  out <- bvpcol(x = s, func = func, bound = bound, xguess = sguess,
                yguess = yguess, leftbc = 1)[, 1:2]
  out[, 1] <- out[, 1] * D / v_val
  out
}

out <- map2(v, phi0, phi_values)

plot(1, 1, xlim = c(0, 6), ylim = 0:1, type = "l",
     xlab = "position of the gut (cm)", ylab = "food concentration (mM)")

walk2(out, colors, function(x, y) lines(x[, 1], x[, 2], col = y, lwd = 2))