library(bvpSolve)
library(ReacTran)
library(deSolve)
N <- 500
v <- 0.5
k <- 0.1/(1/v)
D <- 0.4 # Assigning the parameters values
r <- 0.42
lamb <- (r*D)/(v^2) 
alpha <- 6.13*(10^8)
M0 <-  3.33*10^9
xm <- 1
times <- seq(0, 1, len=25)
xgrid <- setup.grid.1D (x.up = 0, x.down = L, N = N)

func <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    
    dF <- -(r/alpha)*(B+M)*Food/(k+Food)+ advection.1D(C=Food, v=v,cx=xgrid)$dC  + tran.1D(C = Food, D = D, flux.up = v*F_in, flux.down = NULL,  dx = xgrid)$dC
    dB <- r*(M*Food)/(k+Food) + advection.1D(C=B, v=v,cx=xgrid)$dC  + tran.1D(C = B, D = D, dx = xgrid)$dC
    dM <- r*(M*Food)/(k+Food) + advection.1D(C=M, v=v,cx=xgrid)$dC  + tran.1D(C = M, D = D, dx = xgrid)$dC
    
    list(c(dF, dB, dM))})
}

parameters <- list(N = N, v=v , k=k, D=D, r=r, lamb=lamb, alpha=alpha, M0=M0)