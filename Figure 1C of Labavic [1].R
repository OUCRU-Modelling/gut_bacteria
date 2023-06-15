rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L               <- 6
k               <- 0.1
r               <- 0.42
alpha           <- 6.13*(10^8)
### it reached the stationary state
v               <- seq(0,2, len=30)
D               <- seq(0,50, len=30)
lamb            <- matrix(rep(1, len=900),ncol=length(v),nrow=length(D))
conc_profile    <- matrix(rep(1, len=900),ncol=length(v),nrow=length(D))
F_in            <- matrix(rep(1, len=900),ncol=length(v),nrow=length(D))
i               <- 1
j               <- 1
while (i <= length(v)){
while (j <= length(D)){
lamb[i, j]      <- (r*D[j])/(v[i]^2)
xic             <- (L*v[i])/D[j]
tmax            <- 500  
N               <- 6000       ### N must be very large in order to abtain the exact solution compared to stationary system (the old code)
F_in[i, j]      <- 1/v[i]
vec_F_in        <- rep(1,N)*F_in[i, j]               
times           <- seq(0, tmax,len=100)                        ### discretization of times
xgrid           <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x               <- xgrid$x.int                                 ### Our discretization points
F_ini           <- vec_F_in*0.9
B_ini           <- 0.1*alpha*vec_F_in
Sol_system      <- function(t, Y, parms) {
          Food  <- Y[1:N]
          B     <- Y[(N+1):(2*N)]
          dFood <- -(r/alpha)*B*Food/(k+Food) + tran.1D(C = Food, D = D, flux.up = 1     , flux.down = NULL, v=v, dx = xgrid)$dC
          dB    <- r*B*Food/(k+Food)          + tran.1D(C = B   , D = D, flux.up = 0     , flux.down = NULL, v=v, dx = xgrid)$dC
          return(list(c(dFood, dB)))
}
yini <- c(F_ini, B_ini)                                                                           

out  <- ode.1D(y = yini, func = Sol_system, times = times,nspec = 2, names = c("Food","B"), parms = NULL , dimens = N)

conc_profile[i, j]  <- (out[tmax,2]-out[tmax,N+1])/F_in[i,j]
j <- j+1
}
i <- i+1
}
