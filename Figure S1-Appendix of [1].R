rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L         <- 6
v         <- seq(0.05, 2, len=5 )
k         <- 0.1
F_in      <- 1/v
D         <- seq(0.05, 30,len=5)        # Assigning the parameters values
r         <- 0.42
alpha     <- 6.13*(10^8)
tmax      <- 500        ### I choose tmax=500 because in the article, they say when our system pass t=500hrs
                        ### it reached the stationary state
N         <- 6000       ### N must be very large in order to abtain the exact solution compared to stationary system (the old code)
i <- 1
j <- 1
while(i<=length(v)){
while(j<=length(D)){
vec_F_in  <- rep(1,N)*F_in[i]               
times     <- seq(0, tmax,len=100)                        ### discretization of times
xgrid     <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x         <- xgrid$x.int                                 ### Our discretization points
F_ini     <- vec_F_in*0.9
B_ini     <- 0.1*alpha*vec_F_in
Sol_system <- function(t, Y, parms) {
  Food  <- Y[1:N]
  B     <- Y[(N+1):(2*N)]
  dFood <- -(r/alpha)*B*Food/(k+Food) + tran.1D(C = Food, D = D[j], flux.up = 1     , flux.down = NULL, v=v[i], dx = xgrid)$dC
  dB    <- r*B*Food/(k+Food)          + tran.1D(C = B   , D = D[j], flux.up = 0     , flux.down = NULL, v=v[i], dx = xgrid)$dC
  return(list(c(dFood, dB)))
}
yini <- c(F_ini, B_ini)                                                                           
print(system.time(
  out  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 2, names = c("Food","B"), parms = NULL , dimens = N)
))
j <- j+1
}
  i <- i+1
}
#out[,2:(N+1)]   <- out[,2:(N+1)]/F_in
#outtime <- seq(from = tmax-60, to = tmax, by = 20)
#matplot.1D(out, which = "Food", ylim = c(0, 1), las = 1, xlim = c(0,6), subset = time %in% outtime, grid = xgrid$x.mid , xlab="x", ylab='Food', main = "Solution 2", type='l', lwd = 2, col= 'red') ### plot in 2D each of Food; Bacteria or Mutant

