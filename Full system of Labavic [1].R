rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L         <- 6
v         <- 0.5
k         <- 0.1
F_in      <- 1/v
D         <- 0.4       # Assigning the parameters values
r         <- 0.42
alpha     <- 6.13*(10^8)
tmax      <- 700        ### I choose tmax=500 because in the article, they say when our system pass t=500hrs
### it reached the stationary state
N         <- 10000       ### N must be very large in order to abtain the exact solution compared to stationary system (the old code)
vec_F_in  <- rep(1,N)*F_in               
times     <- seq(0, tmax, len=200)                        ### discretization of times
xgrid     <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x         <- xgrid$x.mid                                 ### Our discretization points
F_ini     <- vec_F_in*0.9
B_ini     <- 0.1*alpha*vec_F_in
Sol_system <- function(t, Y, parms) {
  Food  <- Y[1:N]
  B     <- Y[(N+1):(2*N)]
  dFood <- -(r/alpha)*B*Food/(k+Food) + tran.1D(C = Food, D = D, flux.up = 1     , flux.down = NULL, v=v, dx = xgrid)$dC
  dB    <- r*B*Food/(k+Food)          + tran.1D(C = B   , D = D, flux.up = 0     , flux.down = NULL, v=v, dx = xgrid)$dC
  return(list(c(dFood, dB)))
}
yini <- c(F_ini, B_ini)                                                                           
print(system.time(
  out1  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 2, names = c("Food","B"), parms = NULL , dimens = N)
))

#### Solve the system with mutant first appear at xm

vec_F_in  <- rep(1,N)*F_in                ### because 
Food_star <- out1[length(times), 2:(N+1)]*F_in                 ### The initial condition of PDE system, is in Food whereas in without mutant, the solution is F/F_in
B_star    <- alpha*(vec_F_in - Food_star) ### The equation of B_star and F_star in no mutant system (my previous reports)
xm        <- 3                            ### the first location mutant appear
M0        <- 3.33*10^(-9)                ### mutant concentration at this location (xm), initial concentration of Mutant
times     <- seq(0, 1,len=100)            ### discretization of times
xgrid     <- setup.grid.1D (x.up = 0, x.down = L, N = N) ## generating the gird for our solution
x         <- xgrid$x.mid                  ### Our discretization points
F_ini     <- Food_star
B_ini     <- B_star                       ### Our initial condition of the system (at t = 0)
M_ini     <- rep(1,N)
dx        <- L/N
Mini      <- function (x){
  a <- abs(x-xm)
  if (a<=(dx/2)){
    return(M0)
  }                                 ### Initial condition for Mutant
  else{
    return(0)
  }
}
for(i in 1:N){
  M_ini[i] <- Mini(x[i])       ### Asigning the value of function "Mini" for initial condition of Mutant
}
yini <- c(F_ini, B_ini, M_ini) ### Initial condition to use in ode1D

### Solving the full system using ode.1D
func <- function(t, Y, parms) {
  Food  <- Y[1:N]
  B     <- Y[(N+1):(2*N)]
  M     <- Y[((2*N)+1):(3*N)]
  dFood <- -(r/alpha)*(B+M)*Food/(k+Food) + tran.1D(C = Food, D = D, flux.up = v*F_in, flux.down = NULL, v=v, dx = xgrid)$dC
  dB    <- r*(B*Food)/(k+Food)            + tran.1D(C = B   , D = D, flux.up = 0     , flux.down = NULL, v=v, dx = xgrid)$dC   ### tran1D to describe divection diffusion equation
  dM    <- r*(M*Food)/(k+Food)            + tran.1D(C = M   , D = D, flux.up = 0     , flux.down = NULL, v=v, dx = xgrid)$dC
  return(list(c(dFood, dB, dM)))
}

print(system.time(
  out  <- ode.1D(y = yini, func = func, times = times, nspec = 3, names = c("Food","B","M"), parms = NULL , dimens = N)
))
M           <- 10
Food        <- out[, 2:(N+1)]
Bacte       <- out[,(N+2):(2*N+1)]
Mutant      <- out[,(2*(N+1)):(3*N+1)]
par(mfrow=c(3,1))
par(mar = c(4, 5, 2, 7) + 0.05 )
for (i in (length(times)-M):length(times)){
  if(i==length(times)-M){
    plot(x, Food[i,], xlab = 'position x', ylab = 'Food', type='l', col='red', lwd = 1.5)
    plot(x, Bacte[i,], xlab = 'position x', ylab = 'bacteria', type='l', col='blue', lwd = 1.5)
    plot(x, Mutant[i,], type='l', xlab = 'position x', ylab = 'Neutral Mutant', col='purple', lwd = 1.5)
  }
  else{
    lines(x, Food[i,], xlab = 'position x', ylab = 'Food', type='l', col='red', lwd = 1.5)
    lines(x, Bacte[i,], xlab = 'position x', ylab = 'bacteria', type='l', col='blue', lwd = 1.5)
    lines(x, Mutant[i,], type='l', xlab = 'position x', ylab = 'Neutral Mutant', col='purple', lwd = 1.5)
  }
}