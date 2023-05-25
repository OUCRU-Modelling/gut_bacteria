rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L     <- 6
N     <- 10000
v     <- 0.5
k     <- 0.1
F_in  <- 1/v
kappa <- 0.1/(1/v)  # Note that we only use kappa to solve the initial condition for the with mutant system 
D     <- 0.2        # Assigning the parameters values
r     <- 0.42
lamb  <- (r*D)/(v^2)
xic   <- (L*v)/D
alpha <- 6.13*(10^8)
tmax  <- 510
vec_F_in  <- rep(1,N)*F_in               
times     <- seq(0, tmax,len=100)                        ### discretization of times
xgrid     <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x         <- xgrid$x.int                                 ### Our discretization points
F_ini     <- vec_F_in*0.9
B_ini     <- 0.1*alpha*vec_F_in

#### Solving the non-Mutant system using the new approach
func <- function(t, Y, parms) {
  Food  <- Y[1:N]
  B     <- Y[(N+1):(2*N)]
  dFood <- -(r/alpha)*B*Food/(k+Food) + tran.1D(C = Food, D = D, flux.up = 1     , flux.down = NULL, v=v, dx = xgrid)$dC
  dB    <- r*B*Food/(k+Food)          + tran.1D(C = B   , D = D, flux.up = 0     , flux.down = NULL, v=v, dx = xgrid)$dC
  return(list(c(dFood, dB)))
}
yini <- c(F_ini, B_ini)
print(system.time(
  out  <- ode.1D(y = yini, func = func, times = times,nspec = 2, names = c("Food","B"), parms = NULL , dimens = N)
))

#### Solve the system with mutant first appear at xm
tmax      <- 510
Food_star <- out[100,2:(N+1)]                ### The initial condition of PDE system, obtain by extract the solution "out"
B_star    <- out[100,(N+2):(2*N+1)]          ### The equation of B_star and F_star in no mutant system 
xm        <- seq(0,6,by=0.1)                 ###  create the array of xm to drwa the Figure 2
times     <- seq(0, tmax,len=200)            ### discretization of times
xgrid     <- setup.grid.1D (x.up = 0, x.down = L, N = N) ## generating the gird for our solution
x         <- xgrid$x.mid                     ### Our Spatial discretization points
F_ini     <- Food_star
B_ini     <- B_star                          ### Our initial condition of the system (at t = 0)
M_ini     <- rep(1,N)
dx        <- L/(N-1)                         ### Spatial step sizes
M0        <- (3.33*10^(-11))/dx              ### This formular is in the maintext page 3
Rxm       <- rep(1, length(xm))              ### creating the array to contains R at xm, where R is reproduction rate in unit volume and time
                                             ### in maintext, page 4
if(reinitialise_mutant <- TRUE){             ### When you turn it into TRUE, Rstudio will run the code in the bracket
  Mini      <- function (x,x0){
    a <- abs(x-x0)
    if (a<=(dx/2)){
      return(M0)
    }                                        ### Initial condition for Mutant
    else{
      return(0)
    }
  }
}
i  <- 1

#### Initializing the loop
while (i <= length(xm)){
  
  for(j in 1:N){
    M_ini[j] <- Mini(x[j],xm[i]) 
  }                                         ### Asigning the value of function "Mini" for initial condition of Mutant
 
yini <- c(F_ini, B_ini, M_ini)              ### Initial condition to use in ode1D


### Solving the full system using ode.1D
Sol_system <- function(t, Z, parms) {
  Food  <- Z[1:N]
  B     <- Z[(N+1):(2*N)]
  M     <- Z[((2*N)+1):(3*N)]
  dFood <- -(r/alpha)*(M+B)*Food/(k+Food) + tran.1D(C = Food, D = D, flux.up = v*F_in , flux.down = NULL, v=v, dx = xgrid)$dC
  dB    <- r*B*Food/(k+Food)              + tran.1D(C = B   , D = D, flux.up = 0      , flux.down = NULL, v=v, dx = xgrid)$dC   ### tran1D to describe divection diffusion equation
  dM    <- r*M*Food/(k+Food)              + tran.1D(C = M   , D = D, flux.up = 0      , flux.down = NULL, v=v, dx = xgrid)$dC
  return(list(c(dFood, dB, dM)))
}

print(system.time(
  out1  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 3, names = c("Food","B","M"), parms = NULL , dimens = N)
))

FOOD          <- out1[,2:(N+1)]             ### Extracting the food in the full system
Bacte         <- out1[,(N+2):(2*N+1)]       ### Extracting the Bacteria in the full system
Mutant        <- out1[,(2*(N+1)):(3*N+1)]   ### Extracting the Mutant in the full system
ro            <- rep(1,N)                   ### ro, the formular in maintext page 3 bottem line
ro            <- r*(FOOD[length(times),])/(kappa + FOOD[length(times),])
R             <- rep(1,N)
R             <- (Bacte[length(times),])*ro ### R, this formular in the maintext page 4
index_xm      <- (((xm[i])/dx)) + 1         ### index of xm in spatial discretization
Rxm[i]        <- R[index_xm]*10^-8          ### R, at the xm
i <- i +  1
}
plot(xm, Rxm, ylab = 'reproduction rate in 
     unit volume and time R', xlab = 'mutant introduction position xm (cm)',
     type = 'l', col = 'pink',lwd =2)       ### Drar pink line in Fig 2
