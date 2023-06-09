rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L     <- 6
N     <- 10030      ### I also discovered that in order to obtain the exact solution, you also need to set N is odd !!!
v     <- c(1,2,3,4)
D     <- c(1,2,3,4)
k     <- 0.1
F_in  <- 1/v
kappa <- 0.1/(1/v)  # Note that we only use kappa to solve the initial condition for the with mutant system 
        # Assigning the parameters values
r     <- 0.42
lamb  <- (r*D)/(v^2)
xic   <- (L*v)/D
alpha <- 6.13*(10^8)
tmax  <- 510
vec_F_in  <- rep(1,N)*F_in               
times     <- seq(0, tmax,len=100)                        ### discretization of times
xgrid     <- setup.grid.1D (x.up = 0, x.down = L, N = 200) ### generating the gird for our solution
x         <- xgrid$x.mid    ### We should cho x.mid rather than x.int ### Our discretization points
F_ini     <- vec_F_in*0.9
B_ini     <- 0.1*alpha*vec_F_in

Sol_system <- function(t, Y, parms) {
    Food  <- Y[1:N]
    A     <- Y[(N+1):(2*N)]
    B     <- Y[((2*N)+1):(3*N)]
    M     <- Y[((3*N)+1):(4*N)]
    dFood <- -(r_B/alpha_B)*B*Food/(k+Food) + -(r_M/alpha_M)*M*Food/(k+Food) + tran.1D(C = Food, D = D_F, flux.up = 1     , flux.down = NULL, v=v_F, dx = xgrid)$dC
    dB    <- r*B*Food/(k+Food)              + tran.1D(C = B   , D = D, flux.up = 0     , flux.down = NULL, v=v, dx = xgrid)$dC   ### tran1D to describe divection diffusion equation
    dM    <- r*M*Food/(k+Food)              + tran.1D(C = M   , D = D, flux.up = 0     , flux.down = NULL, v=v, dx = xgrid)$dC
    return(list(c(dFood, dB, dM)))
  }
  
  print(system.time(
    out  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 3, names = c("Food","B","M"), parms = NULL , dimens = N)
  ))
  
  

par(mar = c(4, 5, 4, 5) + 0.05 ) 
plot(xm, M_on_B, ylab = 'Ratio of M/B (x10^-19)' 
     , xlab = 'mutant introduction position xm (cm)',
     type = 'l', col = 'black',lwd =2)
par(new = TRUE) 
plot(xm, Rxm, ylab = "", xlab = "", axes = FALSE,
     type = 'l',ylim = c(0, 3.8),las = 1, col = 'pink',lwd =2) 
axis(side = 4, at = pretty(range(Rxm)),col='pink')
mtext("Reproduction rate in 
     unit volume and time R", side = 4, line = 2.9,col="pink")
