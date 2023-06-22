rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L           <- 6
N           <- 10031                  ### Spatial discretization points, must be odd in order to obtain the exact stationary solution
v           <- c(1,2,3,4)             ### flow velocity of Food, Antibiotic, Bacteria and Mutant respectively
D           <- c(1,2,3,4)             ### Diffusion coefficient of Food, Antibiotic, Bacteria and Mutant respectively
k           <- 0.1                    ### Monod constant
F_in        <- 1/v[1]                 ### Food concentration at the entrance of the gut
A_in        <- 1/v[2]                 ### Antibiotic concentration at the entrance of the gut

##### Assigning the parameters values #####
r           <- c(1, 2)                ### Growth rate of Bacteria and Mutants respectively
alpha       <- c(1, 2)                ### Yield of Food to Bacteria and mutants respectively
beta        <- c(1, 2)                ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- c(1, 2)                ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria and Mutant respectively
delta_max   <- 1
tmax        <- 510
times       <- seq(0, tmax,len=100)                        ### discretization of times
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = 200) ### generating the gird for our solution
x           <- xgrid$x.mid    ### We should cho x.mid rather than x.int ### Our discretization points


Sol_system <- function(t, Y, parms) {
  Food     <- Y[1:N]
  A        <- Y[(N+1):(2*N)]
  B        <- Y[((2*N)+1):(3*N)]
  M        <- Y[((3*N)+1):(4*N)]
  dFood    <- -(r[1]/alpha[1])*B*Food/(k+Food) + -(r[2]/alpha[2])*M*Food/(k+Food)   + tran.1D(C = Food, D = D[1] , flux.up = 1   , flux.down = NULL, v=v[1], dx = xgrid)$dC
  dA       <- -(delta_max/beta[1])*B*(A^k)/(A^k+A_50[1]) + -(delta_max/beta[2])*M*(A^k)/(A^k+A_50[2])   + tran.1D(C = Food, D = D[2] , flux.up = 1   , flux.down = NULL, v=v[2], dx = xgrid)$dC
  dB       <- r[1]*B*Food/(k+Food) + -delta_max*B*(A^k)/(A^k+A_50[1]) + tran.1D(C = B , D = D[3], flux.up = 0 , flux.down = NULL, v=v[3] , dx = xgrid)$dC                                         ### tran1D to describe divection diffusion equation
  dM       <- r[2]*M*Food/(k+Food) + -delta_max*M*(A^k)/(A^k+A_50[2]) + tran.1D(C = M , D = D[4], flux.up = 0 , flux.down = NULL, v=v[4] , dx = xgrid)$dC
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