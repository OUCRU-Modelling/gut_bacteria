rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L           <- 6
N           <- 10031                  ### Spatial discretization points, must be odd in order to obtain the exact stationary solution
v           <- c(0.5,0.5,0.5)             ### flow velocity of Food, Antibiotic, Bacteria and Mutant respectively
D           <- c(0.2,0.2,0.2)             ### Diffusion coefficient of Food, Antibiotic, Bacteria and Mutant respectively
k           <- 0.1                    ### Monod constant
F_in        <- 1/v[1]                 ### Food concentration at the entrance of the gut
A_in        <- 1/v[2]                 ### Antibiotic concentration at the entrance of the gut

##### Assigning the parameters values #####
r           <- 0.42                       ### Growth rate of Bacteria 
alpha       <- 6.13*10^8              ### Yield of Food to Bacteria and mutants respectively
beta        <- 0.08208                ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- 4.91                   ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria and Mutant respectively
delta_max   <- 1.13083
tmax        <- 510
times       <- seq(0, tmax,len=100)                        ### discretization of times
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = 200) ### generating the gird for our solution
x           <- xgrid$x.mid    ### We should cho x.mid rather than x.int ### Our discretization points
F_ini       <- rep(1,N)*F_in*0.9
A_ini       <- rep(1,N)*A_in*0.9
B_ini       <- rep(1,N)*0.1*(alpha*F_in-beta*A_in)

Sol_system  <- function(t, Y, parms) {
  Food      <- Y[1:N]
  A         <- Y[(N+1):(2*N)]
  B         <- Y[((2*N)+1):(3*N)]
  dFood     <- -(r/alpha[1])*B*Food/(k+Food)  + tran.1D(C = Food, D = D[1] , flux.up = 1   , flux.down = NULL, v=v[1], dx = xgrid)$dC
  dA        <- -(delta_max/beta[1])*B*(A^k)/(A^k+A_50[1]) + -(delta_max/beta[2])*M*(A^k)/(A^k+A_50[2])   + tran.1D(C = Food, D = D[2] , flux.up = 1   , flux.down = NULL, v=v[2], dx = xgrid)$dC
  dB        <- r*B*Food/(k+Food) + -delta_max*B*(A^k)/(A^k+A_50[1]) + tran.1D(C = B , D = D[3], flux.up = 0 , flux.down = NULL, v=v[3] , dx = xgrid)$dC                                         ### tran1D to describe divection diffusion equation
  
  return(list(c(dFood, dB, dM)))
}
yini <- c(F_ini, A_ini, B_ini)
print(system.time(
  out  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 3, names = c("Food","A","B"), parms = NULL , dimens = N)
))