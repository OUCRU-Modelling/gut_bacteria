rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L           <- 6
N           <- 10031                 ### Spatial discretization points, must be odd in order to obtain the exact stationary solution
v           <- c(0.04, 0.04, 0.04)     ### flow velocity of Food, Antibiotic, Bacteria 
D           <- c(0.2, 0.2, 0.2)     ### Diffusion coefficient of Food, Antibiotic, Bacteria respectively
k           <- 0.1                   ### Monod constant
F_in        <- 1/v[1]                ### Food concentration at the entrance of the gut
A_in        <- 1/v[2]                ### Antibiotic concentration at the entrance of the gut

##### Assigning the parameters values #####
r           <- 0.42                  ### Growth rate of Bacteria 
alpha       <- 6.13*10^8             ### Yield of Food to Bacteria and mutants respectively
beta        <- 6.13*10^7            ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- 0.09                   ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria and Mutant respectively
tmax        <- 700
alpha1      <- (r*F_in)/(k+F_in)
alpha2      <- (A_in)/(A_50 + A_in)
delta_wo1   <- (alpha1 - (v[3]^2)/(4*D[3]))*(1/alpha2)        ### Possible wash out limit
delta_max   <- 0.203
delta_wo2   <- L/(alpha1- delta_max*alpha2)
times       <- seq(0, tmax,len=200)                           ### discretization of times
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N)    ### generating the gird for our solution
x           <- xgrid$x.mid                                    ### We should cho x.mid rather than x.int 
F_ini       <- rep(1,N)*F_in*0.9                              ### Initial condtitions
A_ini       <- rep(1,N)*A_in
B_ini       <- rep(1,N)*0.1*(alpha*F_in)

Sol_system  <- function(t, Y, parms) {
  Food      <- Y[1:N]
  A         <- Y[(N+1):(2*N)]
  B         <- Y[((2*N)+1):(3*N)]
  dFood     <- tran.1D(C = Food, D = D[1], flux.up = 1, flux.down = NULL, v = v[1], dx = xgrid)$dC + -(r/alpha)*B*Food/(k+Food)       
  dA        <- tran.1D(C = A   , D = D[2], flux.up = 1, flux.down = NULL, v = v[2], dx = xgrid)$dC + -(delta_max/beta)*B*(A)/(A+A_50) 
  dB        <- tran.1D(C = B   , D = D[3], flux.up = 0, flux.down = NULL, v = v[3], dx = xgrid)$dC + r*B*Food/(k+Food) + -delta_max*B*(A)/(A+A_50)                                          ### tran1D to describe divection diffusion equation
  
  return(list(c(dFood, dA, dB)))
}
yini <- c(F_ini, A_ini, B_ini)
print(system.time(
  out  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 3, names = c("Food","A","B"), parms = NULL , dimens = N)
))

out[,(2*(N+1)):(3*N+1)] <- out[,(2*(N+1)):(3*N+1)]
Food                    <- out[length(times),2:(N+1)]
Bacte                   <- out[length(times),(2*(N+1)):(3*N+1)]
Antibiotic              <- out[length(times),(N+2):(2*N+1)]
outtime                 <- seq(from = tmax-60, to = tmax, by = 10)
###### Draw the full non-mutant profile
par(mar = c(6, 5, 5, 7) + 0.05)
plot(x, Food, xlab='position', ylab='food', col='red', type='l')
par(new=TRUE)
plot(x, Antibiotic,col="orange",axes=FALSE,xlab="", ylab="", lwd=1.5, type="l")
axis(side = 4, at = pretty(range(Antibiotic)),col="orange", line = 4)
mtext("antibiotic", side = 4, line = 6,col="orange")
par(new = TRUE)
plot(x, Bacte , col = "blue",              
     axes = FALSE, xlab = "",ylab="", lwd=1.5, type = "l")
axis(side = 4, at = pretty(range(Bacte)),col="blue")      

mtext("Bacteria B ( x10 ^9/mL)",col="blue", side = 4,line =2)
