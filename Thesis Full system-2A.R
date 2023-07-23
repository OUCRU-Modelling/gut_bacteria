rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L           <- 6
N           <- 10031                       ### Spatial discretization points, must be odd in order to obtain the exact stationary solution
v           <- c(0.7849518, 0.7849518, 0.6856183, 1.505119)       ### flow velocity of Food, Antibiotic, Bacteria and Mutant respectively
D           <- c(0.2437193, 0.2437193, 0.86, 5.470337)           ### Diffusion coefficient of Food, Antibiotic, Bacteria and Mutant respectively
k           <- 0.1                         ### Monod constant
F_in        <- 1/v[1]                      ### Food concentration at the entrance of the gut
A_in        <- 1/v[2]                      ### Antibiotic concentration at the entrance of the gut

##### Assigning the parameters values #####
r           <- 0.42                        ### Growth rate of Bacteria 
alpha       <- 6.13*10^8                   ### Yield of Food to Bacteria and mutants respectively
beta        <- 6.13*10^7                   ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- 0.09                        ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria 
alpha1      <- (r[1]*F_in)/(k+F_in)
alpha2      <- (A_in)/(A_50[1] + A_in)
delta_wo    <- (alpha1 - (v[3]^2)/(3.2*D[3]))*(1/alpha2)
delta_max1  <- 0.17
tmax        <- 700
times       <- seq(0, tmax,len=200)                        ### discretization of times
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x           <- xgrid$x.mid                 ### We should cho x.mid rather than x.int ### Our discretization points
F_ini       <- rep(1,N)*F_in*0.9
A_ini       <- rep(1,N)*A_in 
B_ini       <- rep(1,N)*0.1*(alpha*F_in)

Sol_system  <- function(t, Y, parms) {
  Food      <- Y[1:N]
  A         <- Y[(N+1):(2*N)]
  B         <- Y[((2*N)+1):(3*N)]
  dFood     <- tran.1D(C = Food, D = D[1], flux.up = 1, flux.down = NULL, v = v[1], dx = xgrid)$dC + -(r/alpha)*B*Food/(k+Food)       
  dA        <- tran.1D(C = A   , D = D[2], flux.up = 1, flux.down = NULL, v = v[2], dx = xgrid)$dC + -(delta_max1/beta)*B*(A)/(A+A_50) 
  dB        <- tran.1D(C = B   , D = D[3], flux.up = 0, flux.down = NULL, v = v[3], dx = xgrid)$dC + r*B*Food/(k+Food) + -delta_max1*B*(A)/(A+A_50)                                          ### tran1D to describe divection diffusion equation
  
  return(list(c(dFood, dA, dB)))
}
yini <- c(F_ini, A_ini, B_ini)
print(system.time(
  out1  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 3, names = c("Food","A","B"), parms = NULL , dimens = N)
))

# Food                  <- out1[(length(times)-10):(length(times)), 2:(N+1)]/F_in
# Antibiotic            <- out1[(length(times)-10):(length(times)), (N+2):(2*N+1)]/A_in
# Bacte                 <- out1[(length(times)-10):(length(times)), (2*(N+1)):(3*N+1)]
# 
# par(mfrow=c(3,1))
# par(mar = c(4, 5, 2, 7) + 0.05)
# for(i in 1:10){
#   if(i==1){
#     plot(x, Food[i, ], xlab='position', ylab='food', col='red', type='l', lwd=2, ylim=c(0, 1))
#   }
#   else{
#     lines(x, Food[i, ], xlab='position', ylab='food', col='red', type='l', lwd=2)
#   }
# }
# 
# for(i in 1:10){
#   if(i==1){
#     plot(x, Antibiotic[i,], col='orange', type='l', lwd=2, ylim=c(0, 1))
#   }
#   else{
#     lines(x, Antibiotic[i,], col='orange', type='l', lwd=2)
#   }
# }
# 
# for(i in 1:10){
#   if(i==1){
#     plot(x, Bacte[i, ], col='blue', type='l', lwd=2)
#   }
#   else{
#     lines(x, Bacte[i, ], col='blue', type='l', lwd=2)
#   }
# }
####### Solving the full system ###########
F_star      <- out1[length(times),2:(N+1)]
B_star      <- out1[length(times),(2*(N+1)):(3*N+1)]
A_star      <- out1[length(times),(N+2):(2*N+1)]
r           <- c(0.42, 0.32)            ### Growth rate of Bacteria and Mutants respectively
alpha       <- c(6.13*10^8, 6.13*10^8)  ### Yield of Food to Bacteria and mutants respectively
beta        <- c(6.13*10^7, 6.13*10^1)      ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- c(0.09, 0.1)             ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria and Mutant respectively     
### Possible wash out limit
delta_max   <- c(delta_max1, 0.00001)  ### Maximum elimination rates of drug killing bacteria and resistant respectively
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x           <- xgrid$x.mid
M0          <- 3.33*10^(-11)  
M_ini       <- rep(1,N)
dx          <- L/N
xm          <- 3
Mini        <- function (x){
  a <- abs(x-xm)
  if (a<=(dx/2)){
    return(M0)
  }                                 ### Initial condition for Mutant
  else{
    return(0)
  }
}
for(i in 1:N){
  M_ini[i] <- Mini(x[i])            ### Asigning the value of function "Mini" for initial condition of Mutant
}

func  <- function(t, Y, parms) {
  Food      <- Y[1:N]
  A         <- Y[(N+1):(2*N)]
  B         <- Y[((2*N)+1):(3*N)]
  M         <- Y[((3*N)+1):(4*N)]
  dFood     <- tran.1D(C = Food, D = D[1], flux.up = 1  , flux.down = NULL, v = v[1], dx = xgrid)$dC + -(r[1]/alpha[1])*B*Food/(k+Food) + -(r[2]/alpha[2])*M*Food/(k+Food)                  
  dA        <- tran.1D(C = A,    D = D[2], flux.up = 1  , flux.down = NULL, v = v[2], dx = xgrid)$dC + -(delta_max[1]/beta[1])*B*(A)/(A+A_50[1]) + -(delta_max[2]/beta[2])*M*(A)/(A+A_50[2]) 
  dB        <- tran.1D(C = B,    D = D[3], flux.up = 0  , flux.down = NULL, v = v[3], dx = xgrid)$dC + r[1]*B*Food/(k+Food) + -delta_max[1]*B*(A)/(A+A_50[1]) 
  dM        <- tran.1D(C = M,    D = D[4], flux.up = 0  , flux.down = NULL, v = v[4], dx = xgrid)$dC + r[2]*M*Food/(k+Food) + -delta_max[2]*M*(A)/(A+A_50[2])
  return(list(c(dFood, dA, dB, dM)))
}
yini <- c(F_star, A_star, B_star, M_ini)
print(system.time(
  out  <- ode.1D(y = yini, func = func, times = times, nspec = 4, names = c("Food","A","B","M"), parms = NULL , dimens = N)
))

Food                    <- out[length(times),2:(N+1)]
Antibiotic              <- out[length(times),(N+2):(2*N+1)]
Bacte                   <- out[length(times),(2*(N+1)):(3*N+1)]
Mutant                  <- out[length(times),(3*N+2):(4*N+1)]
Sum_F <- sum(abs(Food))
Sum_A <- sum(abs(Antibiotic))
Sum_B <- sum(abs(Bacte))
Sum_M <- sum(abs(Mutant))
a <- c(Sum_F, Sum_A, Sum_B, Sum_M)
###### drawing the full profile ########  
par(mar = c(6, 7, 5, 7) + 0.05)
plot(x, Food, xlab='position', ylab='food', col='red', type='l', lwd =2)
axis(side = 2, at = pretty(range(Food)),col="red")
par(new=TRUE)
plot(x, Antibiotic,col="orange",axes=FALSE,xlab="", ylab="", lwd=2, type="l")
axis(side = 2, at = pretty(range(Antibiotic)),col="orange", line = 4)
mtext("antibiotic", side = 2, line = 6,col="orange")
par(new = TRUE)
plot(x, Bacte , col = "blue",
     axes = FALSE, xlab = "",ylab="", lwd=2, type = "l")
axis(side = 4, at = pretty(range(Bacte)),col="blue")
mtext("Bacteria",col="blue", side = 4,line =2)
par(new = TRUE)
plot(x, Mutant , col = "purple",
     axes = FALSE, xlab = "",ylab="", lwd=2, type = "l")
axis(side = 4, at = pretty(range(Mutant)),col="purple", line = 4)
mtext("Mutants",col="purple", side = 4,line = 6)
#image(out, which = "M", legend = TRUE, xlab="times", ylab="position", grid = x)
