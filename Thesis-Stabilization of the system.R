rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
L           <- 6
N           <- 10031                       ### Spatial discretization points, must be odd in order to obtain the exact stationary solution
v           <- c(0.45, 0.45, 0.45, 0.35)         ### flow velocity of Food, Antibiotic, Bacteria and Mutant respectively
D           <- c(0.2, 0.2, 1.5, 1.5)           ### Diffusion coefficient of Food, Antibiotic, Bacteria and Mutant respectively
k           <- 0.1                         ### Monod constant
F_in        <- 1/v[1]                      ### Food concentration at the entrance of the gut
A_in        <- 1/v[2]                      ### Antibiotic concentration at the entrance of the gut
del         <- -0.158
##### Assigning the parameters values #####
r           <- 0.42                        ### Growth rate of Bacteria 
alpha       <- 6.13*10^8                   ### Yield of Food to Bacteria and mutants respectively
beta        <- 10^8                        ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- 0.1                         ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria 
tmax1       <- 900
alpha1      <- (r*F_in)/(k+F_in)
alpha2      <- (A_in)/(A_50 + A_in)
delta_wo    <- (alpha1 - (v[3]^2)/(4*D[3]))*(1/alpha2)      ### Wash out limit
delta_max   <- delta_wo + del 
times       <- seq(0, tmax1,len=200)                        ### discretization of times
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x           <- xgrid$x.mid                                 ### Our discretization points
F_ini       <- rep(1,N)*F_in*0.9
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
  out1  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 3, names = c("Food","A","B"), parms = NULL , dimens = N)
))
M           <- 18
Food        <- out1[, 2:(N+1)]
Anti        <- out1[,(N+2):(2*N+1)]
B           <- out1[,(2*(N+1)):(3*N+1)]
par(mfrow=c(3,1))
par(mar = c(4, 5, 2, 7) + 0.05 )
for (i in (length(times)-M):length(times)){
  if(i==length(times)-M){
    plot(x, Food[i,], xlab = 'position x', ylab = 'Food', type='l', col='red', lwd = 1.5)
    plot(x, Anti[i,], xlab = 'position x', ylab = 'drug', type='l', col='orange', lwd = 1.5)
    plot(x, B[i,], type='l', xlab = 'position x', ylab = 'bacteria', col='blue', lwd = 1.5)
  }
  else{
    lines(x, Food[i,], xlab = 'position x', ylab = 'Food', type='l', col='red', lwd = 1.5)
    lines(x, Anti[i,], xlab = 'position x', ylab = 'drug', type='l', col='orange', lwd = 1.5)
    lines(x, B[i,], type='l', xlab = 'position x', ylab = 'bacteria', col='blue', lwd = 1.5)
  }
}


######### Solving the full system ###########
F_star      <- out1[length(times),2:(N+1)]
B_star      <- out1[length(times),(2*(N+1)):(3*N+1)]
A_star      <- out1[length(times),(N+2):(2*N+1)]
r           <- c(0.42, 0.335)          ### Growth rate of Bacteria and Mutants respectively
alpha       <- c(6.13*10^8, 6.13*10^8)     ### Yield of Food to Bacteria and mutants respectively
beta        <- c(10^8, 10^2)        ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- c(0.15, 0.15)            ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria and Mutant respectively
alpha1      <- (r[1]*F_in)/(k+F_in)
alpha2      <- (A_in)/(A_50[1] + A_in)
delta_wo    <- (alpha1 - (v[3]^2)/(4*D[3]))*(1/alpha2)        ### Possible wash out limit
delta_max   <- c(delta_wo + del, 0.000000001)  ### Maximum elimination rates of drug killing bacteria and resistant respectively
tmax2       <- 9900
times2      <- seq(0, tmax2, len = 300)
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x           <- xgrid$x.mid
M0          <- 3.33*10^(-11)
M_ini       <- rep(1,N)
dx          <- L/N
xm          <- 1
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
  dB        <- tran.1D(C = B,    D = D[3], flux.up = 0  , flux.down = NULL, v = v[3], dx = xgrid)$dC + r[1]*B*Food/(k+Food) - delta_max[1]*B*(A)/(A+A_50[1])
  dM        <- tran.1D(C = M,    D = D[4], flux.up = 0  , flux.down = NULL, v = v[4], dx = xgrid)$dC + r[2]*M*Food/(k+Food) - delta_max[2]*M*(A)/(A+A_50[2])
  return(list(c(dFood, dA, dB, dM)))
}
yini <- c(F_ini, A_ini, B_ini, M_ini)
print(system.time(
  out  <- ode.1D(y = yini, func = func, times = times2, nspec = 4, names = c("Food","A","B","M"), parms = NULL , dimens = N)
))
M                 <- 18
Food              <- out[,2:(N+1)]
Antibiotic        <- out[,(N+2):(2*N+1)]
Bacte             <- out[,(2*(N+1)):(3*N+1)]
Mutant            <- out[,(3*N+2):(4*N+1)]
par(mfrow=c(2,2))
par(mar = c(4, 5, 2, 2.5) + 0.05 )
i                 <- length(times2)-M
plot(x, Food[i,], xlab = 'position x', ylab = 'Food', type='l', col='red', lwd = 1.5)
for (j in (length(times2)-M+1):length(times2)){
  lines(x, Food[j,], xlab = 'position x', ylab = 'Food', type='l', col='red', lwd = 1.5)}

plot(x, Antibiotic[i,], xlab = 'position x', ylab = 'drug', type='l', col='orange', lwd = 1.5)
for (j in (length(times2)-M+1):length(times2)){
  lines(x, Antibiotic[j,], xlab = 'position x', ylab = 'drug', type='l', col='orange', lwd = 1.5)}

plot(x, Bacte[i,], type='l', xlab = 'position x', ylab = 'Wild-type bacteria', col='blue', lwd = 1.5)
for (j in (length(times2)-M+1):length(times2)){
  lines(x, Bacte[j,], xlab = 'position x', ylab = 'Wild-type bacteria',type='l', col='blue', lwd = 1.5)}

plot(x, Mutant[i,], xlab = 'position x', ylab = 'Resistant Bacteria', type='l', col='purple', lwd = 1.5)#, ylim=c(,5.9*10^-74))
for (j in (length(times2)-M+1):length(times2)){
  lines(x, Mutant[j,], xlab = 'position x', ylab = 'Resistant Bacteria', type='l', col='purple', lwd = 1.5)}

print(delta_max)
print(delta_wo)
