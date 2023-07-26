rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
library(parallel)
library(ggplot2)
library(reshape2)
library(tibble)
library(fields)
set.seed(3000)
L           <- 6
N           <- 10031                        ### Spatial discretization points, must be odd in order to obtain the exact stationary solution

### Generating the vector of velocitys and Diffusion coefficients ###
D             <- sort(runif(120, min = 0.01, max = 55), decreasing = FALSE)        
D_B           <- sort(runif(120, min = 0.01, max = 55), decreasing = FALSE)
D_M           <- sort(runif(120, min = 0.01, max = 55), decreasing = FALSE)
v             <- sort(runif(115, min = 0.01, max = 3), decreasing = FALSE)        
v_B           <- sort(runif(115, min = 0.01, max = 3), decreasing = FALSE)
v_M           <- sort(runif(115, min = 0.01, max = 3), decreasing = FALSE)
delta_max_B   <- 0.2
delta_max_M   <- 0.01
numCores      <- detectCores()                            ### Numbers of cores to be using parallel

##### Assigning the parameters values #####
k           <- 0.1                          ### Monod constant

r           <- c(0.42, 0.32)                ### Growth rate of Bacteria and Mutants respectively
alpha       <- c(6.13*10^8, 6.13*10^8)      ### Yield of Food to Bacteria and mutants respectively
beta        <- c(6.13*10^7, 6.13*10)        ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- c(0.09, 0.1)                 ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria and Mutant respectively     
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x           <- xgrid$x.mid
M0          <- 3.33*10^(-11)  
M_ini       <- rep(1,N)
dx          <- L/N
xm          <- 2
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
tmax        <- 700
ite         <- vector('list', length(v)*length(D)) ### Creating the list index to be using in function in parallel
l <- 1
for(i in 1 : length(D)){
  for(j in 1 : length(v)){
    ite[[l]] <- c(i, j)
    l        <- l + 1
  }
}
index <- seq(1, length(v)*length(D), by = 1)       ### Index in parallel computation

#### Function to use in parallel #######
Heat_map <- function(n){                           ### Main function to be compiled in parallel computation
  a           <- array(unlist(ite[n]), dim = c(1, 2))
  i           <- a[1]
  j           <- a[2]
  F_in        <- 1/v[j]                
  A_in        <- 1/v[j] 
  times       <- seq(0, tmax,len = 200)                         ### discretization of times
  
  F_ini1      <- rep(1,N)*F_in*0.9                              ### Initial condtitions
  A_ini1      <- rep(1,N)*A_in
  B_ini1      <- rep(1,N)*0.1*(alpha[1]*F_in)
  
  Sol_system  <- function(t, Y, parms) {
    Food      <- Y[1:N]
    A         <- Y[(N+1):(2*N)]
    B         <- Y[((2*N)+1):(3*N)]
    dFood     <- tran.1D(C = Food, D = D[i],   flux.up = 1, flux.down = NULL, v = v[j], dx = xgrid)$dC + -(r[1]/alpha[1])*B*Food/(k+Food)       
    dA        <- tran.1D(C = A   , D = D[i],   flux.up = 1, flux.down = NULL, v = v[j], dx = xgrid)$dC + -(delta_max_B/beta[1])*B*(A)/(A+A_50[1]) 
    dB        <- tran.1D(C = B   , D = D_B[i], flux.up = 0, flux.down = NULL, v = v_B[j], dx = xgrid)$dC + r[1]*B*Food/(k+Food) + -delta_max_B*B*(A)/(A+A_50[1])                                          ### tran1D to describe divection diffusion equation
    
    return(list(c(dFood, dA, dB)))
  }
  yini1 <- c(F_ini1, A_ini1, B_ini1)
  print(system.time(
    out1  <- ode.1D(y = yini1, func = Sol_system, times = times, nspec = 3, names = c("Food","A","B"), parms = NULL , dimens = N)
  ))
  ##### Solving Full System #######
  F_star      <- out1[length(times), 2:(N+1)]
  B_star      <- out1[length(times), (2*(N+1)):(3*N+1)]
  A_star      <- out1[length(times), (N+2):(2*N+1)]
  
  
  
  func  <- function(t, Y, parms) {
    Food      <- Y[1:N]
    A         <- Y[(N+1):(2*N)]
    B         <- Y[((2*N)+1):(3*N)]
    M         <- Y[((3*N)+1):(4*N)]
    dFood     <- tran.1D(C = Food, D = D[i], flux.up = 1  ,   flux.down = NULL, v = v[j], dx = xgrid)$dC + -(r[1]/alpha[1])*B*Food/(k+Food) + -(r[2]/alpha[2])*M*Food/(k+Food)                  
    dA        <- tran.1D(C = A,    D = D[i], flux.up = 1  ,   flux.down = NULL, v = v[j], dx = xgrid)$dC + -(delta_max_B/beta[1])*B*(A)/(A+A_50[1]) + -(delta_max_M/beta[2])*M*(A)/(A+A_50[2]) 
    dB        <- tran.1D(C = B,    D = D_B[i], flux.up = 0  , flux.down = NULL, v = v_B[j], dx = xgrid)$dC + r[1]*B*Food/(k+Food) + -delta_max_B*B*(A)/(A+A_50[1]) 
    dM        <- tran.1D(C = M,    D = D_M[i], flux.up = 0  , flux.down = NULL, v = v_M[j], dx = xgrid)$dC + r[2]*M*Food/(k+Food) + -delta_max_M*M*(A)/(A+A_50[2])
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
  
  Sum_F   <- sum(abs(Food))
  Sum_A   <- sum(abs(Antibiotic))
  Sum_B   <- sum(abs(Bacte))
  Sum_M   <- sum(abs(Mutant))
  
  return(c(Sum_F, Sum_A, Sum_B, Sum_M))
}
system.time({
  results <- mclapply(index, Heat_map, mc.cores = numCores)
})

conc_profile            <- array(unlist(results), dim = c(4, length(D)*length(v)))
Sum_Food                <- matrix(conc_profile[1, ], ncol = length(delta_max_B), nrow = length(v))
Sum_Antibiotic          <- matrix(conc_profile[2, ], ncol = length(delta_max_B), nrow = length(v))
Sum_Bacteria            <- matrix(conc_profile[3, ], ncol = length(delta_max_B), nrow = length(v))
Sum_Mutant              <- matrix(conc_profile[4, ], ncol = length(delta_max_B), nrow = length(v))