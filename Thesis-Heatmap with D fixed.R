rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
library(parallel)
library(ggplot2)
library(reshape2)
library(tibble)
library(fields)

L           <- 6
N           <- 10031                 ### Spatial discretization points, must be odd in order to obtain the exact stationary solution
### Generating the vector of velocitys and Diffusion coefficients ###
v           <- c(seq(0.05, 0.1, len=30), seq(0.105, 0.2, len=30), 
                 seq(0.201, 0.5, len = 25), seq(0.505, 1, len = 25),
                 seq(1.1, 2.6, len = 13))     
v_B         <- c(seq(0.05, 0.1, len=10), seq(0.105, 0.2, len=10), seq(0.201, 0.5, len = 5),
                 seq(0.505, 1, len = 5), seq(1.1, 2.6, len = 3)) 
D           <- 0.4
delta_max   <- seq(0.01, 0.4, len = 100)
numCores    <- detectCores() - 1     ### Numbers of cores to be using parallel
##### Assigning the parameters values #####
k           <- 0.1                   ### Monod constant
r           <- 0.42                  ### Growth rate of Bacteria 
alpha       <- 6.13*10^8             ### Yield of Food to Bacteria 
beta        <- 6.13*10^7
A_50        <- 0.1
tmax        <- 700
ite         <- vector('list', length(delta_max)*length(v)) ### Creating the list index to be using in function in parallel
l <- 1
for(i in 1 : length(delta_max)){
  for(j in 1 : length(v)){
    ite[[l]] <- c(i, j)
    l        <- l + 1
  }
}
index <- seq(1, length(delta_max)*length(v), by = 1)       ### Index in parallel computation

Heat_map <- function(n){                           ### Main function to be compiled in parallel computation
  a           <- array(unlist(ite[n]), dim = c(1,  2))
  i           <- a[1]
  j           <- a[2]
  F_in        <- 1/v[j]                
  A_in        <- 1/v[j] 
  times       <- seq(0, tmax,len = 160)                         ### discretization of times
  xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N)    ### generating the gird for our solution
  x           <- xgrid$x.mid                                    ### We should cho x.mid rather than x.int 
  F_ini       <- rep(1,N)*F_in*0.9                              ### Initial condtitions
  A_ini       <- rep(1,N)*A_in
  B_ini       <- rep(1,N)*0.1*(alpha*F_in)
  
  Sol_system  <- function(t, Y, parms) {
    Food      <- Y[1:N]
    A         <- Y[(N+1):(2*N)]
    B         <- Y[((2*N)+1):(3*N)]
    dFood     <- tran.1D(C = Food, D = D,   flux.up = 1, flux.down = NULL, v = v[j],   dx = xgrid)$dC + -(r/alpha)*B*Food/(k+Food)       
    dA        <- tran.1D(C = A   , D = D,   flux.up = 1, flux.down = NULL, v = v[j],   dx = xgrid)$dC + -(delta_max[i]/beta)*B*(A)/(A+A_50) 
    dB        <- tran.1D(C = B   , D = D_B, flux.up = 0, flux.down = NULL, v = v_B[j], dx = xgrid)$dC + r*B*Food/(k+Food) + -delta_max[i]*B*(A)/(A+A_50)                                          ### tran1D to describe divection diffusion equation
    
    return(list(c(dFood, dA, dB)))
  }
  yini <- c(F_ini, A_ini, B_ini)
  print(system.time(
    out  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 3, names = c("Food","A","B"), parms = NULL , dimens = N)
  ))
  Bacte             <- out[length(times),(2*(N+1)):(3*N+1)]
  Food_ave          <- (out[length(times), 2] - out[length(times), N+1])/F_in
  Drug_ave          <- (out[length(times),(N+2)] - out[length(times), (2*N+1)])/A_in
  return(c(Food_ave, Drug_ave, mean(abs(Bacte))))
  
}
system.time({
  results <- mclapply(index, Heat_map, mc.cores = numCores)
})

F_in        <- 1/v                
A_in        <- 1/v
alpha1      <- (r*F_in)/(k+F_in)
alpha2      <- (A_in)/(A_50 + A_in)

washout_line_delta       <- (alpha1 - (v_B^2)/(3.5*D))/(alpha2)
washout_line_v           <- L*(alpha1 - delta_max*alpha2)
conc_profile             <- array(unlist(results), dim = c (3, length(delta_max)*length(v)))
conc_profile_Food        <- matrix(conc_profile[1, ], ncol = length(delta_max), nrow = length(v))
conc_profile_Antibiotic  <- matrix(conc_profile[2, ], ncol = length(delta_max), nrow = length(v))
conc_profile_Bacteria    <- matrix(conc_profile[3, ], ncol = length(delta_max), nrow = length(v))
x  <- v
y  <- delta_max
n1 <- length(washout_line_delta[washout_line_delta > 0])
n2 <- length(washout_line_v[washout_line_v > 0])

par(mar = c(6, 5, 6, 6.5) + 0.05 )  

#### Heatmap Food ####
image.plot(conc_profile_Food, x=log(x, 20), y=log(y, 20),# ylim = log(c(0.01, 0.4), 20),
           xlab = 'velocity v', ylab = 'delta max',
           main = 'Spatial dependence regarding Food
           average, D = 0.4', col = NULL, axes = TRUE)
contour(z = conc_profile_Food, x=log(x, 20), y=log(y,20), add = TRUE) 
lines(log(x[1:n2], 20), log(washout_line_delta[washout_line_delta > 0], 20), type = 'l'
      , col = 'magenta', lwd = 3)
lines(log(x[1:n1] ,20), log(washout_line_v[washout_line_v > 0], 20),
      col = 'cyan', lwd = 3)
#### Heatmap drug ####
image.plot(conc_profile_Antibiotic, x=log(x, 20), y=log(y, 20), xlab = 'velocity v'
           , ylab = 'delta max',
           main = 'Spatial dependence regarding Antibiotic
           average, D = 0.4', col = NULL, axes=FALSE)
contour(z = conc_profile_Antibiotic,x=log10(x), y=log10(y),  add=TRUE)
lines(log(x[1:n2], 20), log(washout_line_delta[washout_line_delta > 0], 20), type = 'l'
      , col = 'magenta', lwd = 3)
lines(log(x[1:n1] ,20), log(washout_line_v[washout_line_v > 0], 20),
      col = 'cyan', lwd = 3)

### Heatmap Bacteria ####
image.plot(conc_profile_Bacteria, x=log(x, 20), y=log(y, 20), xlab = 'velocity v'
           , ylab = 'delta max',
           main = 'Spatial dependence regarding Bacteria
           average, D = 0.4', col = NULL, axes=FALSE)
contour(z = conc_profile_Food, x=log(x, 20), y=log(y, 20),  add=TRUE)
lines(log(x[1:n2], 20), log(washout_line_delta[washout_line_delta > 0], 20), type = 'l'
      , col = 'magenta', lwd = 3)
lines(log(x[1:n1] ,20), log(washout_line_v[washout_line_v > 0], 20),
      col = 'cyan', lwd = 3)
