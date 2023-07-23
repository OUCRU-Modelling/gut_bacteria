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
v           <- 0.9 
v_B         <- 1.2
D           <- c(seq(0.01, 0.10, len=30), seq(0.11, 0.50, len=30), seq(0.501, 1, len = 25),
                 seq(1.01, 5, len = 25), seq(5.01, 10, len = 13), seq(11, 50, len = 25))        
D_B         <- c(seq(0.01, 0.10, len=30), seq(0.12, 0.50, len=30), seq(0.501, 0.91, len = 25),
                 seq(1.0, 5, len = 25), seq(5.01, 10, len = 13), seq(11, 50, len = 25))

delta_max   <- seq(0.01, 0.4, len=125)
numCores    <- detectCores()      ### Numbers of cores to be using parallel
##### Assigning the parameters values #####
k           <- 0.1                   ### Monod constant
r           <- 0.42                  ### Growth rate of Bacteria 
alpha       <- 6.13*10^8             ### Yield of Food to Bacteria 
beta        <- 6.13*10^7
A_50        <- 0.9
F_in        <- 1/v                
A_in        <- 1/v
alpha1      <- (r*F_in)/(k+F_in)
alpha2      <- (A_in)/(A_50 + A_in)
tmax        <- 700
ite         <- vector('list', length(delta_max)*length(D)) ### Creating the list index to be using in function in parallel
l <- 1
for(i in 1 : length(D)){
  for(j in 1 : length(delta_max)){
    ite[[l]] <- c(i, j)
    l        <- l + 1
  }
}
index <- seq(1, length(delta_max)*length(D), by = 1)       ### Index in parallel computation

Heat_map <- function(n){                           ### Main function to be compiled in parallel computation
  a           <- array(unlist(ite[n]), dim = c(1,  2))
  i           <- a[1]
  j           <- a[2]
  times       <- seq(0, tmax,len = 180)                         ### discretization of times
  xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N)    ### generating the gird for our solution
  x           <- xgrid$x.mid                                    ### We should cho x.mid rather than x.int 
  F_ini       <- rep(1,N)*F_in*0.9                              ### Initial condtitions
  A_ini       <- rep(1,N)*A_in
  B_ini       <- rep(1,N)*0.1*(alpha*F_in)
  
  Sol_system  <- function(t, Y, parms) {
    Food      <- Y[1:N]
    A         <- Y[(N+1):(2*N)]
    B         <- Y[((2*N)+1):(3*N)]
    dFood     <- tran.1D(C = Food, D = D[i],   flux.up = 1, flux.down = NULL, v = v, dx = xgrid)$dC + -(r/alpha)*B*Food/(k+Food)       
    dA        <- tran.1D(C = A   , D = D[i],   flux.up = 1, flux.down = NULL, v = v, dx = xgrid)$dC + -(delta_max[j]/beta)*B*(A)/(A+A_50) 
    dB        <- tran.1D(C = B   , D = D_B[i], flux.up = 0, flux.down = NULL, v = v_B, dx = xgrid)$dC + r*B*Food/(k+Food) + -delta_max[j]*B*(A)/(A+A_50)                                          ### tran1D to describe divection diffusion equation
    
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

washout_line_D           <- (v_B^2)/((alpha1 - delta_max*alpha2)*(2.3)) 
washout_line_delta       <- (alpha1 - v_B/L)*(1/alpha2)
conc_profile             <- array(unlist(results), dim = c (3, length(delta_max)*length(D)))
conc_profile_Food        <- matrix(conc_profile[1, ], nrow = length(delta_max), ncol = length(D))
conc_profile_Antibiotic  <- matrix(conc_profile[2, ], nrow = length(delta_max), ncol = length(D))
conc_profile_Bacteria    <- matrix(conc_profile[3, ], nrow = length(delta_max), ncol = length(D))
x  <- delta_max
y  <- D
# n1 <- length(washout_line_delta_1[washout_line_delta_1 > 0])
# n2 <- length(washout_line_delta_2[washout_line_delta_2 > 0])
# n3 <- length(washout_line_delta_3[washout_line_delta_3 > 0])

#### Heatmap Food ####
par(mar = c(4, 5, 3.5, 2) + 0.05 )
image.plot(conc_profile_Food, x=log(x, 20), y=log(y, 20),
           xlab = 'delta max', ylab = 'Diffusion Coefficient',
           main = 'Spatial dependence regarding 
           Food average, v=0.9, v_B=1.2
           ', col = NULL, axes = TRUE)
contour(z = conc_profile_Food, x=log(x, 20), y=log(y,20), add = TRUE) 
lines(log(x, 20), log(washout_line_D, 20), type = 'l'
      , col = 'magenta', lwd = 3)
lines(log(y ,20), rep(1, length(D))*log(washout_line_delta, 20),
      col = 'cyan', lwd = 3)

#### Heatmap drug ####
image.plot(conc_profile_Antibiotic, x=log(x, 20), y=log(y, 20), xlab = 'delta max'
           , ylab = 'Diffusion Coefficient',
           main = 'Spatial dependence regarding 
           Antibiotic average, v=0.9, v_B=1.2
           ', col = NULL, axes = TRUE)
contour(z = conc_profile_Antibiotic, x=log(x, 20), y=log(y, 20),  add=TRUE)
lines(log(x, 20), log(washout_line_D, 20), type = 'l'
      , col = 'magenta', lwd = 3)
lines(log(y ,20), rep(1, length(D))*log(washout_line_delta_1, 20),
      col = 'cyan', lwd = 3)

### Heatmap Bacteria ####
image.plot(conc_profile_Bacteria, x=log(x, 20), y=log(y, 20), xlab = 'delta max'
           , ylab = 'delta max',
           main = 'Spatial dependence regarding
           Bacteria average,  v=0.9, v_B=1.2
           ', col = NULL, axes=FALSE)
contour(z = conc_profile_Food, x=log(x, 20), y=log(y, 20),  add=TRUE)
lines(log(x, 20), log(washout_line_D, 20), type = 'l'
      , col = 'magenta', lwd = 3)
lines(log(y ,20), rep(1, length(D))*log(washout_line_delta_1, 20),
      col = 'cyan', lwd = 3)
