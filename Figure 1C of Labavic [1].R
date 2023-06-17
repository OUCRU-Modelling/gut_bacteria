#rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(reshape2)
library(fields)
L         <- 6                                               ### Assigning the parameters values
r         <- 0.42
alpha     <- 6.13*(10^8)
k         <- 0.1
v         <- c(seq(0.05, 1, len=40),  seq(1.1, 5, len=10))     ### v and D scale have diffrent step sizes in it (Focus on Figure C, you can see the different leaps for diffrent interal of values)
F_in      <- 1/v
D         <- c(seq(0.05, 10, len=40), seq(11, 50, len=10))        
tmax      <- 500                                             ### I choose tmax=500 because in the article, they say when our system pass t=500hrs it reached the stationary state
N         <- 6000                                            ### N must be very large in order to abtain the exact solution compared to stationary system, but I still can put it 6000 because the consumption time problem
i         <- 1
j         <- 1
conc_profile  <- matrix(rep(1, len = length(v)*length(D)), ncol=length(v), nrow=length(D))
while(i<=length(D)){                                         ### creating the loop to draw the heatmap
  for(j in 1:length(v)){
    vec_F_in  <- rep(1,N)*F_in[j]               
    times     <- seq(0, tmax,len=100)                        ### discretization of times
    xgrid     <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
    x         <- xgrid$x.int                                 ### Our discretization points
    F_ini     <- vec_F_in*0.9
    B_ini     <- 0.1*alpha*vec_F_in
    Sol_system <- function(t, Y, parms) {
      Food  <- Y[1:N]
      B     <- Y[(N+1):(2*N)]
      dFood <- -(r/alpha)*B*Food/(k+Food) + tran.1D(C = Food, D = D[i], flux.up = 1   , flux.down = NULL, v=v[j], dx = xgrid)$dC
      dB    <- r*B*Food/(k+Food)          + tran.1D(C = B   , D = D[i], flux.up = 0   , flux.down = NULL, v=v[j], dx = xgrid)$dC
      return(list(c(dFood, dB)))
    }
    yini <- c(F_ini, B_ini)                                                                           
    print(system.time(
      out  <- ode.1D(y = yini, func = Sol_system, times = times, nspec = 2, names = c("Food","B"), parms = NULL , dimens = N)
    ))
    #j <- j+1
    conc_profile[i,j] <- (out[length(times),2]-out[length(times),N+1])/F_in[j]
  }
  i <- i+1
}

image.plot(t(conc_profile), x=v, y=D, ylim=c(0.05, 15), xlim=c(0.05, 1.4)) ### two ways to draw heatmap, but this image.plot using 'fields' package is more visualize
#heatmap(conc_profile, Rowv=NA, Colv=NA, col = heat.colors(256))