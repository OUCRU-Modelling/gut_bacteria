rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
library(rgl)
library(plot3D)
library(shape)
L           <- 6
N           <- 900                       ### Spatial discretization points, must be odd in order to obtain the exact stationary solution
v           <- 0.7                       ### flow velocity of Food, Antibiotic, Bacteria and Mutant respectively
D           <- 1.2                        ### Diffusion coefficient of Food, Antibiotic, Bacteria and Mutant respectively
k           <- 0.1                         ### Monod constant
F_in        <- 1/v                      ### Food concentration at the entrance of the gut
A_in        <- 1/v                      ### Antibiotic concentration at the entrance of the gut
F1          <- rep(1, N)
F2          <- rep(1, N)
G1          <- rep(1, N)
G2          <- rep(1, N)
f1_0 <- 0
g1_0 <- 0
##### Assigning the parameters values #####
r           <- 0.42                        ### Growth rate of Bacteria 
alpha       <- 6.13*10^8                   ### Yield of Food to Bacteria and mutants respectively
beta        <- 6.13*10^7                   ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- 0.09                        ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria 
alpha1      <- (r*F_in)/(k+F_in)
alpha2      <- (A_in)/(A_50 + A_in)
delta_wo    <- (alpha1 - (v^2)/(3.2*D))*(1/alpha2)
delta_max   <- 0.23
tmax        <- 700
times       <- seq(0, tmax,len=200)                        ### discretization of times
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x           <- xgrid$x.mid                 ### We should cho x.mid rather than x.int ### Our discretization points
dx          <- N/(2*L)
col         <- c('red3','orange', 'coral', 'blue', 'cyan', 'green', 'brown', 'magenta', 'yellow')
############################
func <- function(t, Y, parms){
  with (as.list(Y), {
    df1 <-  f2
    df2 <- (v*f2 + (r/alpha)*f1*(alpha*(F_in - f1) - beta*(F_in - g1))/(k + f1))/D   # our ODEs system
    dg1 <- g2
    dg2 <- (v*g2 + (delta_max/beta)*g1*(alpha*(A_in - f1) - beta*(A_in - g1))/(A_50 + g1))/D
    return(list(c(df1, df2, dg1, dg2)))
  })
}

#points3d(x=a1, y=0, z=0, col='cyan', size=8)

a1   <- c(0.000001,    0.9,   0.7,   0,  0.3, -0, 2.2, 2.3, 3, 3.1) 
a2   <- c(0, -1.7, -1.5, -0,  -1.1, -00, 0.1, 0.2, 0.5, 0.6)
a3   <- c(0.000001,  0.7,   0.47,  0.1, 0.05, 0, 0.22, 0.55, 0.6)
a4   <- c(0, -1.7, -1.5, 1, -0.9, 0.2, 2.22, 2.5, 2.7)

for(j in 1) {
yini <- c(f1=a1[j], f2=a2[j], g1=a3[j], g2=a4[j])
out  <- ode(y=yini, times=x, func=func, atol = 1e-2)
F1 <- out[,2]
F2 <- out[,3]
G1 <- out[,4]
G2 <- out[,5]

a11 <- (alpha-beta)/(v*alpha)
b11 <- (beta-alpha)/(v*beta)
if(j==1){
plot3d(x=F1, y=F2, z=G1, lwd=2, col=col[j])
plot(F1, F2, type='l', col='red3', lwd=2 
     )}#,xlim=c(0, 1.2), ylim=c(-0.35, 2.5))}
lines3d(x=F1, y=F2, z=G1, lwd=2, col=col[j])
lines(F1, F2, type='l', lwd=2, col=col[j])
lines(0, 0, type='p', lwd=2, col='red')
}

#### Our Solution #######
bound <- function(i, Y, pars) {
  with (as.list(Y), {
    if(i==1) return (-D*f2+v*f1-1) #boundary conditions
    if(i==2) return (-D*g2+v*g1-1)
    if(i==3) return (f2)
    if(i==4) return (g2)
  })
}
xguess <- x
yguess <- matrix(ncol = N,
                 data = (rep(c(f1_0, -(1 - v*f1_0)/D, g1_0, -(1 - v*g1_0)/D), N)))
rownames(yguess) <- c("f1", "f2", "g1", "g2")
Sol <-  bvpcol(x = x, func = func, bound = bound,
               xguess = xguess, yguess = yguess, #Solving
               leftbc = 2, atol = 1e-4)
lines(x=Sol[, 2], y=Sol[, 3],type='l', lwd=2, col='red')
points(x=5, y=0, col='cyan', lwd=8)


  

  
  
  #######################################################

# arrows3D (1, 2, 3, x1 = 1, y1 = 2, z1 = 2,   
#           colvar = NULL, phi = 40, theta = 40,
#           col = NULL, NAcol = "white", breaks = NULL,
#           colkey = NULL, panel.first = NULL,
#           clim = NULL, clab = NULL, bty = "b", type = "triangle", 
#           add = TRUE, plot = TRUE)


