rm(list = ls())
library(bvpSolve)
library(ReacTran)
library(deSolve)
library(rgl)
#library(plot3D)
library(shape)
L           <- 6
N           <- 900                       ### Spatial discretization points, must be odd in order to obtain the exact stationary solution
v           <- 0.2                       ### flow velocity of Food, Antibiotic, Bacteria and Mutant respectively
D           <- 0.4                        ### Diffusion coefficient of Food, Antibiotic, Bacteria and Mutant respectively
k           <- 0.1                         ### Monod constant
F_in        <- 1/v                      ### Food concentration at the entrance of the gut
A_in        <- 1/v                      ### Antibiotic concentration at the entrance of the gut
F1          <- rep(1, N)
F2          <- rep(1, N)
G1          <- rep(1, N)
G2          <- rep(1, N)
f1_0        <- 0
g1_0        <- 0
##### Assigning the parameters values #####
r           <- 0.42                        ### Growth rate of Bacteria 
alpha       <- 6.13*10^8                   ### Yield of Food to Bacteria and mutants respectively
beta        <- 6.13*10^7                   ### Consumption of antibiotic in killing Bacteria and mutants respectively
A_50        <- 0.09                        ### Concentration of Antibiotic corresponding to a half of elimination efficiency on Bacteria 
alpha1      <- (r*F_in)/(k+F_in)
alpha2      <- (A_in)/(A_50 + A_in)
delta_wo    <- (alpha1 - (v^2)/(3.2*D))*(1/alpha2)
delta_max   <- 0.2
tmax        <- 700
times       <- seq(0, tmax,len=200)                        ### discretization of times
xgrid       <- setup.grid.1D (x.up = 0, x.down = L, N = N) ### generating the gird for our solution
x           <- xgrid$x.mid                 ### We should cho x.mid rather than x.int ### Our discretization points
dx          <- N/(2*L)
col         <- c('red3','gray', 'red', 'coral2', 'cyan2', 'green', 
                 'brown', 'chartreuse3', 'blue', 'chartreuse3', 'coral', 'gray2', 'coral2')
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

# #### Our Solution #######
# bound <- function(i, Y, pars) {
#   with (as.list(Y), {
#     if(i==1) return (-D*f2+v*f1-1) #boundary conditions
#     if(i==2) return (-D*g2+v*g1-1)
#     if(i==3) return (f2)
#     if(i==4) return (g2)
#   })
# }
# xguess <- x
# yguess <- matrix(ncol = N,
#                  data = (rep(c(f1_0, -(1 - v*f1_0)/D, g1_0, -(1 - v*g1_0)/D), N)))
# rownames(yguess) <- c("f1", "f2", "g1", "g2")
# Sol <-  bvpcol(x = x, func = func, bound = bound,
#                xguess = xguess, yguess = yguess, #Solving
#                leftbc = 2, atol = 1e-9)
### Phase Plane ######
#par(mfrow=c(2,1))

#### Eigen and axes ######
### (0, 0, 0, 0) ###
A=t(matrix(c(0, 1, 0 , 0, ((1/k)*(1/D)*(r/alpha)*(F_in*(alpha-beta))), v/D, 0, 0, 
             0, 0, 0, 1, 0, 0, ((1/A_50)*(1/D)*(delta_max/beta)*(F_in)*(alpha-beta)), v/D)
           , ncol=4, nrow = 4))
eigen((A))
v1=eigen((A))$values[1]
v2=eigen((A))$values[2]
v3=eigen((A))$values[3]
v4=eigen((A))$values[4]
x1 <- seq(-8, 8, by=2) #generated x corrdinate
x2 <- seq(0, -20, by=-2)
y1 <- v1*x1     #the slopes is eigenvalues
y2 <- v2*x1
y3 <- v3*x1
y4 <- v4*x1
y0 <- -(1-v*x1)/D

par(mar = c(4, 4.7, 3, 3.5) + 0.05 )
a1   <- c(0.2,      5.68,      0.2,   5.43,   1.35,   1.45,  1.65,   1/v,  5.2,  1/v,    5.9,  4.8) 
a2   <- c(-2.45,   -0.5,   0.45,   -0.7,   -2.45,  -2.45, -2.45,  0.34, -0.8, -1.35,  -1,  -0.4)
a3   <- c(0.2,      5.68,      5.8,   5.43,   1.35,   1.45,  1.65,   1/v,  5.2,  4.18,  5.9,  4.8)
a4   <- c(-2.45,   -0.5,   0.45            ,   -0.7,   -2.45,  -2.45, -2.45,  0.34, -0.8, -1.35,  -1,  -0.4)

### F1 & F2 ####
for(j1 in 2) {
  yini <- c(f1=a1[j1], f2=a2[j1], g1=a3[j1], g2=a4[j1])
  out  <- ode(y=yini, times=x, func=func, atol = 1e-9)
  F1 <- out[,2]
  F2 <- out[,3]
  if(j1==1){
    plot(F1, F2, type='l', col='red3', lwd=3, xlim=c(-1.1, 6.2), ylim=c(-2.5, 1.6))
    lines(x=Sol[, 2], y=Sol[, 3], type='l', lwd=3, col='red', xlab='F1', ylab='F2',
        xlim=c(0, 1.2), ylim=c(-2.5, 3))
    lines(x1,y1,pch=20,lty='dashed',col='black',lwd=2.)
    lines(x1,y2,pch=20,lty='dashed',col='blue',lwd=2.)
    lines(x1,y0,pch=20,type='l',col='black',lwd=2)
    lines(x=seq(0, 0, len=15),y=seq(-6.5,6.5,len=15),col='black', lwd=1.5)
    lines(y=seq(0, 0, len=15),x=seq(-6.5,6.5,len=15),col='black', lwd=1.5)
    points(x=1/v, y=0,pch=20, col='black', lwd=4)
    points(0, 0,pch=20, lwd=4, col='black')}
    lines(F1, F2, type='l', lwd=2, col=col[j1],xlim=c(-3, 8), ylim=c(-3, 3))
}

# ### G1 & G2 ###
# for(j2 in 1) {
#   yini <- c(f1=a1[j2], f2=a2[j2], g1=a3[j2], g2=a4[j2])
#   out  <- ode(y=yini, times=x, func=func, atol = 1e-4)
#   G1 <- out[,4]
#   G2 <- out[,5]
#   if(j2==1){
#     plot(G1, G2, type='l', col='red3', lwd=3, xlim=c(-1.1, 6.2), ylim=c(-2.5, 1.6))
#     }
#     lines(G1, G2, type='l', lwd=2, col=col[j2],xlim=c(-3, 8), ylim=c(-3, 3))
# }


lines(x=Sol[,2], y=Sol[,3], type='l', lwd=3, col='red', xlab='G2', ylab='G1',
      xlim=c(0, 1.2), ylim=c(-2.5, 3))



lines(x1,y3,pch=20,lty='dashed',col='black',lwd=2.)
lines(x1,y4,pch=20,lty='dashed',col='blue',lwd=2.)
lines(x1,y0,pch=20,type='l',col='black',lwd=2)
lines(x=seq(0, 0, len=15),y=seq(-6.5,6.5,len=15),col='black', lwd=1.5)
lines(y=seq(0, 0, len=15),x=seq(-6.5,6.5,len=15),col='black', lwd=1.5)
points(x=1/v, y=0,pch=20, col='black', lwd=3)
points(0, 0,pch=20, lwd=4, col='black')


### (1/v, 0, 1/v, 0) ###
  
B=t(matrix(c(0, 1, 0 , 0, ((-r/v)*(1/D)*(k + 1/v)), v/D, (1/D)*(r/alpha)*(beta/(k + 1/v))*(1/(v^2)),0
             , 0, 0, 0, 1, 
             ((1/(A_50 + 1/v))*(1/D)*(delta_max/beta)*(-alpha)*(1/(v^2))), 0,
             ((1/(A_50 + 1/v))*(1/D)*(delta_max)*(1/v)),  v/D), ncol=4, nrow = 4))
eigen((B))
u1=eigen((B))$values[1]
u2=eigen((B))$values[2]
u3=eigen((B))$values[3]
u4=eigen((B))$values[4]
  
  
  #######################################################

# arrows3D (1, 2, 3, x1 = 1, y1 = 2, z1 = 2,   
#           colvar = NULL, phi = 40, theta = 40,
#           col = NULL, NAcol = "white", breaks = NULL,
#           colkey = NULL, panel.first = NULL,
#           clim = NULL, clab = NULL, bty = "b", type = "triangle", 
#           add = TRUE, plot = TRUE)


