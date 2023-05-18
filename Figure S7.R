library(bvpSolve)
library(ReacTran)
library(deSolve)
L     <- 6
N     <- 211
v     <- 0.5
k     <- 0.1
F_in  <- 1/v
kappa <- 0.1/(1/v)
D     <- 0.4             # Assigning the parameters values
r     <- 0.42
lamb  <- (r*D)/(v^2) 
alpha <- 6.13*(10^8)


## Solve for initial condition (i.e the stationary state without mutant)


phi0  <- 0    # our object Numerical solution is phi
              # phi0 is initial guess
              # v is phi's 1st derivative

s     <- seq(0,xic,length.out = N) # independent variable

func <- function(s, Y, pars){
  with (as.list(Y), {
    dphi=v
    dv=v+(lamb*phi*(1-phi))/(kappa+phi) ## our ODEs system
    return(list(c(dphi, dv)))
  })
}

bound <- function(i, Y, pars) {
  with (as.list(Y), {
    if(i==1) return (phi-v-1) # boundary conditions
    if(i==2) return (v)
  })
}
sguess <-  seq(0, xic, length.out = N)
yguess <- matrix(ncol = N,
                 data = (rep(c(phi0,phi0-1), N)))
rownames(yguess) <- c("phi", "v")
Sol <-  bvpcol(x = s, func = func, bound = bound,
               xguess = sguess, yguess = yguess,    # Solving
               leftbc = 1)
Sol[,2] <- Sol[,2]
Sol[,1] <- Sol[,1]*D/v

#### Solve the system with mutant first appear at xm

vec_F_in  <- rep(1,N)
Food_star <- Sol[,2]*F_in                 ### The initial condition of PDE system
B_star    <- alpha*(vec_F_in - Food_star) ### for B and F
xm        <- 3                            ### the first location mutant appear
M0        <- B_star/2                     ### mutant concentration at this location (xm)

times     <- seq(0, 1,len=100)
xgrid     <- setup.grid.1D (x.up = 0, x.down = L, N = N) ## generating the gird for our solution
x         <- seq(0,6,len=N)                              ### Our discretization points
F_ini     <- Food_star
B_ini     <- B_star                                      ### Our initial condition of the system
M_ini     <- rep(1,N)


Mini      <- function (x){
        a <- abs(x-xm)
        if (a<=(0.01/2)){
          return(M0)
        }    ### Initial condition for Mutant
        else{
          return(0)
        }
        }
for(i in 1:N){
  M_ini[i] <- Mini(x[i]) ### Asigning the value of function "Mini" for initial condition of Mutant
}

yini <- c(F_ini, B_ini, M_ini) ## Initial condition to use in ode1D


### Solving the system using ode.1D

Sol_system <- function(t, Y, parms) {
  
     Food  <- Y[1:N]
     B     <- Y[(N+1):(2*N)]
     M     <- Y[((2*N)+1):(3*N)]
     dFood <- -(r/alpha)*(B+M)*Food/(k+Food) + tran.1D(C = Food, D = D, flux.up = v*F_in, flux.down = NULL, v=v, dx = xgrid)$dC
     dB    <- r*(B*Food)/(k+Food)            + tran.1D(C = B   , D = D, flux.up = 0     , flux.down = NULL, v=v, dx = xgrid)$dC   ###tran1D to describe divection diffusion equation
     dM    <- r*(M*Food)/(k+Food)            + tran.1D(C = M   , D = D, flux.up = 0     , flux.down = NULL, v=v, dx = xgrid)$dC
    
    return(list(c(dFood, dB, dM)))
}

print(system.time(
  out  <- ode.1D(y = yini, func = Sol_system, times = times,nspec = 3, names = c("Food","B","M"), parms = NULL , dimens = N)
))

a1     <- xm*rep(1,length(times)) +  sqrt(2*D*times) ## Different lines in the heatmap (from the article)
a2     <- xm*rep(1,length(times)) -  sqrt(2*D*times)
a3     <- xm*rep(1,length(times)) +  v*times

### Visualization

par(mfrow = c(2,2))
par(mar = c(3.7,3.7,4,4)+0.4)
image(out,legend=TRUE, xlab="times", ylab="position", grid = xgrid$x.mid, which = "B")
lines(times, a1, col="white")
lines(times, a2, col="white")
lines(times, a3, col="black")
