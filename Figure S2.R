L <- 6
v <- 0.5
k <- 0.1/(1/v)
D <- 0.2 # Assigning the parameters values
r <- 0.42
lamb <- (r*D)/(v^2)
xic <- (L*v)/D
alpha <- 6.13*(10^8)
phi0 <- 0 #  our object Numerical solution is phi
          #  phi0 is initial guess
          #  v is phi's 1st derivative

N <- 500 #Number of discreteization points

s <- seq(0,xic,length.out = N) # independent variable

func <- function(s, Y, pars){
  with (as.list(Y), {
    dphi=v
    dv=v+(lamb*phi*(1-phi))/(k+phi)   # our ODEs system
    return(list(c(dphi, dv)))
  })
}



bound <- function(i, Y, pars) {
  with (as.list(Y), {
    if(i==1) return (phi-v-1) #boundary conditions
    if(i==2) return (v)
  })
}
sguess <- seq(0, xic, length.out = N)
yguess <- matrix(ncol = N,
                 data = (rep(c(phi0,phi0-1), N)))
rownames(yguess) <- c("phi", "v")
Sol <-  bvpcol(x = s, func = func, bound = bound,
             xguess = sguess, yguess = yguess, #Solving
             leftbc = 1)
Sol[,2] <- Sol[,2]
Sol[,1] <- Sol[,1]*D/v
### I have to solve the solution using bvpSolve in order to draw our solution in Vector field


########## Plot phase plane ####

## function to draw vector field used by phaseR package
Func <- function(t, y, parameters) {
  x <- y[1]
  y <- y[2]
  
  dy <- numeric(2)
  dy[1] <- y
  dy[2] <- y+(lamb*x*(1-x))/(k+x)
  list(dy)
}

## Our eigenvectors directions
x1 <- c(-0.2,-0.147188,0,1) 
y1 <- 1.39509*x1
y2 <- -0.95080465*x1

### Solution vector field & adding other elements
Func.flowField <- flowField(Func, xlim = c(-0.3, 1.3), ylim = c(-0.25, 0.2), points = 17,xlab = "Phi",ylab="Chi", add = FALSE)

lines(Sol[,2], Sol[,3], type='l',lwd=2,col='blue') ## Solution (blue curve)

lines(x=seq(-0.47,1.5,len=70),y=seq(0,0,len=70),col='black') ## two vertical and horizental axis passing two singularity points
lines(x=seq(0,0,len=70),y=seq(-0.3,0.3,len=70),col='black')

lines(x1,y1,pch=20,lty='dashed',col='black',lwd=2) ## two eigenvecotrs direction
lines(x1,y2,pch=20,lty='dashed',col='red',lwd=2)

lines(x=c(1.2,1,0.8),y=c(0.2,0,-0.2),col='black',lwd=2) ## boundary line

points(0,0,lwd=4,pch = 20) ##two sigularity points
points(1,0,lwd=4,pch = 20)


#Func.nullclines <- nullclines(Func, ylim = c(-0.4, 0.4), xlim = c(-0.5, 1.50), points = 500)

