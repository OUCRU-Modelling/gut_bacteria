L=6
v=0.5
k=0.1/(1/v)
D=0.2 # Assigning the parameters values
r=0.42
lamb=(r*D)/(v^2)
xic=(L*v)/D
alpha <- 6.13*(10^8)
phi0=0 #our object Numerical solution is phi
#phi0 is initial guess
# v is phi's 1st derivative



N=50 #Number of discreteization points



s=seq(0,xic,length.out = N) # independent variable



musn <- function(s, Y, pars){
  with (as.list(Y), {
    dphi=v
    dv=v+(lamb*phi*(1-phi))/(k+phi)# our ODEs system
    return(list(c(dphi, dv)))
  })
}



bound <- function(i, Y, pars) {
  with (as.list(Y), {
    if(i==1) return (phi-v-1) #boundary conditions
    if(i==2) return (v)
  })
}
sguess = seq(0, xic, length.out = N)
yguess <- matrix(ncol = N,
                 data = (rep(c(phi0,phi0-1), N)))
rownames(yguess) <- c("phi", "v")



Sol = bvpcol(x = s, func = musn, bound = bound,
             xguess = sguess, yguess = yguess, #Solving
             leftbc = 1)
Sol[,2] <- Sol[,2]*(1/v)
Sol[,1] <- Sol[,1]*D/v

alp_v <- seq(1, 1, le=N)
Bacte <- alpha*(1/v)*alp_v - (alpha*Sol[,2])
Bacte <- Bacte*(10^-9)

ro <- r*(Sol[,2]/(k+Sol[,2]))
### draw 2 plot with 2 y-axis
par(mar = c(6, 5, 5, 7) + 0.05 )              # Additional space for second y-axis
plot(Sol[,1], Sol[,2], lwd=1.5, xlab="Position x (cm)", ylab="Food concentration F (mM) ", col = "green2",type ="l")              # Create first plot
par(new = TRUE)                             # Add new plot
plot(Sol[,1], Bacte , col = "blue",              # Create second plot without axes
     axes = FALSE, xlab = "",ylab="", lwd=1.5, type = "l")
axis(side = 4, at = pretty(range(Bacte)),col="blue")      # Add second axis

mtext("Bacteria B ( x10 ^9/mL)",col="blue", side = 4,line =2)
par(new=TRUE)
plot(Sol[,1],ro,col="red",axes=FALSE,xlab="", ylab="", lwd=1.5, type="l")
axis(side = 4, at = pretty(range(ro)),col="red", line = 4)
mtext("local reproduction p (/h)", side = 4, line = 6,col="red")
#plot(Sol,which="phi",col="green",lwd=1) #Plotting
#lines(Sol[,1],Bacte)