# Jack Lennon s2463109,Denis Zorba s2461643,Yeshwanth Zagabathuni s2319494
# ...Github link...
# Contributions:
# Yeshwanth: Modifications in newt() and hessian(), designing test cases 
# and testing and comments
# Our objective in this project is to optimize a given function using newtons 
# method. Newtons method is a numerical optimization method that begins with
# assuming a set of possible values for a function f(theta). Iteratively it
# proceeds to find an estimate of the optimal values of theta that minimize the
# function 'f' using the gradient of the function, its hessian and thus delta 
# and a few more paramters. But the algorithm also has to reach the optimum 
# theta by going through some conditions such as limiting the number of 
# iterations, limiting the number of times delta is halved, ensuring that the
# gradient does not tend towards infinity and so on. 


# First we need to find the hessian of the function. The function has following 
# parameters:
# theta = the variables of the function to find the Hessian of grad
# grad = the gradient of the function using eps 
# eps = the finite difference intervals to use when a Hessian function is not 
# provided.
# ... are just extra parameters to pass to grad
hessian <- function(theta, grad, eps, ...) {
  
# Start by putting grad into g to make it easier to see what's going on
  g<-grad
  
# Get the number of elements in theta
  n<-length(theta)
  
# Create an empty nxn matrix where we'll put the Hessian
  hess<-matrix(nrow=n,ncol=n)
  
# Loop over each row in the Hessian
  for (i in 1:n){
# Loop over each column in the current row
# This makes the lower triangle of the Hessian
    for (j in 1:i){
# Here rep() is just a function that replicates the value 0 'n' times to
# form a vector of size 'n' filled with zeros.
      e_i<-rep(0,n)
# Make the vector that will add a slight change to the ith element of theta
      e_i[i]<-1
      
# Here rep() is just a function that replicates the value 0 'n' times to
# form a vector of size 'n' filled with zeros.
      e_j<-rep(0,n)
# Make the vector that will add a slight change to the jth element of theta
      e_j[j]<-1
      
# This is the equation for the element H_ij
hess[i, j] <- ((g(theta + eps*e_j,...)[i] - g(theta-eps*e_j,...)[i])/(4*eps) +  
                ((g(theta+eps*e_i,...)[j] - g(theta-eps*e_i,...)[j])/(4*eps)))
      
# Since, the Hessian is symmetric, so we can skip calculating the upper 
# triangle and just use the transpose of the lower triangle 
# for the upper triangle.
      hess[j, i] <- hess[i, j]
    }
  }
# Finally we return the hessian.
  hess
}


# Now that we do have the hessian, we can thus pass it to the Newton's 
# optimization algorithm initialized in newt(). The paramters are as follows:
# theta - vector of initial values for the optimization parameters
# func - The objective function that we need to minimize
# grad - The gradient of the function which is basically a set of parial 
# derivatives of the function with respect to each theta variable.
# hess - The hessian matrix function that we created earlier.
# tol - The convergence tolerance
# fscale - a rough estimate of the magnitude of func near the optimum. This is
# useful at convergence testing.
# maxit - The maximum number of iterations before giving up.
# eps the finite difference intervals to use when a Hessian function 
# is not provided.
newt <- function(theta, func, grad, hess=NULL, ..., tol=1e-8, fscale=1, 
                 maxit=100, max.half=20,eps=1e-6) {
# The gradient of the function is computed using the grad() function with 
# respect to all variables in theta.
  g <- grad(theta,...)
# If the hessian has not been initialized, then is.null(hess) would return 
# true. 
  if (is.null(hess)){
# Thus we create the hessian using the hessian() function that we created  
# created earlier with theta, grad and tolerance value 'eps' as paramters.
    hess<-function(theta,...){hessian(theta,grad,eps,...)}
  }
  
# From here on, we will be using the functions stop() and warning() wherever
# necessary. stop() stops the execution of the rest of the code with a
# error message while warning() only displays the message. So based on our 
# requirement, we will be placing these functions where we desire to check
# the execution progresses as desired.
# If the objective function is not finite at initial theta, it is a
# red flag and thus there is no point of continuing and thus we stop 
# execution with a message.
  if(is.infinite(func(theta, ...)))
    stop('The objective function is not finite at the initial theta')
# If the gradient of the function ends up with values that are not finite,
# then once again there is no point of continuing and thus we stop 
# execution with a message.
  if (all(is.finite(g))==FALSE)
    stop('The gradient vector contains non-finite values at initial theta')
# Next we initialize a variable 'iterations' for counting how many times 
# we iterate Newton's method. This will later be used at the max iterations
# condition.
  iterations<-1

# Convergence can checked by whether all elements of the gradient vector have 
# absolute value less than 'tol' times the absolute value of the objective 
# function plus 'fscale' breaks when there is convergence.
# The opposite of this would be a greater then or equal to condition instead
# of less than and that would be: 
# length(g[abs(g) >= tol*(abs(func(theta, ...))+fscale)]) > 0
# Thus the loop would break if there is convergence.
  while (length(g[abs(g) >= tol*(abs(func(theta, ...))+fscale)]) > 0) { 
    
    iterations<-iterations+1
    
# We exceed the maximum number of iterations maxit, that means that we have 
# reached the limit to find the optimum. 
    if (iterations > maxit)
# And thus we terminate with stop().
      stop("Iterations exceeded maxit ", maxit)
    
# Next step is to check if the hessian is positive definite. If it isn't, 
# then we perturb it
    H<-hess(theta,...)

# If the cholesky decomposition of a matrix fails, then it is not positive
# definite. And thus it needs to be peturbed until it is.
while (inherits(try(chol(H), silent=TRUE), "try-error")) {
# Now we need to check if the hessian has any non-finite values and if it does
# we terminate with stop().
  if (all(is.finite(H))==FALSE)
        stop('The Hessian contains non-finite values')
# Perturb the hessian if it isn't positive definite
      H <- H + diag(length(theta))*(abs(min(eigen(H)$values))+1)
    }
# chol2inv(chol(H)) <-> solve(H) <-> H^-1 because H is positive definite
    Hi <- chol2inv(chol(H))
# Delta is calculated as H^-1*gradient(f) as follows:
    delta <- -Hi %*% grad(theta,...)
    
# step-halving
# we need to check that the function's value is actually being minimised.
# Thus we initialize a counter to count how many times delta is being halved.
# This is governed by a max.half condition.
counter <- 0
# "func(theta + delta,...) > func(theta,...)" is the condition that checks 
# whether the function is actually being minimized. At the same time,
# the function should not contain non-finite values. 
while(func(theta + delta,...) > 
      func(theta,...) || is.infinite(func(theta, ...))) {
# If max.half step halvings are exceeded we terminate with stop().  
  if (counter >= max.half) {
stop("Step halving fails to reduce the objective function despite 
     trying max.half step halvings")
      }
      
# Now half the step size
      delta <- delta/2
      
# Count the number of step halvings
      counter <- counter + 1
    }
    
# Update the value for theta
    theta <- theta + delta
# Now calculate the gradient again
    g <- grad(theta,...)
  }
# Thus we calculate the hessian and eventually check whether it is positive
# definite.
  H<-hess(theta,...)
# If 'H' is not positive definite , we return a warning and return the rest of 
# the values apart from the Hessian, that are:
#  The minimized function func
#  The optimium theta values (let's say "thetahat"), 
#  Number of iterations it took (iter) and 
#  value of gradient at thetahat (g).
  if (inherits(try(chol(H), silent=TRUE), "try-error")) {
warning("Hessian is not positive definite at minimum and 
        thus the inverse cannot be computed")
    Hi <- NULL
    return(list(f=func(theta, ...),
                theta=theta,
                iter=iterations,
                g=grad(theta, ...)))
  } else {
    Hi <- chol2inv(chol(H))
  }
# If 'H' is positive definite, we compute the inverse using 'chol2inv(chol(H))'
# and thus return inverse of Hessian, Hi as well.  
  return(list(f=func(theta, ...),
              theta=theta,
              iter=iterations,
              g=grad(theta, ...),
              Hi=Hi))
}
#test case 1
rb <- function(th,k) {
  k*(th[2]-th[1]^2)^2 + (1-th[1])^2
}
gb <- function(th,k) {
  c(-2*(1-th[1])-k*4*th[1]*(th[2]-th[1]^2),k*2*(th[2]-th[1]^2))
}
hb <- function(th,k) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]^2) - 4*th[1]^2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}

#Test case 2
dummy0<-function(theta){ # (3,0.5)
  x<-theta[1];y<-theta[2]
  return((1.5-x+x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625-x+x*y^3)^2)
}

dummy_grad0<-function(theta){
  return(grad(dummy,theta))
}
dummy_hess0<-function(theta){
  return(hessian_csd(dummy,theta))
}

#Test case 3
dummy1<-function(theta){  # x=(-1), y=1.5
  x<-theta[1];y<-theta[2]
  return(x-y+2*x^2+2*x*y+y^2)
}
dummy_grad1<-function(theta){
  return(grad(dummy1,theta))
}
dummy_hess1<-function(theta){
  return(hessian_csd(dummy1,theta)) # Hessian for complex step derivatives (csd)
}

#Test Case 4
dummy2<-function(theta){  # The Hessian contains non-finite values
  x<-theta[1];y<-theta[2]
  return((exp(x*y)-exp(-x*y))/2)
}
dummy_grad2<-function(theta){
  return(grad(dummy2,theta))
}
dummy_hess2<-function(theta){
  return(hessian_csd(dummy2,theta)) # Hessian for complex step derivatives (csd)
}

# Test Case 5
dummy3<-function(theta){ # This is negative definite
  x<-theta[1];y<-theta[2]
  return(-y-x^2+2*x-2)
}
dummy_grad3<-function(theta){
  return(grad(dummy3,theta))
}
dummy_hess3<-function(theta){
  return(hessian_csd(dummy3,theta)) # Hessian for complex step derivatives (csd)
}

# Test Case 6
dummy4<-function(theta){ # This is positive in-definite
  x<-theta[1];y<-theta[2]
  return(-y+2*x^2-x-5)
}
dummy_grad4<-function(theta){
  return(grad(dummy4,theta))
}
dummy_hess4<-function(theta){
  return(hessian_csd(dummy4,theta)) # Hessian for complex step derivatives (csd)
}

# Test Case 7
dummy5<-function(theta){ # Hess
  x<-theta[1];y<-theta[2];z<-theta[3]
  return(x^2*y^2*z^2+1)
}
dummy_grad5<-function(theta){
  return(grad(dummy5,theta))
}
dummy_hess5<-function(theta){
  return(hessian_csd(dummy5,theta)) # Hessian for complex step derivatives (csd)
}

#Test Case 8
dummy6<-function(theta){ # This is positive definite
  return((theta[1]+theta[2]+theta[3])^2-1)
}
dummy_grad6<-function(theta){
  return(grad(dummy6,theta))
}
dummy_hess6<-function(theta){
  return(hessian_csd(dummy6,theta)) # Hessian for complex step derivatives (csd)
}

#Test Case 9
dummy7<-function(theta){ # This is positive definite
  x<-theta[1];y<-theta[2];z<-theta[3]
  return(100*(x^2+y^2+z^2)+2*x*z+4*y*z)
}
dummy_grad7<-function(theta){
  return(grad(dummy7,theta))
}
dummy_hess7<-function(theta){
  return(hessian_csd(dummy7,theta)) # Hessian for complex step derivatives (csd)
}

#Test Case 9
dummy8<-function(theta){ # This is positive definite
  x<-theta[1];y<-theta[2];z<-theta[3];w<-theta[4]
  return((x-y)^2+(y-z)^2+(z-w)^2+(x-w)^2+x*y)
}
dummy_grad8<-function(theta){
  return(grad(dummy8,theta))
}
dummy_hess8<-function(theta){
  return(hessian_csd(dummy8,theta)) # Hessian for complex step derivatives (csd)
}


findgradient<-function(theta){
  trig.exp <- expression((1.5-x+x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625-x+x*y^3)^2)
  ( dxy <- deriv(trig.exp, c("x", "y")) )
  y <- theta[2]
  x<- theta[1]
  eval(dxy)
}

t80 <- 1:13 ## years since 1980
y1 <- c(12,14,33,50,67,74,123,141,165,204,253,246,240) ## AIDS cases
nll <- function(theta) {
  ## -ve log likelihood for AIDS model y_i ~ Poi(alpha*exp(beta*t_i))
  ## theta = (alpha,beta)
  mu <- theta[1] * exp(theta[2] * t80) ## mu = E(y)
  -sum(dpois(y1,mu,log=TRUE)) ## the negative log likelihood
} ## nll
gll <- function(theta) {
  ## grad of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1];beta <- theta[2] ## enhances readability
  ebt <- exp(beta*t80) ## avoid computing twice
  -c(sum(y1)/alpha - sum(ebt), ## -dl/dalpha
     sum(y1*t80) - alpha*sum(t80*ebt)) ## -dl/dbeta
} ## gll

# 2 Variable Functions
newt(c(400,1),dummy0,dummy_grad0, maxit=500)

newt(c(400,1),dummy1,dummy_grad1, maxit=500)

newt(c(400,1),dummy2,dummy_grad2, maxit=500)

newt(c(400,1),dummy3,dummy_grad3, maxit=500)

newt(c(400,1),dummy4,dummy_grad4, maxit=500)

# 3 Variable Functions
newt(c(10,10,10),dummy5,dummy_grad5, maxit=500)

newt(c(10,10,10),dummy6,dummy_grad6, maxit=500)

newt(c(10,10,10),dummy7,dummy_grad7, maxit=500)

f<-function(x){(x[1]+x[2]+x[3])^2-1}
p=(c(10,10,10))
nlm(f,p)

newt(c(0.01,0.01,0.01,0.01),dummy8,dummy_grad8, maxit=500)

f<-function(x){(x[1]-x[2])^2+(x[2]-x[3])^2+(x[3]-x[4])^2+(x[1]-x[4])^2+x[1]*x[2]}
p=(c(1,1,1,1))
nlm(f,p)