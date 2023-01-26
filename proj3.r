#Yeshwanth Zagabathuni - s2319494

# Our task is to estimate an unknown smoothing function f, using known basis 
# functions bj(x). To avoid overfitting, what we can do is impose a penalty and
# thus penalize the function to ensure smooth flow. Such functions are thus
# referred to as P-splines, where 'P' stands for Penalized. The model will be
# estimated using penalized least squares with unknown parameter lambda
# which we will eventually choose based on a criteria called the generalized
# cross validation criterion (gcv).

# First we define the pspline() function that will first work on computing the
# smoothing parameter lambda based on minimum gcv criterion and will then work
# its way to find required parameters such as the coefficient matrix(coeff), 
# fitted values(fitted), residual variance(sig2) and covariance 
# matrix(cov_matrix). 
# pspline is defined as follows:
# x and y – the vectors of x, y data to smooth.
# k – the number of basis functions to use.
# logsp – the ends of the interval over which to search for the smoothing 
# parameter (log λ scale). If only a single value is provided then no searching 
# is done, and the spline is returned for the given log λ value.
# bord – the B-spline order to use.
# pord – the order of difference to use in the penalty. 
# ngrid – the number of smoothing parameter values to try. 
pspline<-function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100) {
  # First we initialize 'n', which is number of rows in the dataset and it can be 
  # computed using length(x) or even length(y). 
  n<-length(x)
  # Next, we compute the difference 'dk' by which we need to increment every time
  # to get the 'knots'
  dk <- diff(range(x))/(k-bord) 
  # We initialize first value to min(x)-dk*bord. From there the seq() function 
  # increments by 'dk' every time to produce our 'knots'  
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  # Next, we compute the basis matrix 'X' with the 'knots' as computed above
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  # Compute 'D' matrix using the diff() with difference=2 and thus the diff() is
  # applied twice on the identity matrix of order 'k' to thus establish our 'D'
  # matrix.
  D <- diff(diag(k),differences=pord) 
  
  # Now our objective is to compute a list of 'ngrid' possible smoothing parameter 
  # (lambda) values and find the value that minimises gcv. These values should be
  # evenly spaced on the Log Scale.
  # For this, we first compute the difference that we need to add up every time 
  # starting from logsp[1] until we finally get to logsp[2]. This can be 
  # given by (logsp[2]-logsp[1])/(ngrid-1).
  d<-(logsp[2]-logsp[1])/(ngrid-1) 
  # Now we use the seq() function to generate a sequence of linearly spaced values 
  # with common difference 'd' and the result is stored in list 'lst'.
  lst<-list(seq(logsp[1],logsp[2],by=d))
  # Next we calulate exponents of those values to calculate 'ngrid' log-spaced 
  # values staring from logsp[1] to logsp[2]. This can be done easily and 
  # effectively by calling the sapply() function that applies exp() to all the 
  # linearly spaced values computed previously. The results are stored in x1.
  x1<-sapply(lst,function(i) exp(i)) 
  # Now we need to calculate the corresponding 'µhat' values for 'ngrid' lambda 
  # values in x1. We know that µ=X*β and that β=(XT*X + λ*DT*D)^(−1)*XT*y.
  # Thus we use sapply() to apply the formula on each value in list x1. 
  # Note that, t() is used to calculate transpose of a matrix and solve() for its
  # inverse.
  muhat_results<-sapply(x1,function(i) X%*%(solve(t(X)%*%X+i*t(D)%*%D)%*%t(X)%*%y)) 
  
  # From here, our objective is the effective computation of gcv values.
  # First 'X' undergoes QR decomposition. The qr() function is used for the same. 
  QR<-qr(X)
  # We are only interested in the 'R' matrix and thus we extract it into the 
  # variable 'R'.
  R<-qr.R(QR)
  # We need to compute R^(−T)*DT*D*R^(−1), with the result stored in 'Ei'.
  Ei<-solve(t(R))%*%t(D)%*%D%*%(solve(R))
  # Now we find eigen decomposition of the matrix 'Ei' = UΛUT
  E<-eigen(Ei)
  # We are only interested in the diagonal matrix Λ that stores a diagonal of 
  # eigen values and thus we create it (matrix 'V') using diag() function. 
  V<-diag(E$values) 
  # diag(k) creates an identity matrix of order 'k'. 
  I<-diag(k)
  # Now I and V are used to effectively compute 'ngrid' effective degrees of
  # freedom. κ = tr{(I + λi*Λ)^(-1)} for every λi in 'x1'.
  Result<-sapply(x1,function(i) sum(diag(solve(i*V+I))))
  # mu= ∥y − µhat∥^(2)/(n − κ ) and gcv= mu/(n − κ)and thus,
  # gcv =  ∥y − µhat∥^(2)/(n − κ)^(2) 
  # We iterate through every column vector in Gcv and compute 'ngrid' gcv values
  Gcv<-sapply(1:ngrid,function(i) sum((y-muhat_results[,i])^2)/(n-Result[i])^2)
  # Now the min() function is used to find minimum gcv and thus λ,κ values   
  gcv<-min(Gcv)
  # which.min(Gcv) points to the index of the λ value with minimum gcv.
  # x1[which.min(Gcv)] gives required λ value.
  lambda<-x1[which.min(Gcv)] 
  # Result[which.min(Gcv)] gives required κ value.
  K<-Result[which.min(Gcv)]
  # Now we calculate coefficient matrix 'coef' using βhat=(XT X + λDT D)^(−1)*XT*y 
  coef<-solve(t(X)%*%X+lambda*t(D)%*%D)%*%t(X)%*%y
  # Next we calculate µhat = X*βhat or fitted<-X%*%coef
  fitted<-X%*%coef
  # The Residual Variance is calculated below:
  # σˆ2 = ∥y − µˆ∥^2/(n − κ)
  sig2<-sum((y-fitted)^2)/(n-K)
  # The Covariance Matrix is calculated below: 
  # (XT*X + λ*DT*D)^(−1)*σˆ2
  cov_matrix<-solve(t(X)%*%X+lambda*t(D)%*%D)*sig2
  # We finally return required values for further operations in a list as follows:
  res_set<-list(coef=coef,
                fitted=fitted,
                sig2=sig2,
                K=K,
                gcv=gcv,
                bord=bord,
                pord=pord,
                k=k,
                knots=knots,
                n=n,
                x=x,
                y=y,
                cov_matrix=cov_matrix,
                X=X)
  # The class of the list is set of 'pspline'
  class(res_set)<-'pspline'
  # We return the list
  return(res_set)
}

x<-mcycle$times ####
y<-mcycle$accel ####
n<-length(x)    ####
m<-pspline(x,y) ####

# The next part is to create a 'print' function of class 'pspline'.
# For this function, we will be using values returned by the pspline() function 
# to compute required values. These values are in the object 'm'.
print.pspline<-function(m){
  # Extract 'gcv' returned by pspline().
  gcv<-m$gcv
  # Extract 'edf' returned by pspline().
  edf<-m$K
  # Compute r^2 as 1-((n-1)*sig2)/((y-mean(y))^2)
  r2<-1-((m$n-1)*m$sig2)/sum((m$y-mean(m$y))^2)
  # We will next use cat() function to print a formatted output as follows. This
  # output contains a bit of statistical information about the data.
  cat('Order' ,m$bord, 'p-spline with order',m$pord,'penalty')
  cat('\nEffective degrees of freedom:',m$K,' Coefficients:',m$k)
  cat('\nresidual std dev:',sqrt(m$sig2),'  r-squared:',r2,'  GCV:',m$gcv,'\n\n')
  # invisible() helps to return an invisible list containing gcv, edf and r2 
  invisible(list(gcv,edf,r2))
}
m2<-print(m)   ####
m2             ####

# The next part is to create a predict() function of class 'pspline'.
# Predictions from the smooth fit, for new x values within the range of 
# the original data are carried out.
# It has the following arguments:
# m - An object of class 'pspline'
# x - New 'x' values
# se - A boolean value that determines whether we return a list of 2 items
# or just the fitted values
predict.pspline<-function(m,x,se=TRUE){
  # We use the splineDesign() function with previous settings for new 'x' values
  # to generate 'Xp' matrix.
  Xp<-splines::splineDesign(m$knots,x,ord=m$bord+1,outer.ok=TRUE)
  # We find the new fitted values, fit=Xp*coef
  fit<-Xp%*%m$coef
  # If se=TRUE, then we are required to compute Standard Errors. Else we return
  # the vector of predictions.
  if(se){
    # The Standard Errors can be computed with least time complexity as follows:
    se<-rowSums(Xp*(Xp%*%m$cov_matrix))^.5
    # We create a list with first item as ypred 'fit' and second as Standard 
    # Errors 'se'. 
    res_lst<-list(fit=fit,se=se)
    # Return list
    return(res_lst)
  }
  # if se=FALSE we return the predicted 'y' values
  return(fit)
}
coords<-predict(m,x)    ####
class(coords)
# The plot() function under the 'pspline' class. It only takes one argument and
# that is the object 'm'.
plot.pspline<-function(m){
  # We need to define linearly sequenced values of 'x' and thus the increment is 
  # as follows
  d<-(max(m$x)-min(m$x))/(m$n-1)
  # We use seq() with increment of 'd' to generate linearly spaced 'x' values 
  x_new_values<-seq(min(m$x),max(m$x),by=d) 
  # Next we estimate y predictions and thus standard errors 
  predictions<-predict(m,x_new_values,se=TRUE)
  # The lower limit is defined as follows: fit-1.96*se
  ll<-predictions$fit-1.96*predictions$se
  # The upper limit is defined as follows: fit+1.96*se
  ul<-predictions$fit+1.96*predictions$se
  # Next, we plot original x,y data using plot() with labels and title as follows:
  plot(m$x,m$y,xlab='x',ylab='y',sub = 'Original Coordinates Plot')
  # Now we plot the smooth function overlaid as a line using lines() as follows:
  lines(x_new_values,predictions$fit)
  # The upper confidence limit is plotted as follows:
  lines(x_new_values,ul,col='blue')
  # The lower confidence limit is plotted as follows:
  lines(x_new_values,ll,col='red')
  
  # Now let us define a linear interaction between 'fitted' and 'x'
  lin_mod<-lm(m$fitted ~ m$x)
  # We plot linear model using plot(). The which=c(1,2) selects the QQplot and
  # residuals plots alone to be plotted
  plot(lin_mod,which=c(1,2))
  # Finally, we return an invisisble list with lower limit ll, upper limit ul
  # and linearly spaced 'x' values
  return(invisible(list(ll,ul,x_new_values)))
}
plot(m)                 ####

