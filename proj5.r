# Yeshwanth Zagabathuni - s2319494

# The task is to implement a demographic model for England 
# and Wales and to predict the expected number of deaths per week from the 
# beginning of 2020 if death rates had stayed at the same levels as in 2017-19. 
# The model will include an adjustment to allow for seasonal variation in  
# mortality rates. We will then compare this to the actual deaths per week that 
# occurred over this period to obtain an excess deaths time series, and then 
# model this series using a simple Bayesian model in JAGS. This will be used to
# generate 10000 samples for further analysis and visualization.


# First we code the predict_deaths() function that takes the following as input:
# mpop - male starting population
# fpop - female starting population
# mf - male annual death rates for each 1-year age band based on 2017-19 data
# mm - female annual death rates for each 1-year age band based on 2017-19 data 
# d - mortality rate identifier for each week 
# This function returns a vector of values that contain predicted deaths by week

predict_deaths<-function(mpop,fpop,mf,mm,d){
  # First we initialize the result vector
  vec<-c()
  # Here we loop over length(d) weeks and calculate deaths by week
  for(j in 1:length(d)){ 
    # Initialize the counters for male and female deaths
    male_D<-0
    female_D<-0
    # Initialize the Ni* vectors for Male and Female respectively
    mpop_star<-c()
    fpop_star<-c()
    # Now let us loop over the entire set of age classes
    for(i in 1:length(fpop)){ 
      # Male Deaths (Di) value for age class i
      Deaths_male_i <- 0.9885*d[j]*(1-exp(-mm[i]/52))*mpop[i] 
      # Female Deaths (Di) value for age class i
      Deaths_female_i <- 0.9885*d[j]*(1-exp(-mf[i]/52))*fpop[i] 
      # Ni* = Ni-Di for class i for males and then append to vector mpop_star
      mpop_star<-append(mpop_star, mpop[i]-Deaths_male_i) 
      # Ni* = Ni-Di for class i for males and then append to vector fpop_star
      fpop_star<-append(fpop_star, fpop[i]-Deaths_female_i) 
      
      if(i>1){
        # If i > 1 then the female and male populations are updated as 
        # Ni+ = Ni*51/52 + Ni-1*/52
        fpop[i]<-(fpop_star[i]*51+fpop_star[i-1])/52 
        mpop[i]<-(mpop_star[i]*51+mpop_star[i-1])/52 
      }
      # Sum the death value for male for week i for all age classes
      male_D <- male_D+Deaths_male_i   
      # Sum the death values for female for week i for all age classes
      female_D <- female_D+Deaths_female_i 
    }
    
    # Append total predicted deaths (male deaths+female deaths) for week i 
    vec<-append(vec,(male_D+female_D))
    
  }
  # Return vector of predicted deaths
  return(vec)
}

# Now we use the read.table() function to load the .dat files, death1722uk.dat
# and lt1720uk.dat respectively
a<-read.table('death1722uk.dat',header = TRUE)
b<-read.table('lt1720uk.dat',header = TRUE)

# Now let us run the predict_deaths() function from 2020 to end of data
# We pass the week identifier (d) values from 157 to the end as 2020 starts 
# from week 157
res_2020<-predict_deaths(b[['mpop20']],b[['fpop20']],
                         b[['mf']],b[['mm']],a[['d']][157:length(a[['d']])])

# Now let us find out the difference between the total actual deaths and the 
# total predicted deaths from the start of 2020 to the end of the data using
# the sum() function
Excess_deaths<-sum(a[['deaths']][157:length(a[['d']])])-sum(res_2020)
# Print the excess deaths from 2020 to the end of data 
print(Excess_deaths)

# Now let us initialize the weeks
x<-157:305
# We then plot the observed deaths against week with the title displaying the 
# overall Excess deaths from 2020 to end of data using plot() as follows:
plot(x,a[['deaths']][157:length(a[['d']])],
     xlab = 'Week',ylab = 'Observed deaths',
     main=paste('Excess deaths from 2020 to end of data is:',
                Excess_deaths))
# We then overlay a continuous curve showing the predicted deaths using the 
# lines() function
lines(x,res_2020)

# Now let us store a vector of values of excess deaths each week from 2020 to
# end of data
Excess_deaths_2020<-a[['deaths']][157:length(a[['d']])]-res_2020
# Initialize weeks
x<-157:305
# We then use the cumsum() function to calculate cumulative excess 
# deaths by week. Then we plot cumulative excess deaths by week: 
plot(x,cumsum(Excess_deaths_2020),xlab='Week',ylab='Cummulative excess deaths',
     main='Cummulative excess deaths by week')

# The values in weeks 51, 52, 53, 105 and 106 are excluded from distribution 
# due to various recording problems and are all replaced with NA values. But
# before we do let us store those values in a vector 'Excess_deaths_2020_nas'
# for plotting purposes later on.
Excess_deaths_2020_nas<-c(Excess_deaths_2020[51],Excess_deaths_2020[52],
                          Excess_deaths_2020[53],Excess_deaths_2020[105],
                          Excess_deaths_2020[106])
# The corresponding weeks are stored in 'xnas' vector
xnas<-c(x[51],x[52],x[53],x[105],x[106])
# Now let's fill those cells with 'NA' values 
Excess_deaths_2020[51]<-NA
Excess_deaths_2020[52]<-NA
Excess_deaths_2020[53]<-NA
Excess_deaths_2020[105]<-NA
Excess_deaths_2020[106]<-NA

# Now our next objective is to model our data in jags. The jags file model.jags
# contains the model initialization to do the sampling. We use jags.model() 
# function from jags package and pass the jags file, the list of excess deaths 
# (Excess_deaths_2020) and length of the vector to the sampler.
mod1<-jags.model('model.jags',
                 data=list(x=Excess_deaths_2020,N=length(Excess_deaths_2020)))
# We then use the coda.samples function to generate 10000 samples, monitoring
# mu, rho and k.
sam.coda<-coda.samples(mod1,c('mu','rho','k'),n.iter=10000)

# Now let us generate the trace plot for rho using plot() 
plot(sam.coda[[1]][,'rho']) 
# We next insert these 'rho' values into a dataframe rho_plot using the 
# data.frame() function.
rho_plot<-data.frame(rho=sam.coda[[1]][,'rho'])
# The dataframe is then passed to the ggplot() function 
histogra<-ggplot(rho_plot,aes(x=var1))
# Now we add the modifier function geom_histogram() to plot the required 
# histogram. The bars are coloured in black, the background is white and the
# bin-width is specified as 0.01.
histogra<-histogra+geom_histogram(binwidth = 0.01,color="white", fill="black")
# Histogram plot
print(histogra) 

# Now let us generate posterior expected value vector for "mu"
# We can do this by simply computing the mean for every mui list of sample 
# values generated by the sampler previously (We ignore mu[150] as it 
# is unnecessary). This can be the done very quickly using the sapply() function
posterior_mus<-sapply(2:150,function(i)mean(sam.coda[[1]][,i]))
# Let us view the expected posteriors
print(posterior_mus)

# Now we are required to plot every 50th sampled mu vector against week as a 
# grey curve. Next we need to plot expectation for mu overlaid in blue. Finally
# we will add excess deaths as black points and plot the values that we 
# changed to NA previously before the jags modelling as red points. 

# First we replicate weeks for 10000/50 times as there are that many 50th
# sampled mu vectors (50*1 to 50*200 and thus 200 vectors) and each vector has
# to be plotted against week
x<-rep(157:305,200)

# Now let us initialize an iterator that stores indexes of all 50th mu vectors 
# This can easily be done by generating a sequence starting with 50
# and ending with 10000 (50*1 to 50*200) and common difference 50 using 
# seq() function.
iterational<-seq(50,10000,by=50)

# Our task now is to compute a vector of all 50th sampled vector values: 
# The sapply() with the iteration values computed previously can make things
# easier and faster. We use it to quickly convert the coda lists to vectors
m50_s<-sapply(iterational,function(i)unlist(sam.coda[[1]][i,2:150]))

# Now let us plot all the 50th sampled mu values by week. The xlabel, ylabel 
# and main title are set using plot as follows:
plot(x,m50_s,xlab = 'Week',ylab ='50th sampled Muis',
     main='50th sampled Muis by Week')

# Now let us plot a curve coloured in grey using the lines() function and 
# setting col parameter to "grey"
lines(x,m2,col='grey') 

# Let's initialize total weeks from 2020 and beyond
x<-157:305
# Now let us plot the estimated expectation for 'mu' overlaid in blue using 
# the lines() function again setting col='blue'. The line width (lwd) is also
# set to 2 for a much thicker curve and better visibility.
lines(x,posterior_mus,col='blue',lwd=2)

# Next we plot observed Excess deaths by week. 
# This time we will be using the points() function to plot points as 
# black symbols. This can be done by specifying col='black' and the kind of
# symbol we want can be customized by the value of 'pch'. For 21 we get filled
# circle. Finally, line width is set to 2 for better visualization of points.
points(x,Excess_deaths_2020,pch=21,lwd=2,col='black')


# Now let us also visualize the error values that we stored in 
# Excess_deaths_2020_nas previously using points() and specifying col='red'
points(xnas,Excess_deaths_2020_nas,pch=21,lwd=3,col='red')

# Initialize weeks:
x<-157:305
# Plot xi-E(mui) by week using plot() with xlabel = Week, y label = residuals
# and main title being residuals by week.
plot(x,Excess_deaths_2020-posterior_mus,xlab = 'Week',ylab = 'residuals',
     main = 'Residuals by week')