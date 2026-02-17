
# Simultaneously simulate $n$ sample paths
# from a 1-dimensional OU process on the time interval $[0, T]$
# with time-step size $\Delta t$ 
# for user-supplied values of $(n, T, \Delta t, \mu, \nu, \sigma)$

simulate_OU = function(paths, t_end, delta_t, mu, nu, sigma){
  steps = t_end/delta_t
  increments = matrix(rnorm(steps*paths, mean=0, sd=sqrt(delta_t)), nrow=steps, ncol=paths)
  
  X_matrix = matrix(0, nrow=steps+1, ncol=paths)
  
  for (i in (2:(steps+1))){
    X_matrix[i,] = X_matrix[i-1,] + nu*(mu-X_matrix[i-1,])*delta_t + sigma*increments[i-1,]
  }
  return(X_matrix)
}



# plot the transition density of the 1-dimensional OU process
# for user-supplied values of $(n, s, t, x_s, \mu, \nu, \sigma)$

plot_transition_density_OU = function(paths, time_s, time_t, value_s, mu, nu, sigma){
  time_interval = time_t - time_s
  delta_t = 0.01 # default setting
  steps = as.integer(time_interval/delta_t)
  increments = matrix(rnorm(steps*paths, mean=0, sd=sqrt(delta_t)), nrow=steps, ncol=paths)
  
  X_matrix = matrix(value_s, nrow=steps+1, ncol=paths)
  
  for (i in (2:(steps+1))){
    X_matrix[i,] = X_matrix[i-1,] + nu*(mu-X_matrix[i-1,])*delta_t + sigma*increments[i-1,]
  }
  transition_density = X_matrix[nrow(X_matrix),]
  plot(density(transition_density))
}



# Make a plot of the standard normal density
# and superimpose lines for the Hermite polynomial expansion
# with $1, 2, ..., 5$ terms


H0 = function(x){
  return(rep(1, length(x)))
}
H2 = function(x){
  return(x^2 - 1)
}
H4 = function(x){
  return(x^4 - 6*x^2 + 3)
}

H6 = function(x){
  return(x^6 - 15*x^4 + 45*x^2 - 15)
}

H8 = function(x){
  return(x^8 - 28*x^6 + 210*x^4 - 420*x^2 + 105)
}

approx1 = function(x){
  return(1/(2*sqrt(pi))*H0(x))
}

approx3 = function(x){
  return(1/(2*sqrt(pi))*(H0(x) - H2(x)/4))
}

approx5 = function(x){
  return(1/(2*sqrt(pi))*(H0(x) - H2(x)/4 + H4(x)/32))
}

approx7 = function(x){
  return(1/(2*sqrt(pi))*(H0(x) - H2(x)/4 + H4(x)/32 - H6(x)/384))
}

approx9 = function(x){
  return(1/(2*sqrt(pi))*(H0(x) - H2(x)/4 + H4(x)/32 - H6(x)/384 + H8(x)/6144))
}

curve(dnorm(x), from=-4, to=4)
curve(approx1(x), add=TRUE, col="red")
curve(approx3(x), add=TRUE, col="blue")
curve(approx5(x), add=TRUE, col="green")
curve(approx7(x), add=TRUE, col="orange")
curve(approx9(x), add=TRUE, col="purple")

legend("topright", 
       legend=c("dnorm", "He0", "He2", "He4", "He6", "He8"),
       col=c("black", "red", "blue", "green", "orange", "purple"),
       lty=1)
