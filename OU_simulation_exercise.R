
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





