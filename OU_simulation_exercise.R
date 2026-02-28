
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



# Make a plot of OU process of Hermite expansion
# compare to the analytic solution
beta = 1
sigma = 1
delta_t = 0.5
y0 = 0

# mu_Y(y) = -beta*y
# Y_{t+Δ} | Y_t=y0 ~ N(y0*exp(-beta*Δ), (1-exp(-2*beta*Δ))/(2*beta))
# Z = Δ^{-1/2}(Y_{t+Δ} - y0)

mean_y = y0 * exp(-beta * delta_t)
var_y = (1 - exp(-2 * beta * delta_t)) / (2 * beta)

mean_z = delta_t^(-1/2) * (mean_y - y0)
var_z = delta_t^(-1) * var_y
sd_z = sqrt(var_z)


H0 = function(x) rep(1, length(x))
H1 = function(x) x
H2 = function(x) x^2 - 1
H3 = function(x) x^3 - 3*x
H4 = function(x) x^4 - 6*x^2 + 3
H5 = function(x) x^5 - 10*x^3 + 15*x
H6 = function(x) x^6 - 15*x^4 + 45*x^2 - 15


eta0_3 = 1
eta1_3 = beta*y0*delta_t^(1/2) - (1/2)*beta^2*y0*delta_t^(3/2) + (1/6)*beta^3*y0*delta_t^(5/2)  # 修正：第一项是正号
eta2_3 = ((beta*y0^2-1)*beta*delta_t)/2 + ((-3*beta*y0^2+2)*beta^2*delta_t^2)/6 + ((7*beta*y0^2-4)*beta^3*delta_t^3)/24
eta3_3 = -((-beta*y0^2+3)*beta^2*y0*delta_t^(3/2))/6 - ((3*beta*y0^2-7)*beta^3*y0*delta_t^(5/2))/12
eta4_3 = ((beta^2*y0^4-6*beta*y0^2+3)*beta^2*delta_t^2)/24 + ((-beta^2*y0^4+5*beta*y0^2-2)*beta^3*delta_t^3)/12
eta5_3 = -((-beta^2*y0^4+10*beta*y0^2-15)*beta^3*y0*delta_t^(5/2))/120
eta6_3 = ((beta^3*y0^6-15*beta^2*y0^4-15+45*beta*y0^2)*beta^3*delta_t^3)/720


phi = function(x) exp(-x^2/2) / sqrt(2*pi)

approx0 = function(x) phi(x) * eta0_3 * H0(x)
approx1 = function(x) phi(x) * (eta0_3*H0(x) + eta1_3*H1(x))
approx2 = function(x) phi(x) * (eta0_3*H0(x) + eta1_3*H1(x) + eta2_3*H2(x))
approx3 = function(x) phi(x) * (eta0_3*H0(x) + eta1_3*H1(x) + eta2_3*H2(x) + eta3_3*H3(x))
approx4 = function(x) phi(x) * (eta0_3*H0(x) + eta1_3*H1(x) + eta2_3*H2(x) + eta3_3*H3(x) + eta4_3*H4(x))
approx5 = function(x) phi(x) * (eta0_3*H0(x) + eta1_3*H1(x) + eta2_3*H2(x) + eta3_3*H3(x) + eta4_3*H4(x) + eta5_3*H5(x))
approx6 = function(x) phi(x) * (eta0_3*H0(x) + eta1_3*H1(x) + eta2_3*H2(x) + eta3_3*H3(x) + eta4_3*H4(x) + eta5_3*H5(x) + eta6_3*H6(x))



true_density_z = function(z) dnorm(z, mean = mean_z, sd = sd_z)

curve(true_density_z(x), from=-3, to=3, col="black", lwd=2.5, 
      ylab="Density", xlab="z", 
      main=paste0("Hermite Expansion for Z (y0=", y0, ", Δ=", delta_t, ")"))
curve(approx1(x), add=TRUE, col="blue", lwd=1.5, lty=2)
curve(approx2(x), add=TRUE, col="green", lwd=1.5, lty=2)
curve(approx3(x), add=TRUE, col="orange", lwd=1.5, lty=2)
curve(approx6(x), add=TRUE, col="red", lwd=1.5, lty=2)

legend("topright", 
       legend=c("True Z density", "J=1", "J=2", "J=3", "J=6"),
       col=c("black", "blue", "green", "orange", "red"),
       lwd=c(2.5, 1.5, 1.5, 1.5, 1.5),
       lty=c(1, 2, 2, 2, 2))




# time varying coefficients
beta = 1
sigma = 1
delta_t = 0.5
y0 = 0
n_steps = 100
n_paths = 50
J = 3

# Hermite polynomials
H0 = function(x) rep(1, length(x))
H1 = function(x) x
H2 = function(x) x^2 - 1
H3 = function(x) x^3 - 3*x
H4 = function(x) x^4 - 6*x^2 + 3
H5 = function(x) x^5 - 10*x^3 + 15*x
H6 = function(x) x^6 - 15*x^4 + 45*x^2 - 15

# compute eta
compute_eta = function(y_current, delta_t, beta, sigma) {
  eta = numeric(7)
  
  #eta, J=0, K=3
  eta[1] = 1
  
  # eta, J=1, K=3
  eta[2] = beta*y_current*delta_t^(1/2) - 
    (1/2)*beta^2*y_current*delta_t^(3/2) + 
    (1/6)*beta^3*y_current*delta_t^(5/2)
  
  # eta, J=2, K=3
  eta[3] = ((beta*y_current^2-1)*beta*delta_t)/2 + 
    ((-3*beta*y_current^2+2)*beta^2*delta_t^2)/6 + 
    ((7*beta*y_current^2-4)*beta^3*delta_t^3)/24
  
  # eta, J=3, K=3
  eta[4] = -((-beta*y_current^2+3)*beta^2*y_current*delta_t^(3/2))/6 - 
    ((3*beta*y_current^2-7)*beta^3*y_current*delta_t^(5/2))/12
  
  # eta, J=4, K=3
  eta[5] = ((beta^2*y_current^4-6*beta*y_current^2+3)*beta^2*delta_t^2)/24 + 
    ((-beta^2*y_current^4+5*beta*y_current^2-2)*beta^3*delta_t^3)/12
  
  # eta, J=5, K=3
  eta[6] = -((-beta^2*y_current^4+10*beta*y_current^2-15)*beta^3*y_current*delta_t^(5/2))/120
  
  # eta, J=6, K=3
  eta[7] = ((beta^3*y_current^6-15*beta^2*y_current^4-15+45*beta*y_current^2)*beta^3*delta_t^3)/720
  
  return(eta)
}



# hermite_density = function(x_current, y0, delta_t, beta, sigma, J){
#   y = x_current / sigma
#   z = (y - y0)/sqrt(delta_t)
#   
#   phi_z = exp(-z^2/2) / sqrt(2*pi)
#   H_values = c(H0(z), H1(z), H2(z), H3(z), H4(z), H5(z), H6(z))
#   eta = compute_eta(y_current=y0, delta_t=delta_t, beta=beta, sigma=sigma)
#   sum_val = sum(eta[1:(J+1)] * H_values[1:(J+1)])
#   
#   p_z = phi_z * sum_val
#   p_y = p_z / sqrt(delta_t)
#   p_x = p_y / sigma
#   return(p_x)
# }


hermite_density = function(x_next, x_prev, delta_t, beta, sigma, J){
  y_next = x_next / sigma
  y_prev = x_prev / sigma  
  
  z = (y_next - y_prev)/sqrt(delta_t)
  
  phi_z = exp(-z^2/2) / sqrt(2*pi)
  H_values = c(H0(z), H1(z), H2(z), H3(z), H4(z), H5(z), H6(z))
  
  eta = compute_eta(y_current = y_prev, delta_t = delta_t, beta = beta, sigma = sigma)
  
  sum_val = sum(eta[1:(J+1)] * H_values[1:(J+1)])
  
  p_z = phi_z * sum_val
  p_y = p_z / sqrt(delta_t)
  p_x = p_y / sigma
  return(p_x)
}



# test
x_next = 0.5

p_approx = hermite_density(x_next, y0, delta_t, beta, sigma, J)

# exact density
mean_exact = y0 * exp(-beta * delta_t)
var_exact = sigma^2 / (2*beta) * (1 - exp(-2*beta*delta_t))
p_exact = dnorm(x, mean = mean_exact, sd = sqrt(var_exact))

cat("approx density:", p_approx, "\n")
cat("exact density:", p_exact, "\n")
cat("relative error:", abs(p_approx - p_exact) / p_exact * 100, "%\n")



# exact path (for likelihood)
simulate_path_exact = function(y0, beta, sigma, delta_t, n_steps) {
  y_path = numeric(n_steps + 1)
  y_path[1] = y0
  
  for(t in 2:(n_steps + 1)) {
    mean_next = y_path[t-1] * exp(-beta * delta_t)
    var_next = sigma^2 / (2*beta) * (1 - exp(-2*beta*delta_t))
    y_path[t] = rnorm(1, mean = mean_next, sd = sqrt(var_next))
  }
  
  return(y_path)
}