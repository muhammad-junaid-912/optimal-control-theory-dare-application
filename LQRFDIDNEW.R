# Define matrices
A <- matrix(c(
  1, -28.83, 0, 0, -0.2,  
  0, 1, -0.07, 0, -0.01,
  0, -0.003, 1, -0.1, 0.0002,
  0, 0, 0.7, 1, -0.0004,
  0, 0, 0, 0, 1
), nrow=5, byrow=TRUE)

eig <- eigen(A)
cat("Eigenvalues of matrix A:\n")
print(eig$values)

B <- matrix(c(0, 0, 0, 0, 1), nrow = 5, ncol = 1)
Q <- diag(c(1,1,1,1,1)) 
R <- matrix(20, nrow = 1, ncol = 1)  # Scalar as a matrix

# Solve Discrete Algebraic Riccati Equation (DARE)
P <- Q  # Initial guess for P
tolerance <- 1e-6  # Convergence tolerance
max_iter <- 5000  # Maximum iterations

for (i in 1:max_iter) {
  P_prev <- P
  K <- solve(R + t(B) %*% P %*% B) %*% (t(B) %*% P %*% A)
  P <- Q + t(A) %*% P %*% A - t(A) %*% P %*% B %*% K
  
  if (max(abs(P - P_prev)) < tolerance) {
    cat("Converged after", i, "iterations.\n")
    break
  }
}

if (i == max_iter) {
  cat("Warning: Maximum iterations reached without convergence.\n")
}

# Optimal Feedback Gain
K <- solve(R) %*% t(B) %*% P

cat("Matrix P (Solution to DARE):\n")
print(P)
cat("\nOptimal Feedback Gain K:\n")
print(K)

# Define simulation parameters
time_steps <- 5 # Number of simulation steps
forecast_steps <- 5  # Number of forecast steps
x <- matrix(0, nrow = 5, ncol = time_steps + 1)  # State matrix
u <- numeric(time_steps + forecast_steps)  # Control input vector

# Initial state
x[, 1] <- c(71246, 281.92, 13.36, 26.83, 28536.8)  # Example initial condition

# Simulate the system (Control Active)
for (t in 1:time_steps) {
  u[t] <- -K %*% x[, t]
  x[, t + 1] <- A %*% x[, t] + B %*% u[t]
}
u <- numeric(time_steps + forecast_steps)  # Control input vector

# Inside your loop:
u[t] <- -K %*% x[, t]


# Forecast with Control (Now Applied)
forecast_x <- matrix(0, nrow = 5, ncol = forecast_steps + 1)  
forecast_x[, 1] <- x[, time_steps + 1]  

for (t in 1:forecast_steps) {
  u[time_steps + t] <- -K %*% forecast_x[, t]  # Control remains active
  forecast_x[, t + 1] <- A %*% forecast_x[, t] + B %*% u[time_steps + t]
}

# Combine results
combined_x <- cbind(x, forecast_x[, -1])

# Print forecasted values
cat("\nForecasted Values (5 steps ahead with control):\n")
print(forecast_x)
print("Optimal control inputs for next 5 years:")
print(u)
# Plot each state trajectory separately
par(mfrow = c(3, 2))  # Arrange plots in a grid (3 rows, 2 columns)
for (i in 1:5) {
  plot(0:(time_steps + forecast_steps), combined_x[i, ], type = "l", col = "blue", lwd = 2,
       xlab = "Time Steps", ylab = paste("State x", i),
       main = paste("State x", i, "Trajectory"))
  abline(v = time_steps, col = "red", lty = 2)  # Mark transition point
}
par(mfrow = c(1, 1))  # Reset plotting layout


# Combine x and forecast_x to get full trajectory
full_x <- cbind(x[, 1:(time_steps + 1)], forecast_x[, -1])
full_u <- u[1:(time_steps + forecast_steps)]

# Compute cost
J <- 0
for (t in 1:(time_steps + forecast_steps)) {
  xt <- matrix(full_x[, t], ncol = 1)
  ut <- matrix(full_u[t], ncol = 1)
  J <- J + t(xt) %*% Q %*% xt + t(ut) %*% R %*% ut
}
cat("Optimal cost J over 10 steps:", J, "\n")
