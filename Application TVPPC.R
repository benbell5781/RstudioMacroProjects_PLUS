####
## BEGINNING OF CODE
####

## DATA PROCESSING FOR VARIABLES:
#
#inflation
install.packages("eurostat")
library(eurostat)

hicp <- get_eurostat("prc_hicp_midx", 
                     filters = list(geo = "AT", 
                                    coicop = "CP00",
                                    unit = "I15"))
# Convert to year-on-year inflation
library(dplyr)
hicp <- hicp %>%
  arrange(time) %>%
  mutate(inflation = (values / lag(values, 12) - 1) * 100)

# Unemployment
unemp <- get_eurostat("une_rt_m",
                      filters = list(geo = "AT",
                                     s_adj = "SA",
                                     age = "TOTAL",
                                     sex = "T",
                                     unit = "PC_ACT"))
library(readxl)
#Output Gap
output_gap <- read_xlsx("AT_output_gap.xlsx")


library(dplyr)

head(output_gap)
head(hicp)
head(unemp)

# 1. Inflation
hicp <- get_eurostat("prc_hicp_midx",
                     filters = list(geo = "AT",
                                    coicop = "CP00", 
                                    unit = "I15")) %>%
  arrange(time) %>%
  mutate(inflation = (values / lag(values, 12) - 1) * 100) %>%
  mutate(quarter = as.yearqtr(time)) %>%
  group_by(quarter) %>%
  summarise(inflation = mean(inflation, na.rm = TRUE))

# 2. Unemployment
unemp <- get_eurostat("une_rt_m",
                      filters = list(geo = "AT",
                                     s_adj = "SA",
                                     age = "TOTAL",
                                     sex = "T",
                                     unit = "PC_ACT")) %>%
  mutate(quarter = as.yearqtr(time)) %>%
  group_by(quarter) %>%
  summarise(unemp = mean(values, na.rm = TRUE))


# Full table
library(dplyr)
library(zoo)

# --- Clean output gap ---
output_gap_clean <- output_gap %>%
  # Remove first two junk rows
  slice(-c(1,2)) %>%
  # Rename columns
  rename(quarter_str = `Time period`, gap = `Percent per annum, 2015`) %>%
  # Convert date string "1990-Q1" to yearqtr
  mutate(quarter = as.yearqtr(quarter_str, format = "%Y-Q%q"),
         gap = as.numeric(gap)) %>%
  select(quarter, gap)

# --- Merge all three ---
data <- hicp %>%
  left_join(unemp, by = "quarter") %>%
  left_join(output_gap_clean, by = "quarter") %>%
  # Filter to your window
  filter(quarter >= as.yearqtr("1997 Q1") & 
           quarter <= as.yearqtr("2008 Q4")) %>%
  # Make all value columns numeric
  mutate(across(c(inflation, unemp, gap), as.numeric)) %>%
  # Drop any rows with NaN or NA
  filter(if_all(c(inflation, unemp, gap), ~ !is.nan(.) & !is.na(.))) %>%
  arrange(quarter)

# Check result
head(data)
str(data)


# ---- packages -------------------------------------------------------------
suppressPackageStartupMessages({
  require(coda); require(GIGrvg); require(MASS); require(Matrix)
  require(shrinkTVP); require(stochvol); require(zoo); require(mvtnorm)
  require(ggplot2); require(reshape2); require(bayesm); require(scales)
})

dir.create("figures", showWarnings = FALSE)

# ---- data ----------------------------------------------------------------

data <- data %>%
  mutate(inflation_lag = lag(inflation, 1)) %>%
  filter(!is.na(inflation_lag))

# Dependent variable
y <- data$inflation

library(dplyr)
library(zoo)

# ---- Prepare your data ----
y <- as.numeric(data$inflation)
X <- cbind(1,
           as.numeric(data$inflation_lag),
           as.numeric(data$unemp),
           as.numeric(data$gap))

T <- length(y)      # number of time periods
K <- ncol(X)        # number of coefficients (4)

# ---- Priors ----
# Initial state beta_0
beta_0 <- rep(0, K)
P_0    <- diag(K) * 10      # diffuse prior

# Conjugate priors for all model parameters
# beta_0 given a diffuse normal prior centered at zero, with a large prior 
# variance (P_0 = 10*I), this means there is weak information about starting 
# values of TVPs (not imposing strong beliefs about beta_0 beginning)


# Q ~ Inverse-Wishart (innovation covariance)
nu_Q <- K + 1
Q_0  <- diag(K) * 0.01

# Innovation covariance matrix Q, this governs how quickly coefficients evolve
# Its given an inverse- wishart prior with K+1 DoF and a small prior scale
# matrix (Q_0 = 0.01*I), meaning coefficient evolution is gradual (don't want to
# fit noise)

# sigma^2 ~ Inverse-Gamma
a_0 <- 0.01
b_0 <- 0.01

# Observation error variance sigma^2 is given a weakly informative inverse gamma
# prior

# ---- MCMC settings ----
nsave <- 5000
nburn <- 2000
ntotal <- nsave + nburn

# Gibbs is run for 7,000 total iterations, first 2,000 discarded as burn in, 
# this leaves 5,000 retained draws for posterior inference.

# ---- Storage ----
beta_store  <- array(NA, dim = c(nsave, T, K))  # time varying coefficients
sigma_store <- rep(NA, nsave)                     # observation variance
Q_store     <- array(NA, dim = c(nsave, K, K))   # innovation covariance

# ---- Initialise ----
beta_draw  <- matrix(0, T, K)   # T x K matrix of coefficients
sigma2     <- var(y)
Q_draw     <- diag(K) * 0.01

# ================================================================
# Carter-Kohn Forward-Backward Sampler
# Draws the full path {beta_1,...,beta_T} in one block
# ================================================================

# Joint posterior of full coefficient path has no closed form Carter Kohn (1994)
# forward-filtering, backward-sampling algorithim is used draw entire sequence
# of TVPs in a single block at each Gibbs iteration, this happens in two passes

carter_kohn <- function(y, X, sigma2, Q, beta_0, P_0) {
  
  T_len <- length(y)
  k_dim <- length(beta_0)
  
  # --- Forward pass (Kalman filter) ---
  beta_filt <- matrix(NA, T_len, k_dim)
  P_filt    <- array(NA, c(k_dim, k_dim, T_len))
  
  beta_pred <- beta_0
  P_pred    <- P_0
  
  for (t in 1:T_len) {
    Xt <- matrix(X[t, ], nrow = 1)  # 1 x k_dim
    
    # Innovation variance
    v_t <- as.numeric(y[t] - Xt %*% beta_pred)
    S_t <- as.numeric(Xt %*% P_pred %*% t(Xt) + sigma2)
    
    # Kalman gain (k_dim x 1)
    Kgain <- as.vector(P_pred %*% t(Xt)) / S_t
    
    # Update
    beta_filt[t, ] <- beta_pred + Kgain * v_t
    P_filt[,,t]    <- P_pred - outer(Kgain, as.vector(Xt %*% P_pred))
    
    # Predict next period
    beta_pred <- beta_filt[t, ]
    P_pred    <- P_filt[,,t] + Q
  }
  
# Kalman Filter, algorithm steps through time from t=1 to T, producing at each 
# period a filtered estimate of beta_t conditional on data observed up to that
# point. 
# This uses Kalman Filter recursions, computing prediction error, the Kalman
# gain and updating both the mean and variance of the TVP.
  
  # --- Backward pass (simulation smoother) ---
  beta_draw <- matrix(NA, T_len, k_dim)
  
  # Draw at T
  beta_draw[T_len, ] <- MASS::mvrnorm(1, beta_filt[T_len, ], P_filt[,,T_len])
  
  # Draw backwards
  for (t in (T_len - 1):1) {
    P_t  <- P_filt[,,t]
    P_t1 <- P_filt[,,t] + Q
    
    # Smoothing gain
    J_t  <- P_t %*% solve(P_t1)
    
    # Smoothed mean and variance
    m_smooth <- beta_filt[t, ] + J_t %*% (beta_draw[t+1, ] - beta_filt[t, ])
    P_smooth <- P_t - J_t %*% P_t1 %*% t(J_t)
    
    # Ensure P_smooth is symmetric (numerical stability)
    P_smooth <- (P_smooth + t(P_smooth)) / 2
    
    beta_draw[t, ] <- MASS::mvrnorm(1, as.vector(m_smooth), P_smooth)
  }
  
  return(beta_draw)
}

# Simulation Smoother, then the algorithm steps back from t=T to t=1, drawing 
# each beta_t conditional on already drawn beta_t+1 and the filtered estimates 
# from the forward pass. 
# This produces a single draw of the entire coefficient path, correctly 
# accounting for the full joint posterior, rather than treating each period 
# independently

# ================================================================
# Gibbs Sampler
# ================================================================
cat("Starting Gibbs sampler...\n")

for (iter in 1:ntotal) {
  
  # --- Step 1: Draw full beta path via Carter-Kohn ---
  beta_draw <- carter_kohn(y, X, sigma2, Q_draw, beta_0, P_0)

# Conditional on current draws of sigma^2 and Q, seqeunce of TVPs is drawn using
# Carter Kohn algorithim previously described
  
  # --- Step 2: Draw sigma^2 | beta, Q, data ---
  # Compute residuals from observation equation
  resid <- y - rowSums(X * beta_draw)  # T x 1
  
  a_1 <- a_0 + T / 2
  b_1 <- b_0 + 0.5 * sum(resid^2)
  sigma2 <- 1 / rgamma(1, shape = a_1, rate = b_1)

# Conditional on newly drawn coefficient path, residuals from observation 
# equation are computed at each time period.
# Posterior for sigma^2 is inverse gamma with updated shape a_1 = a_0 + T/2
# and updated scale b_1 = b_0 + 1/2∑(residuals^2).
# A draw is obtained by sampling from a gamma distribution and inverting.
  
  # --- Step 3: Draw Q | beta, sigma^2, data ---
  # Compute innovations from state equation
  innov <- diff(beta_draw)  # (T-1) x K, each row is eta_t = beta_t - beta_{t-1}
  
  nu_1 <- nu_Q + (T - 1)
  Q_1  <- Q_0 + t(innov) %*% innov  # K x K
  
  Q_draw <- solve(rWishart(1, df = nu_1, Sigma = solve(Q_1))[,,1])

# Innovations ita_t = beta_t - beta_t-1 are computed by diffrencing the drawn
# coefficient path.
# The posterior for Q is inverse wishart with updated DoF v_1 = v_0 + (T - 1)
# and updated scale matrix Q_1 = Q_0 + ∑ita_t*ita_t'
    
  # --- Store post burn-in ---
  if (iter > nburn) {
    s <- iter - nburn
    beta_store[s, , ]  <- beta_draw
    sigma_store[s]     <- sigma2
    Q_store[s, , ]     <- Q_draw
  }
  
  if (iter %% 500 == 0) cat("Iteration", iter, "of", ntotal, "\n")
}

# Discards first 2000 observations are discarded in burn in, as early draws are 
# influenced by starting values (beta = 0, sigma^2 = var(y), Q = 0.01*I) rather
# than true posterior. 
# Remaining 5,000 draws are stored and treated as samples from the joint 
# posterior distribution of all TVPs.

cat("Gibbs sampler complete.\n")

# ================================================================
# Results — Posterior summaries of beta_t
# ================================================================
b_mean <- apply(beta_store, c(2, 3), mean)
b_lo95 <- apply(beta_store, c(2, 3), quantile, 0.025)
b_hi95 <- apply(beta_store, c(2, 3), quantile, 0.975)

# Time axis
time_index <- as.numeric(data$quarter)

# ================================================================
# Plot time varying coefficients
# ================================================================
library(ggplot2)

coef_names <- c("Intercept", "Lagged Inflation", "Unemployment", "Output Gap")
par(mfrow = c(2, 2), 
    oma = c(5, 0, 0, 0),   # outer bottom margin for legend
    mar = c(4, 4, 3, 1))   # increase top margin on each plot for title space

for (k in 1:K) {
  plot(data$quarter, b_mean[, k], type = "l", col = "navyblue", lwd = 2,
       main = paste("TVP:", coef_names[k]),
       xlab = "Quarter", ylab = "Coefficient",
       ylim = range(c(b_lo95[, k], b_hi95[, k])))
  
  polygon(c(data$quarter, rev(data$quarter)),
          c(b_lo95[, k], rev(b_hi95[, k])),
          col = rgb(0, 0, 1, 0.12), border = NA)
  
  abline(v = as.numeric(as.yearqtr("1999 Q1")),
         col = "red", lty = 2, lwd = 1.5)
}

# Draw single shared legend in the outer bottom margin
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 45, 0), mar = c(0, 0, 0, 0), new = TRUE)
legend("bottom",
       legend = c("Posterior mean", "95% credible band", "Euro adoption 1999 Q1"),
       col    = c("navyblue", rgb(0, 0, 1, 0.3), "red"),
       lty    = c(1, NA, 2),
       lwd    = c(2, NA, 1.5),
       pch    = c(NA, 15, NA),
       pt.cex = 2,
       horiz  = TRUE,      # horizontal legend
       bty    = "n",       # no box around legend
       xpd    = TRUE)
# ================================================================
# Specification 2 — Lagged Inflation + Unemployment only
# ================================================================

X2 <- cbind(1,
            as.numeric(data$inflation_lag),
            as.numeric(data$unemp))

K2 <- ncol(X2)

# Reinitialise priors for new dimension
beta_0_2 <- rep(0, K2)
P_0_2    <- diag(K2) * 10
nu_Q_2   <- K2 + 1
Q_0_2    <- diag(K2) * 0.01

# Storage
beta_store2  <- array(NA, dim = c(nsave, T, K2))
sigma_store2 <- rep(NA, nsave)
Q_store2     <- array(NA, dim = c(nsave, K2, K2))

# Initialise
beta_draw2 <- matrix(0, T, K2)
sigma2_2   <- var(y)
Q_draw2    <- diag(K2) * 0.01

cat("Starting Gibbs sampler - Specification 2...\n")

for (iter in 1:ntotal) {
  
  beta_draw2 <- carter_kohn(y, X2, sigma2_2, Q_draw2, beta_0_2, P_0_2)
  
  resid2 <- y - rowSums(X2 * beta_draw2)
  a_1_2  <- a_0 + T / 2
  b_1_2  <- b_0 + 0.5 * sum(resid2^2)
  sigma2_2 <- 1 / rgamma(1, shape = a_1_2, rate = b_1_2)
  
  innov2 <- diff(beta_draw2)
  nu_1_2 <- nu_Q_2 + (T - 1)
  Q_1_2  <- Q_0_2 + t(innov2) %*% innov2
  Q_draw2 <- solve(rWishart(1, df = nu_1_2, Sigma = solve(Q_1_2))[,,1])
  
  if (iter > nburn) {
    s <- iter - nburn
    beta_store2[s, , ]  <- beta_draw2
    sigma_store2[s]     <- sigma2_2
    Q_store2[s, , ]     <- Q_draw2
  }
  
  if (iter %% 500 == 0) cat("Iteration", iter, "of", ntotal, "\n")
}

cat("Specification 2 complete.\n")

# Posterior summaries
b_mean2 <- apply(beta_store2, c(2, 3), mean)
b_lo95_2 <- apply(beta_store2, c(2, 3), quantile, 0.025)
b_hi95_2 <- apply(beta_store2, c(2, 3), quantile, 0.975)

# Plot specification 2
coef_names2 <- c("Intercept", "Lagged Inflation", "Unemployment")

par(mfrow = c(1, 3),
    mar = c(4, 4, 3, 1),
    oma = c(5, 0, 2, 0))

for (k in 1:K2) {
  plot(data$quarter, b_mean2[, k], type = "l", col = "navyblue", lwd = 2,
       xlab = "Quarter", ylab = "Coefficient",
       ylim = range(c(b_lo95_2[, k], b_hi95_2[, k])))
  title(main = paste("TVP:", coef_names2[k]), line = 1.5, cex.main = 1)
  polygon(c(data$quarter, rev(data$quarter)),
          c(b_lo95_2[, k], rev(b_hi95_2[, k])),
          col = rgb(0, 0, 1, 0.15), border = NA)
  abline(v = as.numeric(as.yearqtr("1999 Q1")),
         col = "red", lty = 2, lwd = 1.5)
}

mtext("Specification 2: Lagged Inflation + Unemployment",
      outer = TRUE, cex = 1.2, font = 2, line = 0.5)



# ================================================================
# Specification 3 — Lagged Inflation + Output Gap only
# ================================================================

X3 <- cbind(1,
            as.numeric(data$inflation_lag),
            as.numeric(data$gap))


K3 <- ncol(X3)

# Reinitialise priors
beta_0_3 <- rep(0, K3)
P_0_3    <- diag(K3) * 10
nu_Q_3   <- K3 + 1
Q_0_3    <- diag(K3) * 0.01

# Storage
beta_store3  <- array(NA, dim = c(nsave, T, K3))
sigma_store3 <- rep(NA, nsave)
Q_store3     <- array(NA, dim = c(nsave, K3, K3))

# Initialise
beta_draw3 <- matrix(0, T, K3)
sigma2_3   <- var(y)
Q_draw3    <- diag(K3) * 0.01

cat("Starting Gibbs sampler - Specification 3...\n")

for (iter in 1:ntotal) {
  
  beta_draw3 <- carter_kohn(y, X3, sigma2_3, Q_draw3, beta_0_3, P_0_3)
  
  resid3 <- y - rowSums(X3 * beta_draw3)
  a_1_3  <- a_0 + T / 2
  b_1_3  <- b_0 + 0.5 * sum(resid3^2)
  sigma2_3 <- 1 / rgamma(1, shape = a_1_3, rate = b_1_3)
  
  innov3 <- diff(beta_draw3)
  nu_1_3 <- nu_Q_3 + (T - 1)
  Q_1_3  <- Q_0_3 + t(innov3) %*% innov3
  Q_draw3 <- solve(rWishart(1, df = nu_1_3, Sigma = solve(Q_1_3))[,,1])
  
  if (iter > nburn) {
    s <- iter - nburn
    beta_store3[s, , ]  <- beta_draw3
    sigma_store3[s]     <- sigma2_3
    Q_store3[s, , ]     <- Q_draw3
  }
  
  if (iter %% 500 == 0) cat("Iteration", iter, "of", ntotal, "\n")
}

cat("Specification 3 complete.\n")

# Posterior summaries
b_mean3  <- apply(beta_store3, c(2, 3), mean)
b_lo95_3 <- apply(beta_store3, c(2, 3), quantile, 0.025)
b_hi95_3 <- apply(beta_store3, c(2, 3), quantile, 0.975)

# Plot specification 3
coef_names3 <- c("Intercept", "Lagged Inflation", "Output Gap")

par(mfrow = c(1, 3),
    mar = c(4, 4, 3, 1),
    oma = c(5, 0, 2, 0))

for (k in 1:K3) {
  plot(data$quarter, b_mean3[, k], type = "l", col = "navyblue", lwd = 2,
       xlab = "Quarter", ylab = "Coefficient",
       ylim = range(c(b_lo95_3[, k], b_hi95_3[, k])))
  title(main = paste("TVP:", coef_names3[k]), line = 1.5, cex.main = 1)
  polygon(c(data$quarter, rev(data$quarter)),
          c(b_lo95_3[, k], rev(b_hi95_3[, k])),
          col = rgb(0, 0, 1, 0.15), border = NA)
  abline(v = as.numeric(as.yearqtr("1999 Q1")),
         col = "red", lty = 2, lwd = 1.5)
}

mtext("Specification 3: Lagged Inflation + Output Gap",
      outer = TRUE, cex = 1.2, font = 2, line = 0.5)

par(fig = c(0, 1, 0, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
legend("bottom",
       legend = c("Posterior mean", "95% credible band", "Euro adoption 1999 Q1"),
       col    = c("navyblue", rgb(0, 0, 1, 0.3), "red"),
       lty    = c(1, NA, 2),
       lwd    = c(2, NA, 1.5),
       pch    = c(NA, 15, NA),
       pt.cex = 2,
       horiz  = TRUE,
       bty    = "n",
       xpd    = TRUE)

# ================================================================
# Posterior comparison: 1999 Q1 vs final quarter
# Specification 2 (Lagged Inflation + Unemployment)
# Coefficient of interest: unemployment (column 3 in X2 -> beta index 3)
# ================================================================

# Find the index for 1999 Q1 and the final quarter in your sample
idx_1999q1 <- which(data$quarter == as.yearqtr("1999 Q1"))
idx_final  <- nrow(data)  # last row = final quarter in sample

cat("1999 Q1 is row:", idx_1999q1, "\n")
cat("Final quarter is row:", idx_final, "(", 
    as.character(data$quarter[idx_final]), ")\n")

# Extract the full posterior of beta_t for unemployment (3rd coefficient)
# beta_store2 has dimensions (nsave x T x K2)
unemp_coef_idx <- 3  # intercept=1, lag_inflation=2, unemployment=3

post_1999q1 <- beta_store2[, idx_1999q1, unemp_coef_idx]
post_final  <- beta_store2[, idx_final,  unemp_coef_idx]

# ---- Plot: histogram + kernel density overlay ----
library(scales)  # for alpha() transparency

# Compute densities
dens_1999q1 <- density(post_1999q1)
dens_final  <- density(post_final)

# Set up plot range to cover both distributions
xlim_range <- range(c(post_1999q1, post_final))
ylim_range <- range(c(dens_1999q1$y, dens_final$y))

par(mar = c(5, 4, 4, 2))

# Histogram for 1999 Q1 (semi-transparent)
hist(post_1999q1, breaks = 40, probability = TRUE,
     col = alpha("navyblue", 0.4), border = NA,
     xlim = xlim_range, ylim = ylim_range,
     main = "Posterior of Unemployment Coefficient:\n1999 Q1 vs Final Quarter",
     xlab = expression(beta[t]^{unemployment}),
     ylab = "Density")

# Histogram for final quarter (semi-transparent, overlaid)
hist(post_final, breaks = 40, probability = TRUE,
     col = alpha("darkorange", 0.4), border = NA,
     add = TRUE)

# Kernel density lines on top
lines(dens_1999q1, col = "navyblue", lwd = 2.5)
lines(dens_final,  col = "darkorange", lwd = 2.5)

# Reference line at zero
abline(v = 0, col = "grey40", lty = 3, lwd = 1.5)

# Legend
legend("topright",
       legend = c(paste0("1999 Q1"), 
                  paste0(as.character(data$quarter[idx_final])),
                  "Zero"),
       fill   = c(alpha("navyblue", 0.4), alpha("darkorange", 0.4), NA),
       border = c(NA, NA, NA),
       lty    = c(NA, NA, 3),
       lwd    = c(NA, NA, 1.5),
       col    = c(NA, NA, "grey40"),
       bty    = "n")

# ---- Print summary statistics ----
cat("\n--- Posterior Summary: Unemployment Coefficient ---\n")
cat("1999 Q1:\n")
cat("  Mean:  ", round(mean(post_1999q1), 4), "\n")
cat("  SD:    ", round(sd(post_1999q1), 4), "\n")
cat("  95% CI: [", round(quantile(post_1999q1, 0.025), 4), ",", 
    round(quantile(post_1999q1, 0.975), 4), "]\n\n")

cat("Final Quarter (", as.character(data$quarter[idx_final]), "):\n")
cat("  Mean:  ", round(mean(post_final), 4), "\n")
cat("  SD:    ", round(sd(post_final), 4), "\n")
cat("  95% CI: [", round(quantile(post_final, 0.025), 4), ",", 
    round(quantile(post_final, 0.975), 4), "]\n\n")

# Posterior probability that coefficient flattened (moved toward zero)
prob_flatten <- mean(abs(post_final) < abs(post_1999q1))
cat("P(|beta_final| < |beta_1999Q1| | data) =", round(prob_flatten, 4), "\n")

### Robustness Check ###
# ================================================================
# Robustness check: Posterior comparison across all 3 specifications
# Coefficient of interest: unemployment (where present) 
# Note: Specification 3 has no unemployment - compare output gap instead,
# or just compare unemployment in Spec 1 and Spec 2, and output gap in
# Spec 1 and Spec 3 separately
# ================================================================

idx_1999q1 <- which(data$quarter == as.yearqtr("1999 Q1"))
idx_final  <- nrow(data)

# ---- Unemployment coefficient: Specification 1 vs Specification 2 ----
# Spec 1 (full model): intercept=1, lag_inf=2, unemp=3, gap=4
# Spec 2 (lag+unemp):  intercept=1, lag_inf=2, unemp=3

post_1999q1_spec1_u <- beta_store[, idx_1999q1, 3]   # unemployment in full model
post_final_spec1_u  <- beta_store[, idx_final, 3]

post_1999q1_spec2_u <- beta_store2[, idx_1999q1, 3]  # unemployment in spec 2
post_final_spec2_u  <- beta_store2[, idx_final, 3]

# ---- Output gap coefficient: Specification 1 vs Specification 3 ----
# Spec 1 (full model): gap = column 4
# Spec 3 (lag+gap):    intercept=1, lag_inf=2, gap=3

post_1999q1_spec1_g <- beta_store[, idx_1999q1, 4]   # output gap in full model
post_final_spec1_g  <- beta_store[, idx_final, 4]

post_1999q1_spec3_g <- beta_store3[, idx_1999q1, 3]  # output gap in spec 3
post_final_spec3_g  <- beta_store3[, idx_final, 3]

# ================================================================
# 2x2 panel: rows = unemployment / output gap, cols = with vs without
# ================================================================

library(scales)

par(mfrow = c(2, 2), mar = c(4.5, 4, 3, 1), oma = c(5, 0, 2, 0))

plot_posterior_compare <- function(post_early, post_late, title_text, 
                                   label_early = "1999 Q1", 
                                   label_late  = "Final Quarter") {
  
  dens_early <- density(post_early)
  dens_late  <- density(post_late)
  
  xlim_range <- range(c(post_early, post_late))
  ylim_range <- range(c(dens_early$y, dens_late$y))
  
  hist(post_early, breaks = 40, probability = TRUE,
       col = alpha("navyblue", 0.4), border = NA,
       xlim = xlim_range, ylim = ylim_range,
       main = title_text, xlab = "Coefficient value", ylab = "Density")
  hist(post_late, breaks = 40, probability = TRUE,
       col = alpha("darkorange", 0.4), border = NA, add = TRUE)
  
  lines(dens_early, col = "navyblue", lwd = 2.5)
  lines(dens_late,  col = "darkorange", lwd = 2.5)
  abline(v = 0, col = "grey40", lty = 3, lwd = 1.5)
}

# Panel 1: Unemployment, Specification 1 (full model)
plot_posterior_compare(post_1999q1_spec1_u, post_final_spec1_u,
                       "Unemployment Coef.\n(Full Model - Spec 1)")

# Panel 2: Unemployment, Specification 2 (robustness)
plot_posterior_compare(post_1999q1_spec2_u, post_final_spec2_u,
                       "Unemployment Coef.\n(Lag + Unemp - Spec 2)")

# Panel 3: Output Gap, Specification 1 (full model)
plot_posterior_compare(post_1999q1_spec1_g, post_final_spec1_g,
                       "Output Gap Coef.\n(Full Model - Spec 1)")

# Panel 4: Output Gap, Specification 3 (robustness)
plot_posterior_compare(post_1999q1_spec3_g, post_final_spec3_g,
                       "Output Gap Coef.\n(Lag + Gap - Spec 3)")

mtext("Robustness Check: Posterior Comparison Across Specifications",
      outer = TRUE, cex = 1.2, font = 2, line = 0.5)

par(fig = c(0, 1, 0, 1), oma = c(0,9,0,0), mar = c(0,0,0,0), new = TRUE)
legend("bottom",
       legend = c("1999 Q1", "Final Quarter", "Zero"),
       fill   = c(alpha("navyblue", 0.4), alpha("darkorange", 0.4), NA),
       border = c(NA, NA, NA),
       lty    = c(NA, NA, 3),
       lwd    = c(NA, NA, 1.5),
       col    = c(NA, NA, "grey40"),
       horiz  = TRUE,
       bty    = "n",
       xpd    = TRUE)

# Spec 1 — Unemployment coefficient
prob_flatten_spec1_u <- mean(abs(post_final_spec1_u) < abs(post_1999q1_spec1_u))
cat("Spec 1 - Unemployment P(flattened | data) =", round(prob_flatten_spec1_u, 4), "\n")

# Spec 1 — Output Gap coefficient
prob_flatten_spec1_g <- mean(abs(post_final_spec1_g) < abs(post_1999q1_spec1_g))
cat("Spec 1 - Output Gap P(flattened | data) =", round(prob_flatten_spec1_g, 4), "\n")

# For specification 2 (unemployment)
prob_flatten_u <- mean(abs(post_final_spec2_u) < abs(post_1999q1_spec2_u))
cat("P(flattened | data) =", round(prob_flatten_u, 4), "\n")

# For specification 3 (output gap)
prob_flatten_g <- mean(abs(post_final_spec3_g) < abs(post_1999q1_spec3_g))
cat("P(flattened | data) =", round(prob_flatten_g, 4), "\n")

# ================================================================
# Summary table across specifications
# ================================================================

cat("\n=== ROBUSTNESS CHECK SUMMARY ===\n\n")

cat("--- Unemployment Coefficient ---\n")
cat("Spec 1 (Full Model):\n")
cat("  1999 Q1 mean:", round(mean(post_1999q1_spec1_u), 4),
    " | Final mean:", round(mean(post_final_spec1_u), 4), "\n")
cat("  P(flattened) =", round(mean(abs(post_final_spec1_u) < abs(post_1999q1_spec1_u)), 4), "\n\n")


cat("Spec 2 (Lag + Unemployment):\n")
cat("  1999 Q1 mean:", round(mean(post_1999q1_spec2_u), 4),
    " | Final mean:", round(mean(post_final_spec2_u), 4), "\n")
cat("  P(flattened) =", round(mean(abs(post_final_spec2_u) < abs(post_1999q1_spec2_u)), 4), "\n\n")

cat("--- Output Gap Coefficient ---\n")
cat("Spec 1 (Full Model):\n")
cat("  1999 Q1 mean:", round(mean(post_1999q1_spec1_g), 4),
    " | Final mean:", round(mean(post_final_spec1_g), 4), "\n")
cat("  P(flattened) =", round(mean(abs(post_final_spec1_g) < abs(post_1999q1_spec1_g)), 4), "\n\n")

cat("Spec 3 (Lag + Output Gap):\n")
cat("  1999 Q1 mean:", round(mean(post_1999q1_spec3_g), 4),
    " | Final mean:", round(mean(post_final_spec3_g), 4), "\n")
cat("  P(flattened) =", round(mean(abs(post_final_spec3_g) < abs(post_1999q1_spec3_g)), 4), "\n\n")

####
##END OF CODE
####
