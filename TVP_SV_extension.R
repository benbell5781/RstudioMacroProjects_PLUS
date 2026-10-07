# ================================================================
# TVP-SV Phillips Curve — Stochastic Volatility Extension
# Extends the baseline TVP model by allowing sigma_t^2 to evolve
# over time via a log-random walk, following Kim, Shepherd & Chib (1998)
#
# Model:
#   Observation: pi_t = x_t' beta_t + eps_t,  eps_t ~ N(0, exp(h_t))
#   State:       beta_t = beta_{t-1} + eta_t,  eta_t ~ N(0, Q)
#   Volatility:  h_t = h_{t-1} + u_t,          u_t  ~ N(0, phi^2)
#
# Additional parameter vs baseline:
#   h_t = log(sigma_t^2) — log-volatility path (T x 1)
#   phi^2 — variance of log-volatility innovations
#
# Gibbs steps:
#   1. Draw {beta_t} | h, Q, data       via Carter-Kohn (same as before)
#   2. Draw Q        | {beta_t}         via Inverse-Wishart (same as before)
#   3. Draw {h_t}    | {beta_t}, phi^2  via Kim-Shepherd-Chib mixture approx
#   4. Draw phi^2    | {h_t}            via Inverse-Gamma
# ================================================================

# ----------------------------------------------------------------
# NOTE: Run AFTER TVPPC section (i.e. Below END OF CODE) code so that:
#   - y, X, X2, X3 are already defined
#   - T, K, K2, K3 are defined
#   - carter_kohn() function is defined
#   - idx_1999q1, idx_final are defined
#   - data$quarter is defined
# ----------------------------------------------------------------



cat("=== TVP-SV Model ===\n\n")

# ================================================================
# Kim-Shepherd-Chib (1998) mixture approximation constants
# Approximates log(chi^2(1)) as a 7-component Normal mixture
# ================================================================
KSC_p  <- c(0.0073, 0.1056, 0.6000, 0.2575, 0.0294, 0.0056, 0.0000)  # mixture weights
KSC_m  <- c(-10.12999, -3.97281, -8.56686,  2.77786,  0.61942,  1.79518, -1.08819) # means
KSC_v2 <- c(5.79596,   2.61369,  5.17950,  0.16735,  0.64009,  0.34023,  1.26261)  # variances

# ================================================================
# Helper: draw log-volatility path via KSC mixture approximation
# This is a single Carter-Kohn call on the transformed data
# ================================================================
draw_log_vol <- function(resid2_log, phi2) {
  # resid2_log = log(eps_t^2 + 0.001) — the "data" for the vol equation
  # This is a univariate state space model:
  #   resid2_log_t = h_t + xi_t,  xi_t ~ mixture Normal (KSC approx)
  #   h_t = h_{t-1} + u_t,        u_t  ~ N(0, phi^2)
  
  T_v <- length(resid2_log)
  h   <- rep(0, T_v)
  
  # --- Sample mixture indicators s_t ---
  s <- rep(1, T_v)
  for (t in 1:T_v) {
    # Compute posterior weight for each mixture component
    log_prob <- log(KSC_p) - 0.5 * log(KSC_v2) -
      0.5 * (resid2_log[t] - h[t] - KSC_m)^2 / KSC_v2
    log_prob <- log_prob - max(log_prob)
    prob     <- exp(log_prob) / sum(exp(log_prob))
    s[t]     <- sample(1:7, 1, prob = prob)
  }
  
  # --- Given indicators, run Carter-Kohn on univariate h state space ---
  # Observation: resid2_log_t - KSC_m[s_t] = h_t + noise ~ N(0, KSC_v2[s_t])
  y_tilde  <- resid2_log - KSC_m[s]
  obs_var  <- KSC_v2[s]   # time-varying observation variance
  
  # Forward pass
  h_filt <- rep(NA, T_v)
  P_filt <- rep(NA, T_v)
  h_pred <- 0
  P_pred <- 10   # diffuse initialisation
  
  for (t in 1:T_v) {
    v_t    <- y_tilde[t] - h_pred
    S_t    <- P_pred + obs_var[t]
    Kg     <- P_pred / S_t
    h_filt[t] <- h_pred + Kg * v_t
    P_filt[t] <- P_pred - Kg * P_pred
    h_pred <- h_filt[t]
    P_pred <- P_filt[t] + phi2
  }
  
  # Backward pass
  h_draw <- rep(NA, T_v)
  h_draw[T_v] <- rnorm(1, h_filt[T_v], sqrt(P_filt[T_v]))
  
  for (t in (T_v - 1):1) {
    P_t   <- P_filt[t]
    P_t1  <- P_filt[t] + phi2
    J_t   <- P_t / P_t1
    m_s   <- h_filt[t] + J_t * (h_draw[t+1] - h_filt[t])
    V_s   <- P_t - J_t^2 * P_t1
    h_draw[t] <- rnorm(1, m_s, sqrt(max(V_s, 1e-10)))
  }
  
  return(h_draw)
}

# ================================================================
# TVP-SV Gibbs Sampler — Specification 2 (Lag + Unemployment)
# Used as baseline for direct comparison with TVP baseline
# ================================================================

# Priors
beta_0_sv <- rep(0, K2)
P_0_sv    <- diag(K2) * 10
nu_Q_sv   <- K2 + 1
Q_0_sv    <- diag(K2) * 0.01

# phi^2 prior: Inverse-Gamma(a_phi, b_phi)
a_phi <- 0.01
b_phi <- 0.01

# MCMC settings (same as baseline)
nsave_sv <- 5000
nburn_sv <- 2000
ntotal_sv <- nsave_sv + nburn_sv

# Storage
beta_store_sv  <- array(NA, dim = c(nsave_sv, T, K2))
h_store_sv     <- matrix(NA, nrow = nsave_sv, ncol = T)  # log-vol path
phi2_store_sv  <- rep(NA, nsave_sv)
Q_store_sv     <- array(NA, dim = c(nsave_sv, K2, K2))

# Initialise
beta_draw_sv <- matrix(0, T, K2)
h_draw_sv    <- rep(0, T)          # log-volatility initialised at 0
phi2_sv      <- 0.1               # initial phi^2
Q_draw_sv    <- diag(K2) * 0.01

cat("Starting TVP-SV Gibbs sampler (Spec 2)...\n")

for (iter in 1:ntotal_sv) {
  
  # --- Step 1: Draw beta path | h, Q via Carter-Kohn ---
  # sigma_t^2 = exp(h_t) — pass as vector, Carter-Kohn handles scalar sigma
  # We use the mean exp(h) as a scalar approximation here for simplicity
  # (full SV Carter-Kohn uses time-varying obs variance per period)
  
  # Time-varying observation variance
  sigma_t2 <- exp(h_draw_sv)
  
  # Run Carter-Kohn with time-varying sigma at each t
  # Modified to handle vector sigma_t^2
  T_len <- T
  k_dim <- K2
  
  beta_filt_sv <- matrix(NA, T_len, k_dim)
  P_filt_sv    <- array(NA, c(k_dim, k_dim, T_len))
  beta_pred_sv <- beta_0_sv
  P_pred_sv    <- P_0_sv
  
  for (t in 1:T_len) {
    Xt <- matrix(X2[t, ], nrow = 1)
    v_t_sv <- as.numeric(y[t] - Xt %*% beta_pred_sv)
    S_t_sv <- as.numeric(Xt %*% P_pred_sv %*% t(Xt) + sigma_t2[t])  # time-varying!
    Kg_sv  <- as.vector(P_pred_sv %*% t(Xt)) / S_t_sv
    beta_filt_sv[t, ] <- beta_pred_sv + Kg_sv * v_t_sv
    P_filt_sv[,,t]    <- P_pred_sv - outer(Kg_sv, as.vector(Xt %*% P_pred_sv))
    beta_pred_sv <- beta_filt_sv[t, ]
    P_pred_sv    <- P_filt_sv[,,t] + Q_draw_sv
  }
  
  # Backward pass
  beta_draw_sv_new <- matrix(NA, T_len, k_dim)
  beta_draw_sv_new[T_len, ] <- MASS::mvrnorm(1, beta_filt_sv[T_len,], P_filt_sv[,,T_len])
  
  for (t in (T_len - 1):1) {
    P_t   <- P_filt_sv[,,t]
    P_t1  <- P_filt_sv[,,t] + Q_draw_sv
    J_t   <- P_t %*% solve(P_t1)
    m_s   <- beta_filt_sv[t,] + J_t %*% (beta_draw_sv_new[t+1,] - beta_filt_sv[t,])
    P_s   <- P_t - J_t %*% P_t1 %*% t(J_t)
    P_s   <- (P_s + t(P_s)) / 2
    beta_draw_sv_new[t,] <- MASS::mvrnorm(1, as.vector(m_s), P_s)
  }
  
  beta_draw_sv <- beta_draw_sv_new
  
  # --- Step 2: Draw Q | beta ---
  innov_sv <- diff(beta_draw_sv)
  nu_1_sv  <- nu_Q_sv + (T - 1)
  Q_1_sv   <- Q_0_sv + t(innov_sv) %*% innov_sv
  Q_draw_sv <- solve(rWishart(1, df = nu_1_sv, Sigma = solve(Q_1_sv))[,,1])
  
  # --- Step 3: Draw log-volatility path {h_t} | beta, phi^2 ---
  resid_sv     <- y - rowSums(X2 * beta_draw_sv)
  resid2_log   <- log(resid_sv^2 + 0.001)   # offset avoids log(0)
  h_draw_sv    <- draw_log_vol(resid2_log, phi2_sv)
  
  # --- Step 4: Draw phi^2 | {h_t} ---
  h_innov  <- diff(h_draw_sv)
  a_phi_1  <- a_phi + (T - 1) / 2
  b_phi_1  <- b_phi + 0.5 * sum(h_innov^2)
  phi2_sv  <- 1 / rgamma(1, shape = a_phi_1, rate = b_phi_1)
  
  # --- Store ---
  if (iter > nburn_sv) {
    s <- iter - nburn_sv
    beta_store_sv[s, , ] <- beta_draw_sv
    h_store_sv[s, ]      <- h_draw_sv
    phi2_store_sv[s]     <- phi2_sv
    Q_store_sv[s, , ]    <- Q_draw_sv
  }
  
  if (iter %% 500 == 0) cat("Iteration", iter, "of", ntotal_sv, "\n")
}

cat("TVP-SV complete.\n\n")

# ================================================================
# Posterior summaries — TVP-SV
# ================================================================
b_mean_sv  <- apply(beta_store_sv, c(2, 3), mean)
b_lo95_sv  <- apply(beta_store_sv, c(2, 3), quantile, 0.025)
b_hi95_sv  <- apply(beta_store_sv, c(2, 3), quantile, 0.975)

h_mean_sv  <- apply(h_store_sv, 2, mean)
h_lo95_sv  <- apply(h_store_sv, 2, quantile, 0.025)
h_hi95_sv  <- apply(h_store_sv, 2, quantile, 0.975)

# ================================================================
# Plot 1: Side-by-side coefficient comparison
# Baseline TVP vs TVP-SV for Specification 2
# ================================================================
coef_names2 <- c("Intercept", "Lagged Inflation", "Unemployment")

par(mfrow = c(3, 2),
    mar = c(3, 4, 3, 1),
    oma = c(5, 0, 3, 0))

for (k in 1:K2) {
  
  # --- Baseline TVP (left column) ---
  ylim_k <- range(c(b_lo95_2[, k], b_hi95_2[, k],
                    b_lo95_sv[, k], b_hi95_sv[, k]))
  
  plot(data$quarter, b_mean2[, k], type = "l",
       col = "navyblue", lwd = 2,
       ylim = ylim_k,
       xlab = "", ylab = "Coefficient")
  title(main = paste("Baseline TVP:", coef_names2[k]), line = 1)
  polygon(c(data$quarter, rev(data$quarter)),
          c(b_lo95_2[, k], rev(b_hi95_2[, k])),
          col = rgb(0, 0, 1, 0.15), border = NA)
  abline(v = as.numeric(as.yearqtr("1999 Q1")),
         col = "red", lty = 2, lwd = 1.5)
  abline(h = 0, col = "grey70", lty = 3)
  
  # --- TVP-SV (right column) ---
  plot(data$quarter, b_mean_sv[, k], type = "l",
       col = "darkgreen", lwd = 2,
       ylim = ylim_k,
       xlab = "", ylab = "Coefficient")
  title(main = paste("TVP-SV:", coef_names2[k]), line = 1)
  polygon(c(data$quarter, rev(data$quarter)),
          c(b_lo95_sv[, k], rev(b_hi95_sv[, k])),
          col = rgb(0, 0.5, 0, 0.15), border = NA)
  abline(v = as.numeric(as.yearqtr("1999 Q1")),
         col = "red", lty = 2, lwd = 1.5)
  abline(h = 0, col = "grey70", lty = 3)
}

mtext("Coefficient Paths: Baseline TVP (left) vs TVP-SV (right) — Spec 2",
      outer = TRUE, cex = 1.1, font = 2, line = 1)

par(fig = c(0, 1, 0, 1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
legend("bottom",
       legend = c("Baseline TVP mean", "Baseline 95% CI",
                  "TVP-SV mean",       "TVP-SV 95% CI",
                  "Euro 1999 Q1",      "Zero line"),
       col    = c("navyblue", rgb(0,0,1,0.3),
                  "darkgreen", rgb(0,0.5,0,0.3),
                  "red", "grey70"),
       lty    = c(1, NA, 1, NA, 2, 3),
       lwd    = c(2, NA, 2, NA, 1.5, 1),
       pch    = c(NA, 15, NA, 15, NA, NA),
       pt.cex = 1.5,
       horiz  = TRUE, bty = "n", xpd = TRUE, cex = 0.85)

# ================================================================
# Plot 2: Stochastic volatility path — what sigma_t^2 is doing
# ================================================================
par(mfrow = c(1, 1), mar = c(4, 4, 3, 2))

plot(data$quarter, exp(h_mean_sv), type = "l",
     col = "darkred", lwd = 2,
     main = "Estimated Time-Varying Observation Variance (TVP-SV)",
     xlab = "Quarter", ylab = expression(hat(sigma)[t]^2))
polygon(c(data$quarter, rev(data$quarter)),
        c(exp(h_lo95_sv), rev(exp(h_hi95_sv))),
        col = rgb(1, 0, 0, 0.15), border = NA)
abline(v = as.numeric(as.yearqtr("1999 Q1")),
       col = "navy", lty = 2, lwd = 1.5)
legend("topright",
       legend = c("Posterior mean sigma_t^2", "95% CI", "Euro 1999 Q1"),
       col    = c("darkred", rgb(1,0,0,0.3), "navy"),
       lty    = c(1, NA, 2), lwd = c(2, NA, 1.5),
       pch    = c(NA, 15, NA), pt.cex = 1.5, bty = "n")

# ================================================================
# Plot 3: Kernel density comparison — Baseline vs TVP-SV
# Unemployment coefficient: 1999 Q1 vs Final Quarter
# ================================================================
post_1999q1_sv <- beta_store_sv[, idx_1999q1, 3]
post_final_sv  <- beta_store_sv[, idx_final,  3]

par(mfrow = c(1, 2), mar = c(5, 4, 4, 2), oma = c(0, 0, 3, 0))

# --- Baseline ---
xlim_all <- range(c(post_1999q1_spec2_u, post_final_spec2_u,
                    post_1999q1_sv,       post_final_sv))

dens_b_99  <- density(post_1999q1_spec2_u)
dens_b_fin <- density(post_final_spec2_u)
ylim_b <- range(c(dens_b_99$y, dens_b_fin$y))

hist(post_1999q1_spec2_u, breaks = 40, probability = TRUE,
     col = alpha("navyblue", 0.4), border = NA,
     xlim = xlim_all, ylim = ylim_b,
     main = "Baseline TVP",
     xlab = expression(beta[t]^{unemployment}), ylab = "Density")
hist(post_final_spec2_u, breaks = 40, probability = TRUE,
     col = alpha("darkorange", 0.4), border = NA, add = TRUE)
lines(dens_b_99,  col = "navyblue",   lwd = 2.5)
lines(dens_b_fin, col = "darkorange", lwd = 2.5)
abline(v = 0, col = "grey40", lty = 3, lwd = 1.5)
legend("topright",
       legend = c("1999 Q1", "Final Quarter"),
       fill = c(alpha("navyblue", 0.4), alpha("darkorange", 0.4)),
       border = NA, bty = "n")

# --- TVP-SV ---
dens_sv_99  <- density(post_1999q1_sv)
dens_sv_fin <- density(post_final_sv)
ylim_sv <- range(c(dens_sv_99$y, dens_sv_fin$y))

hist(post_1999q1_sv, breaks = 40, probability = TRUE,
     col = alpha("navyblue", 0.4), border = NA,
     xlim = xlim_all, ylim = ylim_sv,
     main = "TVP-SV",
     xlab = expression(beta[t]^{unemployment}), ylab = "Density")
hist(post_final_sv, breaks = 40, probability = TRUE,
     col = alpha("darkorange", 0.4), border = NA, add = TRUE)
lines(dens_sv_99,  col = "navyblue",   lwd = 2.5)
lines(dens_sv_fin, col = "darkorange", lwd = 2.5)
abline(v = 0, col = "grey40", lty = 3, lwd = 1.5)
legend("topright",
       legend = c("1999 Q1", "Final Quarter"),
       fill = c(alpha("navyblue", 0.4), alpha("darkorange", 0.4)),
       border = NA, bty = "n")

mtext("Unemployment Coefficient Posterior: Baseline TVP vs TVP-SV (Spec 2)",
      outer = TRUE, cex = 1.1, font = 2, line = 1)

# ================================================================
# Posterior probabilities — comparison table
# ================================================================
prob_sv <- mean(abs(post_final_sv) < abs(post_1999q1_sv))

cat("===========================================\n")
cat("POSTERIOR PROBABILITY OF FLATTENING\n")
cat("P(|beta_final| < |beta_1999Q1| | data)\n")
cat("===========================================\n")
cat(sprintf("Baseline TVP  - Spec 2 (Unemployment): %.4f\n", prob_flatten_u))
cat(sprintf("TVP-SV        - Spec 2 (Unemployment): %.4f\n", prob_sv))
cat("===========================================\n\n")

# Summary statistics comparison
cat("--- Posterior means at 1999 Q1 ---\n")
cat(sprintf("  Baseline TVP:  %.4f\n", mean(post_1999q1_spec2_u)))
cat(sprintf("  TVP-SV:        %.4f\n", mean(post_1999q1_sv)))

cat("\n--- Posterior means at Final Quarter ---\n")
cat(sprintf("  Baseline TVP:  %.4f\n", mean(post_final_spec2_u)))
cat(sprintf("  TVP-SV:        %.4f\n", mean(post_final_sv)))

cat("\n--- Posterior SD at 1999 Q1 ---\n")
cat(sprintf("  Baseline TVP:  %.4f\n", sd(post_1999q1_spec2_u)))
cat(sprintf("  TVP-SV:        %.4f\n", sd(post_1999q1_sv)))

cat("\n--- Posterior SD at Final Quarter ---\n")
cat(sprintf("  Baseline TVP:  %.4f\n", sd(post_final_spec2_u)))
cat(sprintf("  TVP-SV:        %.4f\n", sd(post_final_sv)))

cat("\n--- 95% Credible Intervals ---\n")
cat("1999 Q1:\n")
cat(sprintf("  Baseline TVP:  [%.4f, %.4f]\n",
            quantile(post_1999q1_spec2_u, 0.025),
            quantile(post_1999q1_spec2_u, 0.975)))
cat(sprintf("  TVP-SV:        [%.4f, %.4f]\n",
            quantile(post_1999q1_sv, 0.025),
            quantile(post_1999q1_sv, 0.975)))
cat("Final Quarter:\n")
cat(sprintf("  Baseline TVP:  [%.4f, %.4f]\n",
            quantile(post_final_spec2_u, 0.025),
            quantile(post_final_spec2_u, 0.975)))
cat(sprintf("  TVP-SV:        [%.4f, %.4f]\n",
            quantile(post_final_sv, 0.025),
            quantile(post_final_sv, 0.975)))

cat("\n--- Phi^2 (log-volatility innovation variance) ---\n")
cat(sprintf("  Posterior mean: %.4f\n", mean(phi2_store_sv)))
cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
            quantile(phi2_store_sv, 0.025),
            quantile(phi2_store_sv, 0.975)))
##
##END
##
