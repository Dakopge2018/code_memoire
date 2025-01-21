# Charger la bibliothèque nécessaire
library(depmixS4)

# Définir les paramètres du modèle
set.seed(123)
n <- 365 * 3  # Nombre de jours simulés (3 ans de données)
states <- 2  # Nombre d'états
degree_trans_pol <- 2 # Degree of trigonometric covariates
degree_obs_pol <- 2 # Degree of trigonometric covariates
time = 1:n  # Temps
period <- 365  # Période saisonnière

# Simuler les états selon une chaîne de Markov avec transition saisonnière
simulate_states <- function(n, period, beta, K) {
  Q <- function(t) {
    P <- matrix(0, nrow = K, ncol = K)
    for (i in 1:K) {
      for (j in 1:(K-1)) {
        P[i, j] <- (exp(beta[(i-1)*(K-1) + j] + beta[(i-1)*(K-1) + j + K] * cos(2 * pi * t / period) + beta[(i-1)*(K-1) + j + 2*K] * sin(2 * pi * t / period)))
      }
    }
    Q <- P / (1 + rowSums(P))
    Q[, K] <- 1 / (1 + rowSums(P))
    return(Q)
  }
  states <- numeric(n)
  states[1] <- sample(1:2, size = 1)
  for (t in 2:n) {
    trans_probs <- Q(t %% period)[states[t - 1], ]
    states[t] <- sample(1:2, size = 1, prob = trans_probs)
  }
  return(states)
}

# Paramètres pour la transition saisonnière
# beta <- c(1, 0.7, 0.5, -1, -0.6, 0.7)
beta <- c(1, 0.7, 0.5, -1.2, 0.4, 0.2, -0.8, 0.3, -0.2, -1.0)
states <- simulate_states(n, period, beta,2)

# Simuler les observations avec \( m_k(t) \)
simulate_observations <- function(states, n, period, mu, delta, sd, degree) {
  observations <- numeric(n)
  for (t in 1:n) {
    k <- states[t]
    mean_t <- mu[k]
    for (d in 1:degree) {
      mean_t <- mean_t + delta[k, 2*d-1] * cos(2 * pi * d * t / period) + delta[k, 2*d] * sin(2 * pi * d * t / period)
    }
    observations[t] <- rnorm(1, mean = mean_t, sd = sd[k])
  }
  return(observations)
}

# Generalized function to compute trigonometric covariates
generate_trig_covariates <- function(time, period, degree) {
  trig_covs <- data.frame(time = time)
  for (d in 1:degree) {
    trig_covs[[paste0("sin_", d)]] <- sin(2 * pi * d * time / period)
    trig_covs[[paste0("cos_", d)]] <- cos(2 * pi * d * time / period)
  }
  return(trig_covs)
}

# Paramètres pour les émissions
mu <- c(-1, 2)  # Moyennes des états
delta <- matrix(c(5, 4, -3, -2,3,8,3,5), nrow = 2, byrow = TRUE)  # Coefficients pour cos et sin
sd <- c(1, 0.5)  # Écarts-types

# Ajuster un modèle SHMM avec DepMixS4
trig_covs <- generate_trig_covariates(time, period, degree_trans_pol)


# Simulated response data
observations <- simulate_observations(states, n, period, mu, delta, sd, degree_obs_pol)
data <- cbind(
  obs = observations,
  trig_covs
)

# Define transition formula dynamically
#transition_formula <- as.formula(
#  paste("~", paste(names(trig_covs)[-1], collapse = " + "))
# )
transition_formulas <- lapply(1:states, function(i) {
  as.formula(paste("~", paste(names(trig_covs)[-1], collapse = " + ")))
})
# Define obs formula dynamically
obs_formula <- as.formula(
  paste("obs ~", paste(names(trig_covs)[-1], collapse = " + "))
)


mod <- depmix(
  response = obs_formula,
  data = data,
  nstates = 2,
  family = gaussian(),
  transition = transition_formula,
)

fitted_mod <- fit(mod, verbose = TRUE)

# Résultats
summary(fitted_mod, which='transition')
summary(fitted_mod, which='response')

transition_formulas <- lapply(1:states, function(i) {
  as.formula(paste("~", paste(names(trig_covs)[-1], collapse = " + ")))
})