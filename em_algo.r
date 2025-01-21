# Charger les bibliothèques nécessaires
library(MASS)

# Définir les paramètres du modèle
set.seed(123)
n <- 365 * 3  # Nombre de jours simulés (3 ans de données)
states <- 2  # Nombre d'états
degree_trans_pol <- 2 # Degree of trigonometric covariates
degree_obs_pol <- 2 # Degree of trigonometric covariates
time <- 1:n  # Temps
period <- 365  # Période saisonnière

# Simuler les états selon une chaîne de Markov avec transition saisonnière
simulate_states <- function(n, period, beta, K, degree) {
  Q <- function(t) {
    P <- matrix(0, nrow = K, ncol = K)
    for (i in 1:K) {
      for (j in 1:(K-1)) {
        P[i, j] <- exp(beta[(i-1)*(K-1)*degree*2 + (j-1)*degree*2 + 1] + 
                       sum(sapply(1:degree, function(d) {
                         beta[(i-1)*(K-1)*degree*2 + (j-1)*degree*2 + 2*d] * cos(2 * pi * d * t / period) + 
                         beta[(i-1)*(K-1)*degree*2 + (j-1)*degree*2 + 2*d + 1] * sin(2 * pi * d * t / period)
                       })))
      }
    }
    P[, K] <- 1  # Set the last column to 1 for normalization
    Q <- P / rowSums(P)  # Normalize the rows to sum to 1
    return(Q)
  }
  
  states <- numeric(n)
  states[1] <- sample(1:K, size = 1)
  for (t in 2:n) {
    trans_probs <- Q(t %% period)[states[t - 1], ]
    if (any(is.na(trans_probs))) {
      stop("NA values in transition probabilities")
    }
    states[t] <- sample(1:K, size = 1, prob = trans_probs)
  }
  return(states)
}

# Simuler les observations
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

# Initialisation des paramètres
K <- states
degree <- degree_trans_pol
beta <- runif(states * (states - 1) * degree * 2)  # Coefficients aléatoires pour l'exemple
mu <- runif(K, -1, 1)
delta <- matrix(runif(K * degree_obs_pol * 2, -1, 1), nrow = K)
sd <- runif(K, 0.5, 1.5)

# Simuler les états et les observations
simulated_states <- simulate_states(n, period, beta, states, degree)
observations <- simulate_observations(simulated_states, n, period, mu, delta, sd, degree_obs_pol)

# Algorithme EM pour inférer les paramètres
em_algorithm <- function(observations, K, degree, period, max_iter = 100, tol = 1e-6) {
  n <- length(observations)
  time <- 1:n
  
  # Initialisation des paramètres
  beta <- runif(K * (K - 1) * degree * 2)
  mu <- runif(K, -1, 1)
  delta <- matrix(runif(K * degree * 2, -1, 1), nrow = K)
  sd <- runif(K, 0.5, 1.5)
  
  log_likelihood <- function() {
    # Calculer la log-vraisemblance
    ll <- 0
    for (t in 1:n) {
      k <- states[t]
      mean_t <- mu[k]
      for (d in 1:degree) {
        mean_t <- mean_t + delta[k, 2*d-1] * cos(2 * pi * d * t / period) + delta[k, 2*d] * sin(2 * pi * d * t / period)
      }
      ll <- ll + dnorm(observations[t], mean = mean_t, sd = sd[k], log = TRUE)
    }
    return(ll)
  }
  
  prev_ll <- -Inf  # Initialiser prev_ll à une valeur très basse
  
  for (iter in 1:max_iter) {
    # E-step: Calculer les probabilités a posteriori des états cachés
    gamma <- matrix(0, nrow = n, ncol = K)
    for (t in 1:n) {
      for (k in 1:K) {
        mean_t <- mu[k]
        for (d in 1:degree) {
          mean_t <- mean_t + delta[k, 2*d-1] * cos(2 * pi * d * t / period) + delta[k, 2*d] * sin(2 * pi * d * t / period)
        }
        gamma[t, k] <- dnorm(observations[t], mean = mean_t, sd = sd[k])
      }
      gamma[t, ] <- gamma[t, ] / sum(gamma[t, ])
    }
    
    # M-step: Mettre à jour les paramètres
    for (k in 1:K) {
      mu[k] <- sum(gamma[, k] * observations) / sum(gamma[, k])
      for (d in 1:degree) {
        delta[k, 2*d-1] <- sum(gamma[, k] * observations * cos(2 * pi * d * time / period)) / sum(gamma[, k])
        delta[k, 2*d] <- sum(gamma[, k] * observations * sin(2 * pi * d * time / period)) / sum(gamma[, k])
      }
      sd[k] <- sqrt(sum(gamma[, k] * (observations - mu[k])^2) / sum(gamma[, k]))
    }
    
    # Vérifier la convergence
    ll <- log_likelihood()
    if (iter > 1 && abs(ll - prev_ll) < tol) {
      break
    }
    prev_ll <- ll
  }
  
  return(list(mu = mu, delta = delta, sd = sd, beta = beta))
}

# Appliquer l'algorithme EM
params <- em_algorithm(observations, K, degree, period)

# Afficher les paramètres inférés
print(params)