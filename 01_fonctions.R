
# 1. simuler_trajectoires()   : simulation du sous-jacent (exact / Euler)
# 2. prix_put_bs()            : formule fermée du put européen (Q1)
# 3. prix_lookback_formule()  : formule fermée du lookback put (Q2)
# 4. prix_put_mc()            : prix Monte-Carlo du put européen (Q1)
# 5. prix_lookback_mc()       : prix Monte-Carlo du lookback put (Q2)
# =========================================================================================

# 1) Simulation des trahectoires
# ---------------------------------------------------------------------------------------
# On simule N trajectoires du sous-jacent entre 0 et T, en n pas de temps.
#
# Le résultat est une matrice de taille N x (n+1) :
#   - la ligne i    = la trajectoire numéro i
#   - la colonne j  = l'instant t_{j-1} = (j-1) * T/n
# La colonne 1 correspond donc à l'instant 0 et vaut S0 partout avec S0 le prix de S(t) en t=0
# la colonne n+1 correspond à l'instant T et contient les S_T
#
# L'argument schema vaut "exact" ou "euler" :
#
#   - "exact" : on utilise la solution exacte de l'EDS du mouvement
#     brownien géométrique. Entre deux instants séparés de dt :
#         S_{j+1} = S_j * exp( (r - sigma^2/2)*dt + sigma*sqrt(dt)*epsilon )
#     avec epsilon ~ N(0,1). Cette formule est exacte : elle ne fait aucune
#     approximation sur la loi de S_{j+1} sachant S_j.
#
#   - "euler" : on discrétise directement l'EDS dS = r*S*dt + sigma*S*dW
#     en remplaçant dt par un petit pas et dW par sqrt(dt)*epsilon :
#         S_{j+1} = S_j + r*S_j*dt + sigma*S_j*sqrt(dt)*epsilon
#                 = S_j * ( 1 + r*dt + sigma*sqrt(dt)*epsilon )
#     C'est une approximation : elle introduit une erreur qui diminue quand n augmente (donc quand dt diminue et qu'on approche de la version continue).

simuler_trajectoires <- function(S0, r, sigma, T, n, N, schema) {

  # longueur d'un pas de temps
  dt <- T / n

  # matrice vide qui va contenir toutes les trajectoires
  trajectoires <- matrix(0, nrow = N, ncol = n + 1)

  # à l'instant 0, toutes les trajectoires partent de S0
  trajectoires[, 1] <- S0

  # boucle sur les n pas de temps : on remplit la matrice colonne
  # par colonne, de la gauche vers la droite
  for (j in 1:n) {

    # un tirage gaussien par trajectoire, pour le pas de temps courant
    epsilon <- rnorm(N)

    if (schema == "exact") {
      # S_{j+1} = S_j * exp( (r - sigma^2/2)*dt + sigma*sqrt(dt)*epsilon )
      trajectoires[, j + 1] <- trajectoires[, j] *
        exp((r - sigma^2 / 2) * dt + sigma * sqrt(dt) * epsilon)
    } else {
      # S_{j+1} = S_j * ( 1 + r*dt + sigma*sqrt(dt)*epsilon )
      trajectoires[, j + 1] <- trajectoires[, j] *
        (1 + r * dt + sigma * sqrt(dt) * epsilon)
    }
  }

  return(trajectoires)
}


# ---------------------------------------------------------------------
# 2)Prix du Put européen par B&S  (Q1)

# Formule du put européen de prix d'exercice K et de maturité T :
#
#   P = K * exp(-r*T) * N(-d2) - S0 * N(-d1)
#
# avec
#   d1 = ( log(S0/K) + (r + sigma^2/2)*T ) / ( sigma*sqrt(T) )
#   d2 = d1 - sigma*sqrt(T)
#
# et N(.) la fonction de répartition de la loi normale centrée réduite, qui s'obtient en R avec pnorm().

prix_put_bs <- function(S0, r, sigma, K, T) {

  d1 <- (log(S0 / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)

  prix <- K * exp(-r * T) * pnorm(-d2) - S0 * pnorm(-d1)

  return(prix)
}


# -----------------------------------------------------------------------------------------------------------
# 3) Prix du LP par la formuel  (Q2)
# Le lookback put considéré paie, en T, la différence entre le maximum du sous-jacent sur [0, T] et sa valeur finale : max(S) - S_T.
#
# Formule fermée donnée dans l'énoncé :
#
#   LP = S0 * ( exp(-r*T)*N(-c2) - N(-c1) ) + S0 * sigma^2/(2*r) * ( N(c1) - exp(-r*T)*N(-c2) )

# avec
#   c1 = (r + sigma^2/2)*T / ( sigma*sqrt(T) )
#   c2 = (r - sigma^2/2)*T / ( sigma*sqrt(T) )
#
# Cette formule suppose un maximum continu sur [0, T] -> quand on va discrétiser on aura un soucis

prix_lookback_formule <- function(S0, r, sigma, T) {

  c1 <- (r + sigma^2 / 2) * T / (sigma * sqrt(T))
  c2 <- (r - sigma^2 / 2) * T / (sigma * sqrt(T))

  prix <- S0 * (exp(-r * T) * pnorm(-c2) - pnorm(-c1)) +
    S0 * sigma^2 / (2 * r) * (pnorm(c1) - exp(-r * T) * pnorm(-c2))

  return(prix)
}


# ---------------------------------------------------------------------
# 4) Prix du put européen par monte-carlo  (Q1)

# La fonction reçoit en input une matrice de trajectoires déjà simulée.
#
# le prix est l'espérance risque-neutre du payoff actualisé,
# que l'on estime par la moyenne empirique sur les N trajectoires :
#
#   prix = moyenne de  exp(-r*T) * max(K - S_T, 0)
#
# Intervalle de confiance à 95 % (TCL)
#
#   prix +/- 1.96 * ecart_type_empirique / sqrt(N)
#
# L'écart-type est calculé sur les payoffs actualisé 
# (c'est bien de leur moyenne qu'on mesure la précision), 
# avec sd() donc un diviseur N-1.

prix_put_mc <- function(trajectoires, r, K, T) {

  # nombre de trajectoires et nombre de colonnes, 
  # déduits de la matrice
 
  N <- nrow(trajectoires)
  nb_colonnes <- ncol(trajectoires)

  # valeurs finales du sous-jacent en dernière colonne de la matrice
  S_T <- trajectoires[, nb_colonnes]

  # payoff du put européen en T, puis actualisation en 0
  payoffs <- pmax(K - S_T, 0)
  payoffs_actualises <- exp(-r * T) * payoffs

  # estimation Monte-Carlo et IC à 95 %
  prix <- mean(payoffs_actualises)
  demi_largeur <- 1.96 * sd(payoffs_actualises) / sqrt(N)

  borne_inf_ic <- prix - demi_largeur
  borne_sup_ic <- prix + demi_largeur

  resultat <- list(prix = prix,
                   borne_inf_ic = borne_inf_ic,
                   borne_sup_ic = borne_sup_ic)

  return(resultat)
}


# ---------------------------------------------------------------------
# 5) Prix du LP par Monte-Carlo (Q2)

# Même principe que pour le put européen, mais le payoff dépend de toute la trajectoire :
#
#   prix = moyenne de  exp(-r*T) * ( max(S) - S_T )
#
# ATTENTION : le max est ici un max sur l'un des point de la grille temporelle qui est granulée de pas, 
# pris sur les n+1 instants simulés seulement. Il est donc inférieur ou égal au vrai maximum continu, 
# donc le prix Monte-Carlo est systématiquement plus bas que celui de la formule fermée. 

prix_lookback_mc <- function(trajectoires, r, T) {

  N <- nrow(trajectoires)
  nb_colonnes <- ncol(trajectoires)

  S_T <- trajectoires[, nb_colonnes]

  # maximum de chaque trajectoire : on parcourt les lignes une par une
  # (la colonne 1, qui vaut S0, est incluse dans le maximum)-> si jamais l'action ne fait que chuter
  maximums <- rep(0, N)
  for (i in 1:N) {
    maximums[i] <- max(trajectoires[i, ])
  }

  # payoff du lookback put en T.
  # le maximum est pris sur toute la trajectoire, dernière colonne S_T comprise, 
  # donc on a toujours maximums >= S_T et le payoff est automatiquement >= 0.
  payoffs <- maximums - S_T
  payoffs_actualises <- exp(-r * T) * payoffs

  prix <- mean(payoffs_actualises)
  demi_largeur <- 1.96 * sd(payoffs_actualises) / sqrt(N)

  borne_inf_ic <- prix - demi_largeur
  borne_sup_ic <- prix + demi_largeur

  resultat <- list(prix = prix,
                   borne_inf_ic = borne_inf_ic,
                   borne_sup_ic = borne_sup_ic)

  return(resultat)
}
