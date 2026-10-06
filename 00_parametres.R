
# Paramètres du groupe 15 et du produit
# =====================================================================

# Paramètres du modèle de B&S (Groupe 15)
# Dynamique du sous-jacent sous la mesure risque-neutre :
#   dS = r * S * dt + sigma * S * dW

S0 <- 105      # prix initial du sous-jacent
r <- 0.08      # rfr
sigma <- 0.35  # volatilité

# --- Paramètres des produits -----------------------------------------

K <- 110       # prix d'exercice du put européen
T <- 0.5       # maturité en années (6 mois)
