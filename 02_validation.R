# Ce script vérifie que l'engin de simulation et les formules du fichier 01_fonctions.R sont corrects, avant de répondre aux questions
# =====================================================================================================================================

# chargement des paramètres du groupe et des fonctions
source("00_parametres.R")
source("01_fonctions.R")

# seed fixée 
set.seed(123)

# taille de l'expérience de validation
N <- 100000   # nombre de trajectoires
n <- 26       # nombre de pas de temps (26 semaines sur 6 mois)

cat("=====================================================\n")
cat(" Validation pricing engine - GROUPE 15\n")
cat("=====================================================\n")
cat("S0 =", S0, "  r =", r, "  sigma =", sigma, "  K =", K, "  T =", T, "\n")
cat("N  =", format(N, scientific = FALSE), "trajectoires,  n =", n, "pas de temps\n\n")


# TEST 1 : est-ce qu'on a bien une dynamique de prix sous Q martingales ?
# -------------------------------------------------------------------------------
# Sous la mesure risque-neutre, le prix actualisé du sous-jacent est une martingale, donc  E[ exp(-r*T) * S_T ] = S0.
# Si la simulation est correcte, la moyenne empirique doit être proche de S0, à l'erreur Monte-Carlo près.

cat("--- TEST 1 : check martingale (schéma exact) ---------\n")

trajectoires_exact <- simuler_trajectoires(S0, r, sigma, T, n, N, "exact")
S_T_exact <- trajectoires_exact[, n + 1]

# valeurs actualisées du sous-jacent en T
valeurs_actualisees <- exp(-r * T) * S_T_exact

moyenne_actualisee <- mean(valeurs_actualisees)
ecart <- moyenne_actualisee - S0

# erreur-type empirique de la moyenne : elle mesure l'ampleur normal du bruit Monte-Carlo sur cette estimation
erreur_type <- sd(valeurs_actualisees) / sqrt(N)

# on compare l'écart à cette erreur-type : le rapport doit être petit
ratio <- ecart / erreur_type

cat("moyenne de exp(-r*T)*S_T :", sprintf("%.4f", moyenne_actualisee), "\n")
cat("valeur attendue S0       :", sprintf("%.4f", S0), "\n")
cat("écart                    :", sprintf("%.4f", ecart), "\n")
cat("erreur-type Monte-Carlo  :", sprintf("%.4f", erreur_type), "\n")
cat("écart / erreur-type      :", sprintf("%.4f", ratio), "\n")

if (abs(ratio) < 4) {
  cat("=> OK : l'écart tient dans le bruit Monte-Carlo.\n\n")
} else {
  cat("=> ECHEC : l'écart est trop grand pour être du hasard.\n\n")
}


# ---------------------------------------------------------------------
# TEST 2 : Formules théoriques
# On vérifie que les deux formules redonnent bien les valeurs de référence du groupe 15.

cat("------------------------- TEST 2 : formules fermées ------------------------\n")

put_bs <- prix_put_bs(S0, r, sigma, K, T)
lookback_formule <- prix_lookback_formule(S0, r, sigma, T)

cat("put européen (Black-Scholes) :", sprintf("%.4f", put_bs), " | attendu : 10.7214\n")
if (abs(put_bs - 10.7214) < 1e-4) {
  cat("=> OK\n")
} else {
  cat("=> ECHEC\n")
}

cat("lookback put (formule)       :", sprintf("%.4f", lookback_formule), " | attendu : 19.9832\n")
if (abs(lookback_formule - 19.9832) < 1e-4) {
  cat("=> OK\n\n")
} else {
  cat("=> ECHEC\n\n")
}


# ---------------------------------------------------------------------
# TEST 3 : Prix MC put européen (schéma exact)

# On réutilise les trajectoires du test 1. 
# Le prix Monte-Carlo est une estimation du vrai prix 
# son intervalle de confiance à 95 % doit contenir la valeur exacte donnée par Black-Scholes.

cat("------ TEST 3 : put européen Monte-Carlo (exact) --------\n")

put_mc_exact <- prix_put_mc(trajectoires_exact, r, K, T)

cat("prix Monte-Carlo :", sprintf("%.4f", put_mc_exact$prix), "\n")
cat("IC à 95 %        : [", sprintf("%.4f", put_mc_exact$borne_inf_ic), ";",
    sprintf("%.4f", put_mc_exact$borne_sup_ic), "]\n")
cat("prix Black-Scholes :", sprintf("%.4f", put_bs), "\n")

if (put_bs >= put_mc_exact$borne_inf_ic && put_bs <= put_mc_exact$borne_sup_ic) {
  cat("=> OK : la valeur Black-Scholes est dans l'IC à 95 %.\n\n")
} else {
  cat("=> ECHEC : la valeur Black-Scholes est hors de l'IC.\n\n")
}


# ---------------------------------------------------------------------
# TEST 4 : PRIX MONTE-CARLO DU LOOKBACK PUT (schéma exact)
# Toujours les mêmes trajectoires, avec n = 26 pas.
# Le maximum simulé n'est observé qu'aux 27 instants de la grille, 
# alors que la formule fermée suppose un maximum continu 
# le prix Monte-Carlo doit donc être inférieur à 19.9832. 
# -->>> Ce biais de discrétisation est attendu >> pas corriger

cat("----- TEST 4 : lookback put Monte-Carlo (exact) --------\n")

lookback_mc_exact <- prix_lookback_mc(trajectoires_exact, r, T)

cat("prix Monte-Carlo :", sprintf("%.4f", lookback_mc_exact$prix), "\n")
cat("IC à 95 %        : [", sprintf("%.4f", lookback_mc_exact$borne_inf_ic), ";", sprintf("%.4f", lookback_mc_exact$borne_sup_ic), "]\n")
cat("prix formule fermée (maximum continu) :", sprintf("%.4f", lookback_formule), "\n")

if (lookback_mc_exact$prix < 19.98) {
  cat("=> OK : prix inférieur à 19.98, biais du maximum discret.\n\n")
} else {
  cat("=> ECHEC : le biais de discrétisation attendu n'apparaît pas.\n\n")
}


# ---------------------------------------------------------------------------------------------------------------------------
# Comparaison Exact§ vs Euler discrétisatrion

# On refait les tests 3 et 4 avec le schéma d'Euler.
#
# On remet volontairement set.seed(123) juste avant : les deux schémas consomment ainsi les mêmes alés. 
# L'intérêt est que la différence entre les deux prix ne vient alors que de la discrétisation du schéma, et non du bruit MC

cat("--- COMPARAISON : schéma exact contre schéma d'Euler -\n")
cat("(mêmes tirages gaussiens pour les deux schémas)\n\n")

set.seed(123)
trajectoires_euler <- simuler_trajectoires(S0, r, sigma, T, n, N, "euler")

put_mc_euler <- prix_put_mc(trajectoires_euler, r, K, T)
lookback_mc_euler <- prix_lookback_mc(trajectoires_euler, r, T)

cat("PUT EUROPÉEN\n")
cat("  schéma exact :", sprintf("%.4f", put_mc_exact$prix), " IC [", sprintf("%.4f", put_mc_exact$borne_inf_ic), ";",
    sprintf("%.4f", put_mc_exact$borne_sup_ic), "]\n")
cat("  schéma Euler :", sprintf("%.4f", put_mc_euler$prix), " IC [", sprintf("%.4f", put_mc_euler$borne_inf_ic), ";",
    sprintf("%.4f", put_mc_euler$borne_sup_ic), "]\n")
cat("  Black-Scholes:", sprintf("%.4f", put_bs), "\n")
cat("  écart Euler - exact :",
    sprintf("%.4f", put_mc_euler$prix - put_mc_exact$prix), "\n")

if (put_bs >= put_mc_euler$borne_inf_ic && put_bs <= put_mc_euler$borne_sup_ic) {
  cat("  => OK : la valeur Black-Scholes est dans l'IC du schéma d'Euler.\n\n")
} else {
  cat("  => ECHEC : la valeur Black-Scholes est hors de l'IC d'Euler.\n\n")
}

cat("LOOKBACK PUT\n")
cat("  schéma exact :", sprintf("%.4f", lookback_mc_exact$prix), " IC [", sprintf("%.4f", lookback_mc_exact$borne_inf_ic), ";",
    sprintf("%.4f", lookback_mc_exact$borne_sup_ic), "]\n")
cat("  schéma Euler :", sprintf("%.4f", lookback_mc_euler$prix), " IC [", sprintf("%.4f", lookback_mc_euler$borne_inf_ic), ";",
    sprintf("%.4f", lookback_mc_euler$borne_sup_ic), "]\n")
cat("  formule fermée:", sprintf("%.4f", lookback_formule), "\n")
cat("  écart Euler - exact :",
    sprintf("%.4f", lookback_mc_euler$prix - lookback_mc_exact$prix), "\n")

if (lookback_mc_euler$prix < 19.98) {
  cat("  => OK : prix inférieur à 19.98, même biais de discrétisation.\n\n")
} else {
  cat("  => ECHEC : le biais de discrétisation attendu n'apparaît pas.\n\n")
}

cat("=====================================================\n")
cat(" END validation\n")
cat("=====================================================\n")
