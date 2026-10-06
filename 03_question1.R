
# QUESTION 1 : le PUT EUROPÉEN (K = 110, T = 0.5)
#   Q1a : impact des paramètres S0, sigma et r sur le prix
#   Q1b : comparaison Monte-Carlo / formule fermée aux paramètres du groupe
#   Q1c : convergence, en nombre de trajectoires N puis en pas de temps dt
#
# schéma "exact", n = 26, N = 100000.

source("00_parametres.R")
source("01_fonctions.R")

set.seed(123)

# réglages communs à tout le script
N <- 100000   # nombre de trajectoires
n <- 26       # nombre de pas de temps (26 semaines sur 6 mois)

# prix de référence donné par la formule fermée, aux paramètres du groupe
prix_bs <- prix_put_bs(S0, r, sigma, K, T)

cat("=====================================================\n")
cat(" QUESTION 1 - PUT européen - GROUPE 15\n")
cat("=====================================================\n")
cat("S0 =", S0, "  r =", r, "  sigma =", sigma, "  K =", K, "  T =", T, "\n\n")


# =====================================================================
# Q1a : l'impact des para sur le prix du put

# on fix à chaque fois les paramètres sauf celui qu'on test et puis on run MC et on compare le prix obtenu avec la formule théorique
#


cat("--- Q1a : impact des paramètres -----------------------\n")


# Q1a (i) : le prix initial S0 ---------------------------------------------------------------------------------

valeurs_S0 <- seq(80, 130, by = 5)

prix_mc_S0 <- rep(0, length(valeurs_S0))
inf_S0 <- rep(0, length(valeurs_S0))
sup_S0 <- rep(0, length(valeurs_S0))

for (i in 1:length(valeurs_S0)) {

  set.seed(123)   # aléa commun ->> mêmes tirages pour tous les points
  trajectoires <- simuler_trajectoires(valeurs_S0[i], r, sigma, T, n, N, "exact")

  resultat <- prix_put_mc(trajectoires, r, K, T)
  prix_mc_S0[i] <- resultat$prix
  inf_S0[i] <- resultat$borne_inf_ic
  sup_S0[i] <- resultat$borne_sup_ic
}

# courbe théorique sur une grille assez fine pour obtenir un trait lisse.
# //prix_put_bs accepte directement un vecteur : log(), pnorm() et exp()
# puis vectorisées en R, le calcul est donc fait pour les 200 valeurs.
grille_S0 <- seq(80, 130, length.out = 200)
courbe_bs_S0 <- prix_put_bs(grille_S0, r, sigma, K, T)

# graph prix start
png("figures/q1a_prix_vs_S0.png", width = 900, height = 650)
plot(valeurs_S0, prix_mc_S0, pch = 19,
     ylim = range(inf_S0, sup_S0, courbe_bs_S0),
     xlab = "Prix initial du sous-jacent S0",
     ylab = "Prix du put européen",
     main = "Q1a - Prix du put en fonction de S0 (sigma = 0.35, r = 0.08)")
arrows(valeurs_S0, inf_S0, valeurs_S0, sup_S0,
       angle = 90, code = 3, length = 0.03)
lines(grille_S0, courbe_bs_S0, lwd = 2)
abline(v = S0, lty = 2)
legend("topright",
       legend = c("Monte-Carlo (IC à 95 %)",
                  "Formule de Black-Scholes",
                  "Valeur du groupe 15 (S0 = 105)"),
       pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 1),
       bty = "n")
invisible(dev.off())   # ///invisible() évite l'affichage de "null device"

cat("figure enregistrée : figures/q1a_prix_vs_S0.png\n")
cat("  le prix du put DÉCROÎT quand S0 augmente : plus le sous-jacent\n")
cat("  est haut, moins le droit de vendre à K = 110 a de valeur.\n")


# --- Q1a (ii) : la volatilité sigma ----------------------------------------------------

valeurs_sigma <- seq(0.10, 0.60, by = 0.05)

prix_mc_sigma <- rep(0, length(valeurs_sigma))
inf_sigma <- rep(0, length(valeurs_sigma))
sup_sigma <- rep(0, length(valeurs_sigma))

for (i in 1:length(valeurs_sigma)) {

  set.seed(123)   # aléa commun
  trajectoires <- simuler_trajectoires(S0, r, valeurs_sigma[i], T, n, N, "exact")

  resultat <- prix_put_mc(trajectoires, r, K, T)
  prix_mc_sigma[i] <- resultat$prix
  inf_sigma[i] <- resultat$borne_inf_ic
  sup_sigma[i] <- resultat$borne_sup_ic
}

grille_sigma <- seq(0.10, 0.60, length.out = 200)
courbe_bs_sigma <- prix_put_bs(S0, r, grille_sigma, K, T)

#graph sigma
png("figures/q1a_prix_vs_sigma.png", width = 900, height = 650)
plot(valeurs_sigma, prix_mc_sigma, pch = 19,
     ylim = range(inf_sigma, sup_sigma, courbe_bs_sigma),
     xlab = "Volatilité sigma",
     ylab = "Prix du put européen",
     main = "Q1a - Prix du put en fonction de sigma (S0 = 105, r = 0.08)")
arrows(valeurs_sigma, inf_sigma, valeurs_sigma, sup_sigma,
       angle = 90, code = 3, length = 0.03)
lines(grille_sigma, courbe_bs_sigma, lwd = 2)
abline(v = sigma, lty = 2)
legend("topleft",
       legend = c("Monte-Carlo (IC à 95 %)",
                  "Formule de Black-Scholes",
                  "Valeur du groupe 15 (sigma = 0.35)"),
       pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 1),
       bty = "n")
invisible(dev.off())   # invisible() évite l'affichage de "null device"

cat("figure enregistrée : figures/q1a_prix_vs_sigma.png\n")
cat("  le prix du put CROÎT avec sigma : plus d'incertitude sur S_T\n")
cat("  augmente la valeur de l'option, dont le payoff est asymétrique.\n")


# --- Q1a (iii) : le taux sans risque r -------------------------------

valeurs_r <- seq(0.00, 0.16, by = 0.01)

prix_mc_r <- rep(0, length(valeurs_r))
inf_r <- rep(0, length(valeurs_r))
sup_r <- rep(0, length(valeurs_r))

for (i in 1:length(valeurs_r)) {

  set.seed(123)   # aléa commun
  trajectoires <- simuler_trajectoires(S0, valeurs_r[i], sigma, T, n, N, "exact")

  # attention : le taux sert ici DEUX fois, dans la dérive de la
  # simulation et dans l'actualisation du payoff
  resultat <- prix_put_mc(trajectoires, valeurs_r[i], K, T)
  prix_mc_r[i] <- resultat$prix
  inf_r[i] <- resultat$borne_inf_ic
  sup_r[i] <- resultat$borne_sup_ic
}

grille_r <- seq(0.00, 0.16, length.out = 200)
courbe_bs_r <- prix_put_bs(S0, grille_r, sigma, K, T)
# graph r
png("figures/q1a_prix_vs_r.png", width = 900, height = 650)
plot(valeurs_r, prix_mc_r, pch = 19,
     ylim = range(inf_r, sup_r, courbe_bs_r),
     xlab = "Taux sans risque r",
     ylab = "Prix du put européen",
     main = "Q1a - Prix du put en fonction de r (S0 = 105, sigma = 0.35)")
arrows(valeurs_r, inf_r, valeurs_r, sup_r,
       angle = 90, code = 3, length = 0.03)
lines(grille_r, courbe_bs_r, lwd = 2)
abline(v = r, lty = 2)
legend("topright",
       legend = c("Monte-Carlo (IC à 95 %)",
                  "Formule de Black-Scholes",
                  "Valeur du groupe 15 (r = 0.08)"),
       pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 1),
       bty = "n")
invisible(dev.off())   # invisible() évite l'affichage de "null device"

cat("figure enregistrée : figures/q1a_prix_vs_r.png\n")
cat("  le prix du put DÉCROÎT quand r augmente : la dérive risque-neutre\n")
cat("  pousse le sous-jacent vers le haut et le strike est plus actualisé.\n\n")


# =====================================================================
# Q1b : MC vs formules
# On simule une seule fois, aux paramètres du groupe 15, et again on confronte l'estimation Monte-Carlo à la valeur exacte de Black-Scholes.

cat("--- Q1b : Monte-Carlo versus formule fermée -----------\n")

set.seed(123)
trajectoires_q1b <- simuler_trajectoires(S0, r, sigma, T, n, N, "exact")
put_mc_q1b <- prix_put_mc(trajectoires_q1b, r, K, T)

# l'IC vaut prix +/- 1.96 * erreur_type : on retrouve l'erreur-type en divisant la demi-largeur par 1.96
demi_largeur_q1b <- put_mc_q1b$borne_sup_ic - put_mc_q1b$prix
erreur_type_q1b <- demi_largeur_q1b / 1.96

ecart_q1b <- put_mc_q1b$prix - prix_bs
ecart_en_erreurs_types <- ecart_q1b / erreur_type_q1b

cat("N =", format(N, scientific = FALSE), " n =", n, " schéma exact\n")
cat("  prix Monte-Carlo      :", sprintf("%.4f", put_mc_q1b$prix), "\n")
cat("  IC à 95 %             : [", sprintf("%.4f", put_mc_q1b$borne_inf_ic),
    ";", sprintf("%.4f", put_mc_q1b$borne_sup_ic), "]\n")
cat("  erreur-type           :", sprintf("%.4f", erreur_type_q1b), "\n")
cat("  prix Black-Scholes    :", sprintf("%.4f", prix_bs), "\n")
cat("  écart MC - BS         :", sprintf("%+.4f", ecart_q1b), "\n")
cat("  écart en erreurs-types:", sprintf("%+.2f", ecart_en_erreurs_types), "\n")

# Conclusion : l'écart entre l'estimation Monte-Carlo et la valeur exacte
# est de l'ordre d'une demi erreur-type, en dessous du seuil de 1.96 qui délimite l'IC à 95 %. La formule fermée tombe donc dans l'IC : 
# les deux méthodes se rejoignent ->> l'écart résiduel est bruit MC
if (prix_bs >= put_mc_q1b$borne_inf_ic && prix_bs <= put_mc_q1b$borne_sup_ic) {
  cat("  => la valeur Black-Scholes est DANS l'IC à 95 % : les deux\n")
  cat("     méthodes concordent, l'écart n'est que du bruit Monte-Carlo.\n\n")
} else {
  cat("  => la valeur Black-Scholes est HORS de l'IC à 95 %.\n\n")
}


# =====================================================================
# Q1c (1) : Convergence N
# =====================================================================
# On augmente N et on observe deux choses : 
# le prix se rapproche de la valeur exacte, 
# l'intervalle de confiance se resserre.

cat("--- Q1c (1) : convergence en N ------------------------\n")

valeurs_N <- c(1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)

prix_selon_N <- rep(0, length(valeurs_N))
inf_selon_N <- rep(0, length(valeurs_N))
sup_selon_N <- rep(0, length(valeurs_N))
demi_largeur_selon_N <- rep(0, length(valeurs_N))

cat(sprintf("%10s %10s %22s %14s\n", "N", "prix MC", "IC à 95 %", "demi-largeur"))

for (i in 1:length(valeurs_N)) {

  set.seed(123)
  trajectoires <- simuler_trajectoires(S0, r, sigma, T, n, valeurs_N[i], "exact")

  resultat <- prix_put_mc(trajectoires, r, K, T)
  prix_selon_N[i] <- resultat$prix
  inf_selon_N[i] <- resultat$borne_inf_ic
  sup_selon_N[i] <- resultat$borne_sup_ic
  demi_largeur_selon_N[i] <- resultat$borne_sup_ic - resultat$prix

  cat(sprintf("%10d %10.4f   [%8.4f ;%8.4f] %14.4f\n",
              valeurs_N[i], prix_selon_N[i], inf_selon_N[i], sup_selon_N[i],
              demi_largeur_selon_N[i]))
}

# --- figure 1 : le prix converge vers la valeur Black-Scholes ---------
# L'axe des abscisses est logarithmique

png("figures/q1c_convergence_N.png", width = 900, height = 650)
plot(valeurs_N, prix_selon_N, pch = 19, log = "x",
     ylim = range(inf_selon_N, sup_selon_N, prix_bs),
     xlab = "Nombre de trajectoires N (échelle logarithmique)",
     ylab = "Prix du put européen",
     main = "Q1c - Convergence du prix Monte-Carlo quand N augmente")
arrows(valeurs_N, inf_selon_N, valeurs_N, sup_selon_N,
       angle = 90, code = 3, length = 0.03)
abline(h = prix_bs, lwd = 2)
legend("topright",
       legend = c("Monte-Carlo (IC à 95 %)", "Prix Black-Scholes (10.7214)"),
       pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2),
       bty = "n")
invisible(dev.off())   # pour éviter l'erreur 

cat("figure enregistrée : figures/q1c_convergence_N.png\n")

# --- figure 2 : la demi-largeur de l'IC décroît en 1/sqrt(N) ------------------------------------------------------------------------------
# On construit une droite de référence qui passe par le premier point et suit exactement 1/sqrt(N), pour vérifier visuellement cette pente.

demi_largeur_reference <- demi_largeur_selon_N[1] * sqrt(valeurs_N[1] / valeurs_N)

png("figures/q1c_demi_largeur_ic.png", width = 900, height = 650)
plot(valeurs_N, demi_largeur_selon_N, pch = 19, log = "xy",
     xlab = "Nombre de trajectoires N (échelle logarithmique)",
     ylab = "Demi-largeur de l'IC à 95 % (échelle logarithmique)",
     main = "Q1c - La précision Monte-Carlo décroît en 1/sqrt(N)")
lines(valeurs_N, demi_largeur_selon_N)
lines(valeurs_N, demi_largeur_reference, lty = 2, lwd = 2)
legend("topright",
       legend = c("Demi-largeur observée",
                  "Référence de pente -1/2 (1/sqrt(N))"),
       pch = c(19, NA), lty = c(1, 2), lwd = c(1, 2),
       bty = "n")
invisible(dev.off())   # invisible() évite l'affichage de "null device" //warning

cat("figure enregistrée : figures/q1c_demi_largeur_ic.png\n")


# ====================================================================================================
# Q1c (2) : convergence n, Schéma exact vs Euler
# On fait varier le nombre de pas n (donc dt = T/n) et on compare les deux schémas de simulation.


cat("------------- Q1c (2) : convergence en dt, exact versus Euler ---------------------\n")

liste_n <- c(2, 6, 13, 26, 52, 180) # les pas n différents

prix_exact_selon_n <- rep(0, length(liste_n))
prix_euler_selon_n <- rep(0, length(liste_n))
moyenne_diff <- rep(0, length(liste_n))
inf_diff <- rep(0, length(liste_n))
sup_diff <- rep(0, length(liste_n))

cat(sprintf("%6s %10s %10s %12s %22s %10s\n",
            "n", "exact", "euler", "moy(d)", "IC apparié à 95 %", "moy(d)*n"))

for (i in 1:length(liste_n)) {

  n_courant <- liste_n[i]

  # schéma exact
  set.seed(123)
  trajectoires_exact <- simuler_trajectoires(S0, r, sigma, T, n_courant, N, "exact")
  S_T_exact <- trajectoires_exact[, n_courant + 1]
  payoffs_exact <- pmax(K - S_T_exact, 0) * exp(-r * T)

  # schéma d'Euler, tjs avec les mêmes tirages gaussiens
  set.seed(123)
  trajectoires_euler <- simuler_trajectoires(S0, r, sigma, T, n_courant, N, "euler")
  S_T_euler <- trajectoires_euler[, n_courant + 1]
  payoffs_euler <- pmax(K - S_T_euler, 0) * exp(-r * T)

  # prix de chaque schéma (mean(payoffs) est exactement ce que renvoie prix_put_mc, on le recalcule ici car on a besoin des payoffs eux-mêmes pas seulement des means ça on fait après
  prix_exact_selon_n[i] <- mean(payoffs_exact)
  prix_euler_selon_n[i] <- mean(payoffs_euler)

  # les différences trajectoire par trajectoire
  d <- payoffs_euler - payoffs_exact
  moyenne_diff[i] <- mean(d)
  demi_largeur_diff <- 1.96 * sd(d) / sqrt(N)
  inf_diff[i] <- moyenne_diff[i] - demi_largeur_diff
  sup_diff[i] <- moyenne_diff[i] + demi_largeur_diff

  cat(sprintf("%6d %10.4f %10.4f %12.5f   [%8.5f ;%8.5f] %10.4f\n",
              n_courant, prix_exact_selon_n[i], prix_euler_selon_n[i],
              moyenne_diff[i], inf_diff[i], sup_diff[i],
              moyenne_diff[i] * n_courant))
}

png("figures/q1c_schemas_exact_euler.png", width = 900, height = 650)
plot(liste_n, prix_exact_selon_n, type = "b", pch = 19,
     ylim = range(prix_exact_selon_n, prix_euler_selon_n, prix_bs),
     xlab = "Nombre de pas de temps n",
     ylab = "Prix du put européen",
     main = "Q1c - Schéma exact versus schéma d'Euler (N = 100000)")
lines(liste_n, prix_euler_selon_n, type = "b", pch = 17, lty = 3)
abline(h = prix_bs, lwd = 2)
legend("bottomright",
       legend = c("Schéma exact", "Schéma d'Euler",
                  "Prix Black-Scholes (10.7214)"),
       pch = c(19, 17, NA), lty = c(1, 3, 1), lwd = c(1, 1, 2),
       bty = "n")
invisible(dev.off())   # invisible() évite l'affichage de "null device"

cat("figure enregistrée : figures/q1c_schemas_exact_euler.png\n\n")
cat("=====================================================\n")
cat(" Fin Q1\n")
cat("=====================================================\n")
