
# QUESTION 2 : le LOOKBACK PUT (T = 0.5), de payoff  max(S) - S_T
#   Q2a : prix Monte-Carlo au pas hebdomadaire
#   Q2b : stabilité et convergence, en N puis en n
#   Q2c : taille d'échantillon justifiée par une précision cible
#   Q2d : comparaison avec la formule fermée
#   Q2e : comparaison avec le put européen
#
# SCHÉMA "exact" partout : le biais étudié dans cette question n'est pas celui d'Euler, mais celui du max discret. 
# La formule fermée suppose un maximum continu sur [0, T], alors que la simulation ne l'observe qu'aux n+1 instants de la grille temporelle.
#



# =====================================================================

source("00_parametres.R")
source("01_fonctions.R")

set.seed(123)

# prix de référence donné par la formule fermée de l'énoncé
prix_formule <- prix_lookback_formule(S0, r, sigma, T)

cat("=====================================================\n")
cat(" QUESTION 2 - Lookback Put \n")
cat("=====================================================\n")
cat("S0 =", S0, "  r =", r, "  sigma =", sigma, "  T =", T, "\n")
cat("formule fermée (maximum continu) :", sprintf("%.4f", prix_formule), "\n\n")


# =====================================================================
# ----------Q2a : Prix MC pas hebdo
# 26 pas sur 6 mois, soit un relevé par semaine.

cat("--- Q2a : prix au pas hebdomadaire --------------------\n")

N <- 100000
n <- 26

set.seed(123)
trajectoires_q2a <- simuler_trajectoires(S0, r, sigma, T, n, N, "exact")
lookback_q2a <- prix_lookback_mc(trajectoires_q2a, r, T)

cat("N =", format(N, scientific = FALSE), " n =", n, " schéma exact\n")
cat("  prix Monte-Carlo :", sprintf("%.4f", lookback_q2a$prix), "\n")
cat("  IC à 95 %        : [", sprintf("%.4f", lookback_q2a$borne_inf_ic), ";", sprintf("%.4f", lookback_q2a$borne_sup_ic), "]\n")

# proportion de trajectoires dont le payoff est exactement nul.
# On recalcule les maximums ligne par ligne, comme dans prix_lookback_mc.
S_T_q2a <- trajectoires_q2a[, n + 1]
maximums_q2a <- rep(0, N)
for (i in 1:N) {
  maximums_q2a[i] <- max(trajectoires_q2a[i, ])
}
proportion_nulle_lookback <- mean((maximums_q2a - S_T_q2a) == 0)

# >CHECK : ordre de grandeur donné par la loi arcsinus discrète pour la position du maximum d'une marche symétrique. 
# La probabilité que le maximum soit atteint au tout dernier instant est d'ordre 1/sqrt(pi*n).
controle_arcsinus <- 1 / sqrt(pi * n)

cat("  proportion de payoff NUL :", sprintf("%.4f", proportion_nulle_lookback), "\n")
cat("  contrôle 1/sqrt(pi*n)    :", sprintf("%.4f", controle_arcsinus), "\n")
cat("  => le lookback paie sur", sprintf("%.2f", 100*(1-proportion_nulle_lookback)), "% des scénarios.\n\n")

# ///
# Le payoff max - S_T est toujours >= 0, car le maximum est pris sur touite la trajectoire, 
# dernière colonne S_T comprise. Il n'y a donc pas d'exercice de l'option.
#
# Mais max atteint au dernier point pour environ 11.5 % des trajectoires (celles qui finissent sur leur all time high), et le payoff y vaut exactement 0.
# En temps continu, la probabilité que le maximum tombe exactement en T serait pratiquement nulle
# L'ordre de grandeur observé est bien celui de la loi arcsinus discrète, en 1/sqrt(pi*26) = 0.1106. 
# Le léger excès observé vient de la dérive positive (r - sigma^2/2) > 0, qui pousse les trajectoires à finir plus haut, donc plus souvent sur leur maximum.


# ======================================================================================================================================================================
# Q2b (1) : Convergence N

cat("--- Q2b (1) : convergence en N (n = 26) ---------------\n")

valeurs_N <- c(1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000)

prix_selon_N <- rep(0, length(valeurs_N))
inf_selon_N <- rep(0, length(valeurs_N))
sup_selon_N <- rep(0, length(valeurs_N))

cat(sprintf("%10s %10s %22s %14s\n", "N", "prix MC", "IC à 95 %", "demi-largeur"))

for (i in 1:length(valeurs_N)) {

  set.seed(123)
  trajectoires <- simuler_trajectoires(S0, r, sigma, T, 26, valeurs_N[i], "exact")

  resultat <- prix_lookback_mc(trajectoires, r, T)
  prix_selon_N[i] <- resultat$prix
  inf_selon_N[i] <- resultat$borne_inf_ic
  sup_selon_N[i] <- resultat$borne_sup_ic

  cat(sprintf("%10d %10.4f   [%8.4f ;%8.4f] %14.4f\n",
              valeurs_N[i], prix_selon_N[i], inf_selon_N[i], sup_selon_N[i],
              sup_selon_N[i] - prix_selon_N[i]))
}

# vers quoi le prix se stabilise-t-il réellement, à n = 26 ?
limite_empirique <- mean(prix_selon_N[7:9])
cat("limite empirique (moyenne des 3 derniers points) :", sprintf("%.4f", limite_empirique), "\n")

png("figures/q2b_convergence_trajectoires.png", width = 900, height = 650)
plot(valeurs_N, prix_selon_N, pch = 19, log = "x",
     ylim = range(inf_selon_N, sup_selon_N, prix_formule),
     xlab = "Nombre de trajectoires N (échelle logarithmique)",
     ylab = "Prix du lookback put",
     main = "Q2b - Augmenter N réduit la variance, pas le biais")
arrows(valeurs_N, inf_selon_N, valeurs_N, sup_selon_N,
       angle = 90, code = 3, length = 0.03)
abline(h = prix_formule, lwd = 2)
abline(h = limite_empirique, lty = 2, lwd = 2)
legend("right",
       legend = c("Monte-Carlo (IC à 95 %)",
                  "Formule fermée (19.9832)",
                  "Limite empirique à n = 26"),
       pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 2, 2),
       bty = "n")
invisible(dev.off())

cat("figure enregistrée : figures/q2b_convergence_trajectoires.png\n")

# L'intervalle de confiance rétrécit régulièrement quand N augmente
# de +/- 0.84 à N = 1000 jusqu'à +/- 0.04 à N = 500000.
#  Mais le prix ne converge pas vers la formule fermée 20
# se stabilise vers 16.73. 
#
# Augmenter N réduit la variance de l'estimateur, pas le biais 
# Or ce que l'on estime ici, avec une précision de plus en plus grande, 
# c'est le prix du lookback à mex discret sur 26 pas, 
# ce qui n'est pas le prix du lookback à maximum continu. 
# Pour corriger cela, il faut augmenter n pas N.
cat("  => le prix converge vers ~16.73, PAS vers 19.98.\n\n")


# =====================================================================
# Q2b (2) : Convergence n

cat("--- Q2b (2) : convergence en n (N = 100000) -----------\n")

liste_n <- c(6, 13, 26, 52, 104, 180)

prix_selon_n <- rep(0, length(liste_n))
inf_selon_n <- rep(0, length(liste_n))
sup_selon_n <- rep(0, length(liste_n))
biais <- rep(0, length(liste_n))

cat(sprintf("%6s %10s %22s %10s %14s %12s\n",
            "n", "prix MC", "IC à 95 %", "biais", "biais*sqrt(n)", "biais*n"))

for (i in 1:length(liste_n)) {

  set.seed(123)
  trajectoires <- simuler_trajectoires(S0, r, sigma, T, liste_n[i], 100000, "exact")

  resultat <- prix_lookback_mc(trajectoires, r, T)
  prix_selon_n[i] <- resultat$prix
  inf_selon_n[i] <- resultat$borne_inf_ic
  sup_selon_n[i] <- resultat$borne_sup_ic

  # le biais est l'écart entre la formule (maximum continu) et le prix simulé (maximum discret)
  biais[i] <- prix_formule - prix_selon_n[i]

  cat(sprintf("%6d %10.4f   [%8.4f ;%8.4f] %10.4f %14.4f %12.4f\n",
              liste_n[i], prix_selon_n[i], inf_selon_n[i], sup_selon_n[i],
              biais[i], biais[i] * sqrt(liste_n[i]), biais[i] * liste_n[i]))
}

png("figures/q2b_convergence_pas.png", width = 900, height = 650)
plot(liste_n, prix_selon_n, type = "b", pch = 19,
     ylim = range(inf_selon_n, sup_selon_n, prix_formule),
     xlab = "Nombre de pas de temps n",
     ylab = "Prix du lookback put",
     main = "Q2b - Le prix rejoint la formule fermée par en dessous")
arrows(liste_n, inf_selon_n, liste_n, sup_selon_n,
       angle = 90, code = 3, length = 0.03)
abline(h = prix_formule, lwd = 2)
legend("bottomright",
       legend = c("Monte-Carlo (IC à 95 %)", "Formule fermée (19.9832)"),
       pch = c(19, NA), lty = c(1, 1), lwd = c(1, 2),
       bty = "n")
invisible(dev.off())

cat("figure enregistrée : figures/q2b_convergence_pas.png\n")

# ---sous question:  le biais décroît-il bien en 1/sqrt(n) ? -------------------------
# Droite de référence de pente -1/2, ancrée sur le premier point.
biais_reference <- biais[1] * sqrt(liste_n[1] / liste_n)

png("figures/q2b_biais_loglog.png", width = 900, height = 650)
plot(liste_n, biais, pch = 19, log = "xy",
     xlab = "Nombre de pas de temps n (échelle logarithmique)",
     ylab = "Biais = formule - prix MC (échelle logarithmique)",
     main = "Q2b - Le biais du maximum discret décroît en O(sqrt(dt))")
lines(liste_n, biais)
lines(liste_n, biais_reference, lty = 2, lwd = 2)
legend("topright",
       legend = c("Biais observé", "Référence de pente -1/2 (1/sqrt(n))"),
       pch = c(19, NA), lty = c(1, 2), lwd = c(1, 2),
       bty = "n")
invisible(dev.off())

cat("figure enregistrée : figures/q2b_biais_loglog.png\n\n")

# Le prix monte vers la formule fermée par le dessous: le maximum observé sur une grille discrète est toujours inférieur ou égal au maximum continu, donc le payoff simulé est systématiquement sous-estimé.
#
# Le taux de décroissance est en O(sqrt(dt)), soit 1/sqrt(n). 
# Les deux colonnes du tableau le montrent :
#   - biais*sqrt(n) ne varie que de 15.26 à 17.69, soit +16 % pour un facteur 30 sur n 
#   - biais*n varie, lui, de 37.4 à 237.4, soit un facteur 6.3.
# >> C'est donc bien sqrt(n), et non n, qui est le bon ordre de grandeur.
#
# Lien avec Q1c : Pour le put européen, dont le payoff ne dépend que de S_T,ça avait gardé le biais en dt 
# S_T est tiré dans sa loi exacte, quel que soit le nombre de pas. 
# Pour le lookback, dont le payoff dépend de toute la trajectoire, le biais en dt existe même avec le schéma exact. 
# Chaque S_{t_j} est pourtant simulé exactement ! >> c'est le maximum qui n'est observé qu'aux points de grille et rate systématiquement les pics survenus entre deux relevés.
#
# ---------------------------------------------------------------------
# BONUS: pourquoi biais*sqrt(n) dérive-t-il encore de +16 % au lieu d'être constant ?
#
# Convergence lente vers une constante théorique connue. La correction de continuité du maximum discret 
#
#   ça donne: biais * sqrt(n)  ->  beta * sigma * sqrt(T) * (LP + S0)
#
# où beta ≈ 0.5826 est la constante de la théorie du maximum observé en temps discret, et où LP + S0 = exp(-r*T) * E[max continu]. 
# Cette seconde égalité vient de la propriété de martingale par définition
# LP = exp(-r*T)*E[max - S_T], et exp(-r*T)*E[S_T] = S0.
# Numériquement : LP + S0 = 124.98.
constante_asymptotique <- 0.5826 * sigma * sqrt(T) * (prix_formule + S0)
cat("BONUS - constante asymptotique théorique de biais*sqrt(n) :",
    sprintf("%.4f", constante_asymptotique), "\n")
cat("        valeurs observées à n = 104 et n = 180 :",
    sprintf("%.4f", biais[5]*sqrt(liste_n[5])), "et",
    sprintf("%.4f", biais[6]*sqrt(liste_n[6])), "\n\n")
# Les valeurs observées 17.67 et 17.69 convergent vers 18.02 pas en bas : les +16 % de dérive sont cette convergence, et non un défaut du taux O(sqrt(dt)).


# =====================================================================
# Q2c : Quelle taille d'échantillon ?
# =====================================================================
# Objectif : obtenir un IC à 95 % de demi-largeur au plus 0.05.
#
# La demi-largeur vaut 1.96 * sd / sqrt(N). On veut donc
#   1.96 * sd / sqrt(N) <= 0.05,  c'est-à-dire  N >= (1.96*sd/0.05)^2.
# Comme sd est inconnu, on l'estime par une petite simulation test

cat("--- Q2c : taille d'échantillon justifiée --------------\n")

N_test <- 10000

set.seed(123)
trajectoires_test <- simuler_trajectoires(S0, r, sigma, T, 26, N_test, "exact")
resultat_test <- prix_lookback_mc(trajectoires_test, r, T)

# on retrouve l'écart-type des payoffs à partir de la demi-largeur de l'IC, sans avoir à recalculer les payoffs :
#   demi_largeur = 1.96 * sd / sqrt(N)  donc  sd = demi_largeur*sqrt(N)/1.96
demi_largeur_test <- resultat_test$borne_sup_ic - resultat_test$prix
sd_test <- demi_largeur_test * sqrt(N_test) / 1.96

precision_visee <- 0.05
N_brut <- ceiling((1.96 * sd_test / precision_visee)^2)
N_requis <- ceiling(N_brut / 1000) * 1000   # arrondi au millier supérieur

cat("simulation test : N =", format(N_test, scientific = FALSE), ", n = 26\n")
cat("  prix test            :", sprintf("%.4f", resultat_test$prix), "\n")
cat("  demi-largeur test    :", sprintf("%.4f", demi_largeur_test), "\n")
cat("  écart-type estimé      :", sprintf("%.4f", sd_test), "\n")
cat("précision visée          :", precision_visee, "\n")
cat("  N brut = (1.96*sd/0.05)^2 =", N_brut, "\n")
cat("  N requis (arrondi)       =", format(N_requis, scientific = FALSE), "\n")

set.seed(123)
trajectoires_q2c <- simuler_trajectoires(S0, r, sigma, T, 26, N_requis, "exact")
resultat_q2c <- prix_lookback_mc(trajectoires_q2c, r, T)
demi_largeur_q2c <- resultat_q2c$borne_sup_ic - resultat_q2c$prix

cat("simulation finale : N =", format(N_requis, scientific = FALSE), ", n = 26\n")
cat("  prix Monte-Carlo :", sprintf("%.4f", resultat_q2c$prix), "\n")
cat("  IC à 95 %        : [", sprintf("%.4f", resultat_q2c$borne_inf_ic), ";", sprintf("%.4f", resultat_q2c$borne_sup_ic), "]\n")
cat("  demi-largeur     :", sprintf("%.5f", demi_largeur_q2c), "\n")

if (demi_largeur_q2c <= precision_visee) {
  cat("  => OK : la précision visée de 0.05 est atteinte.\n")
} else {
  cat("  => ECHEC : la précision visée n'est pas atteinte.\n")
}

# Cette précision cible ne gouverne que l'erreur statistique inhérente au MC, c'est-à-dire la largeur de l'intervalle de confiance. 
# c'est pas l'écart avec la vraie valeur : le prix obtenu reste biaisé d'environ -3.25 par le maximum discret, puisqu'on travaille toujours à n = 26.
# Pour se rapprocher de la vraie valeur, il faut augmenter n, pas N et c'est exactement ce que montre Q2b.
cat("  Warning : IC serré ne veut pas dire prix exact. Le biais du\n")
cat("  maximum discret (~", sprintf("%.2f", prix_formule - resultat_q2c$prix),
    ") reste entier : précision n'est pas exactitude.\n\n")


# =====================================================================
# Q2d : On copare avec la formule
# On compare le pas hebdomadaire (n = 26) au pas quotidien (n = 180).

cat("--- Q2d : comparaison avec la formule fermée ----------\n")

set.seed(123)
trajectoires_hebdo <- simuler_trajectoires(S0, r, sigma, T, 26, 100000, "exact")
lookback_hebdo <- prix_lookback_mc(trajectoires_hebdo, r, T)

set.seed(123)
trajectoires_quotidien <- simuler_trajectoires(S0, r, sigma, T, 180, 100000, "exact")
lookback_quotidien <- prix_lookback_mc(trajectoires_quotidien, r, T)

ecart_hebdo <- prix_formule - lookback_hebdo$prix
ecart_quotidien <- prix_formule - lookback_quotidien$prix

cat("N = 100000, schéma exact\n")
cat(sprintf("%-26s %10s %10s %10s\n", "", "prix", "écart", "écart %"))
cat(sprintf("%-26s %10.4f %10s %10s\n", "formule fermée", prix_formule, "-", "-"))
cat(sprintf("%-26s %10.4f %10.4f %9.2f%%\n", "MC hebdomadaire (n = 26)", lookback_hebdo$prix, ecart_hebdo, 100 * ecart_hebdo / prix_formule))
cat(sprintf("%-26s %10.4f %10.4f %9.2f%%\n", "MC quotidien (n = 180)", lookback_quotidien$prix, ecart_quotidien, 100 * ecart_quotidien / prix_formule))

# l'écart est-il compatible avec du bruit statistique ? Non : on le mesure en nombre d'erreurs-types de l'estimateur.
erreur_type_hebdo <- (lookback_hebdo$borne_sup_ic - lookback_hebdo$prix) / 1.96
cat("  écart hebdomadaire en erreurs-types :",
    sprintf("%.1f", ecart_hebdo / erreur_type_hebdo), "\n")
cat("  réduction du biais 26 -> 180 : rapport",
    sprintf("%.2f", ecart_hebdo / ecart_quotidien),
    " (attendu sqrt(180/26) =", sprintf("%.2f", sqrt(180/26)), ")\n\n")

# COMMENTAIRE.
# L'écart à la formule fermée n'est pas une erreur statistique
# il vaut plusieurs dizaines d'erreurs-types, et les intervalles de confiance excluent très largement 20 
# C'est donc le fameux biais d'estimation du maximum discret, déjà tracé sur figures/q2b_biais_loglog.png.
# Il se réduit de 3.25 (16.3 %) à 1.32 (6.6 %) en passant du pas hebdomadaire au pas quotidien, soit un rapport de 2.47, 
# à comparer au rapport de pas sqrt(180/26) = 2.63 >> c'est bien le taux O(sqrt(dt)).
# Même au pas quotidien, il reste 6.6 % d'écart : le maximum discret est un biais qui doit être dur de supprimer.


# =====================================================================
# Q2e : comparaison avec le put européen

cat("--- Q2e : lookback put contre put européen ------------\n")

set.seed(123) # pour que l'aléa soit comparable
trajectoires_q2e <- simuler_trajectoires(S0, r, sigma, T, 26, 100000, "exact")

put_q2e <- prix_put_mc(trajectoires_q2e, r, K, T)
lookback_q2e <- prix_lookback_mc(trajectoires_q2e, r, T)

# proportions de trajectoires à payoff nul
S_T_q2e <- trajectoires_q2e[, 27]
maximums_q2e <- rep(0, 100000)
for (i in 1:100000) {
  maximums_q2e[i] <- max(trajectoires_q2e[i, ])
}
proportion_nulle_put <- mean(pmax(K - S_T_q2e, 0) == 0)
proportion_nulle_lb <- mean((maximums_q2e - S_T_q2e) == 0)

# valeur théorique de la proportion de payoffs nuls du put
# le payoff est nul quand S_T >= K, ce qui arrive avec probabilité pnorm(d2) sous la mesure risque-neutre. C'est le même d2 que dans prix_put_bs.
d2 <- (log(S0 / K) + (r - sigma^2 / 2) * T) / (sigma * sqrt(T))
proportion_theorique_put <- pnorm(d2)

cat("mêmes trajectoires : N = 100000, n = 26, schéma exact\n")
cat("  prix du put européen  :", sprintf("%.4f", put_q2e$prix), "\n")
cat("  prix du lookback put  :", sprintf("%.4f", lookback_q2e$prix), "\n")
cat("  rapport (Monte-Carlo) :", sprintf("%.4f", lookback_q2e$prix / put_q2e$prix), "\n")
cat("  rapport (formules fermées) :",
    sprintf("%.4f", prix_formule / prix_put_bs(S0, r, sigma, K, T)), "\n")
cat("  proportion de payoff nul, put      :",
    sprintf("%.4f", proportion_nulle_put),
    " (théorie pnorm(d2) =", sprintf("%.4f", proportion_theorique_put), ")\n")
cat("  proportion de payoff nul, lookback :",
    sprintf("%.4f", proportion_nulle_lb), "\n\n")

cat("=====================================================\n")
cat(" FIN Q2\n")
cat("=====================================================\n")
