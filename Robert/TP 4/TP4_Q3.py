import numpy as np
import matplotlib.pyplot as plt
try:
    from Mission import *
except:
    pass

################# MONTEE #################
# Calcul du poids a 1500 ft (avant la montée)
ZFW = 31500 + (20*225) 
RW = ZFW + 15000 # Poids à la rampe
TOW = RW - 200 # Poids au décollage
ETO = TOW - 250 # Poids à la fin du ségment de décollage

Hpi = 1500
Hpf = 41000

VWIND = 0
Mach = 0.78
Wi = ETO
ISA = 0
regime = "MCL"
no_moteurs = "AEO"

values_montee = Montee_ou_Descente(Hpi, Hpf, 275, VWIND, Mach, Wi, 
                            ISA, "Montee", regime, no_moteurs, 100, 300)


distance_montee = values_montee[1]

TOC = ETO - values_montee[2] # Poids au sommet de la montée

################# DESCENTE #################
# Calcul du poids au début de l'approche et de l'atterissage
LW = ZFW + 2000
BAL = LW + 200

Hpi = values_montee[10][-1] # Altitude de croisière déterminée durant la phase de montée
Hpf = 1500
    
VWIND = 0
Mach = 0.78
Wi = ETO
ISA = 0
regime = "Idle"
no_moteurs = "AEO"

# Calcul itératif pour déterminer TOD (poids au début de la descente) et les valeurs pour la descente

# TOD_var = np.linspace(BAL, TOC, 1000)
# valeurs_descente = []
# for i in range(len(TOD_var)):
#     valeurs_descente = Montee_ou_Descente(Hpi, Hpf, 275, VWIND, Mach, TOD_var[i], 
#                                 ISA, "Descente", regime, no_moteurs, 1000, 300)
    
#     carburant_descente = valeurs_descente[2]
    
#     if abs(TOD_var[i] - BAL - carburant_descente) < 10:
#         TOD = TOD_var[i]
#         distance_descente = valeurs_descente[1]
#         carburant_descente = valeurs_descente[2]

# Valeurs calculées avec la boucle ci-haut, mais temps de calcul très long, donc boucle for en comment
# Seulement les résultats qui nous intéressent sont gardés

TOD = 38306.89741098786
distance_descente = 86.71139976385341
carburant_descente = 113.12234557417227

################# CROISIERE #################
# Il est nécessaire de définir un delta_t = 1 secondes,
# temps dans lequel l'avion franchit une distance donnée
# ceci est dans le but de calculer une delta_d à partir de Vg
# à chaque itération. Pour le critère d'arrêt, on vérifie que
# le poids ne se situe plus entre TOC et TOD avec comme variation
# le carburant consommé durant la croisière.

delta_t = 1/3600
Valeurs_SAR = 0
carburant_croisiere = 0
W_croisiere = 0.001
while abs(TOC - TOD - carburant_croisiere) > 1:
    valeurs_croisiere = Croisiere(values_montee[10][-1], ISA, VWIND, W_croisiere, "Mach", Mach)
    carburant_croisiere += valeurs_croisiere[-1]*delta_t
    W_croisiere = TOC - carburant_croisiere
    Valeurs_SAR += valeurs_croisiere[1]*delta_t
    
# Calcul de la distance totale et du carburant consommée

distance_totale = distance_montee + distance_descente + Valeurs_SAR

# Pour le calcul du carburant consommé, on assume que l'avion ait
# consommé tout le carburant pendant le vol

OWE = 31500 # Poids à vide
carburant_consommée = RW - (OWE + (20*225) + (2000-100))

# Affichage des résultats

print("\nRésultats :")
print("-" * 60)
print("Distance de montée (nm) : " + str("{:.5g}".format(distance_montee)))
print("Altitude de croisière (ft) : " + str("{:.5g}".format(values_montee[10][-1])))
print("Distance de croisière (nm) : " + str("{:.5g}".format(Valeurs_SAR)))
print("Distance de descente (nm) : " + str("{:.5g}".format(distance_descente)))
print("Distance totale (nm) : " + str("{:.5g}".format(distance_totale)))
print("Carburant total consommé (lb) : " + str("{:.5g}".format(carburant_consommée)))
print("-" * 60)
print("\n")