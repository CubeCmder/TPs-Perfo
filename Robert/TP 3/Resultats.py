import numpy as np
try:
    from Mission import *
    from Calcul_SAR import *
except:
    pass

# Hpi = input("Entrer l'hauteur initiale (ft) : ")
# Hpf = input("Entrer l'hauteur finale (ft) : ")

# print("-------------------------------------")

# if Hpf < Hpi:
#     print("Vous êtes en DESCENTE")
#     choix = "Descente"
# else:
#     print("Vous êtes en MONTÉE")
#     choix = "Montee"
    
# print("-------------------------------------")

# VKCAS = input("Entrer la valeur du VKCAS (kts) : ")
# VWIND = input("Entrer la valeur du VWIND (kts) : ")
# Mach = input("Entrer la valeur du nombre de Mach (-) : ")
# Wi = input("Entrer la valeur de la masse initiale (lbs) : ")
# ISA = input("Entrer la valeur de la déviation par rapport à ISA (C) : ")

# regime = input("Entrer le choix du régime (MCR, MTO, GA, Idle, MCT, MCL ou Poussee_totale) : ")
# while regime not in ("MCR", "MTO", "GA", "Idle", "MCT", "MCL", "Poussee_totale"):
#     regime = input("Entrer le choix du régime (MCR, MTO, GA, Idle, MCT, MCL ou Poussee_totale) : ")

# no_moteurs = input("Entrer le nombre de moteurs (AEO, OEI) - Mettre n'importe lequel si vous avez choisi Poussee_totale : ")
# while no_moteurs not in ("AEO", "OEI"):
#     no_moteurs = input("Entrer le nombre de moteurs (AEO, OEI) - Mettre n'importe lequel si vous avez choisi Poussee_totale : ")

# values = Montee_ou_Descente(float(Hpi), float(Hpf), float(VKCAS), 
#                               float(VWIND), float(Mach), float(Wi), 
#                               float(ISA), choix, regime, no_moteurs, 10)

# print("\nRésultats montée et descente :")
# print("-" * 60)
# print("Temps total (s) : " + str("{:.5g}".format(values[0])))
# print("Distance totale (NM) : " + str("{:.5g}".format(values[1])))
# print("Carburant total consommé (lbs) : " + str("{:.5g}".format(values[2])))
# print("Altitude transition (ft) : " + str("{:.5g}".format(values[3])))
# print("Accélération : (ft/s^2) : " + str("{:.5g}".format(values[4])))
# print("Vitesse vraie moyenne(ft/s) : " + str("{:.5g}".format(values[5])))
# print("Temps d'accélération (s) : " + str("{:.5g}".format(values[6])))
# print("Distance d'accélération (NM) : " + str("{:.5g}".format(values[7])))
# print("Carburant consommé pendant l'accélération' (lbs) : " + str("{:.5g}".format(values[8])))
# print("Masse initiale au début de l'accélération (lbs) : " + str("{:.5g}".format(values[9])))
# print("-" * 60)
# print("\n")

# values = Montee_ou_Descente(1500, 41000, 275, 
#                               20, 0.74, 52800, 
#                               10, "Montee", "MCL", "AEO", 100)

# values = Montee_ou_Descente(38000, 1500, 275, 
#                               -20, 0.74, 47000, 
#                               -10, "Descente", "Idle", "AEO", 100)

# values = Montee_ou_Descente(9000, 10000, 275, 
#                               20, 0.74, 50000, 
#                               10, "Montee", "MCL", "AEO", 100)

# SAR_func = SAR(15000, 0, 40000, 50)

# # Verification de l'altitude de croisiere
# Hp = input("Veuillez choisir une altitude qui se situe entre 2000 ft et 41000 ft : ")
# while float(Hp) < 2000 or float(Hp) > 41000:
#     Hp = input("Veuillez choisir une altitude qui se situe entre 2000 ft et 41000 ft : ")

a = Croisiere(15000, 0, 0, 40000, "LRC")