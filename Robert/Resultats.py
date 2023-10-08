import numpy as np
try:
    from Atmosphere import *
    from Forces import *
    from Transformations import *
except:
    pass

#Atmospheric values
#%% ------------------------------------------------------------------ #%%

#Prompts the user for values_atm
hauteur_pression = input("Entrer la valeur de l'hauteur pression (ft) : ")  
choix_T_ou_ISA = input("Choisir le paramètre à utiliser pour les calculs (T ou ISA): ")

while choix_T_ou_ISA not in ('T', 'ISA'):
    choix_T_ou_ISA = input("Veuillez choisir entre les valeurs suivantes : T ou ISA : ")

if choix_T_ou_ISA == "T": 
    
    temp = input("Entrer la valeur de la température (C) : ")
    while temp.isnumeric() == False:
        temp = input("Entrer la valeur de la température (C) : ")
        
    values_atm = conditions_atm(float(hauteur_pression), float(temp), choix_T_ou_ISA)

    print("\nRésultats atmosphériques :")
    print("-" * 60)
    print("Température (K) : " + str("{:.4g}".format(values_atm[0])))
    print("ISA (C) : " + str("{:.4g}".format(values_atm[1])))
    print("Theta : " + str("{:.4g}".format(values_atm[2])))
    print("Delta : " + str("{:.4g}".format(values_atm[3])))
    print("Sigma : " + str("{:.4g}".format(values_atm[4])))
    print("Pression (lb/ft^2) : " + str("{:.4g}".format(values_atm[5])))
    print("Rho (slugs/ft^3) : " + str("{:.4g}".format(values_atm[6]))) 
    print("-" * 60)
    print("\n")

    isa = str(values_atm[1]) # pour lorsqu'on donne une temperature T au lieu de ISA
    
elif choix_T_ou_ISA == "ISA":
    
    isa = input("Entrer la valeur de la déviation par rapport à delta ISA (C) : ")
    while isa.isnumeric() == False:
        isa = input("Entrer la valeur de la déviation par rapport à delta ISA (C) : ")
        
    values_atm = conditions_atm(float(hauteur_pression), float(isa), choix_T_ou_ISA)
    
    print("\nRésultats atmosphériques :")
    print("-" * 60)
    print("Température (K) : " + str("{:.4g}".format(values_atm[0])))
    print("Température (C) : " + str("{:.4g}".format(values_atm[1])))
    print("Theta : " + str("{:.4g}".format(values_atm[2])))
    print("Delta : " + str("{:.4g}".format(values_atm[3])))
    print("Sigma : " + str("{:.4g}".format(values_atm[4])))
    print("Pression (lb/ft^2) : " + str("{:.4g}".format(values_atm[5])))
    print("Rho (slugs/ft^3) : " + str("{:.4g}".format(values_atm[6])))
    print("-" * 60)
    print("\n")
    
#Speed values
#%% ------------------------------------------------------------------ #%%

#Prompts the user for values_forces

choix_vitesse = input("Entrer le choix de vitesse (Vc, Ve, V, Mach ou Vsr) : ")
while choix_vitesse not in ("Vc", "Ve", "V", "Mach", "Vsr"):
    choix_vitesse = input("Entrer le choix de vitesse (Vc, Ve, V, Mach ou Vsr) : ")

vitesse = input("Entrer la vitesse basée sur votre choix précédent (kts ou sans unité) : ")

regime = input("Entrer le choix du régime (MCR, MTO, GA, Idle, MCT, MCL ou Poussee_totale) : ")
while regime not in ("MCR", "MTO", "GA", "Idle", "MCT", "MCL", "Poussee_totale"):
    regime = input("Entrer le choix du régime (MCR, MTO, GA, Idle, MCT, MCL ou Poussee_totale) : ")

no_moteurs = input("Entrer le nombre de moteurs (AEO, OEI) - Mettre n'importe lequel si vous avez choisi Poussee_totale : ")
while no_moteurs not in ("AEO", "OEI"):
    no_moteurs = input("Entrer le nombre de moteurs (AEO, OEI) - Mettre n'importe lequel si vous avez choisi Poussee_totale : ")
   
regime_climb = input("Entrer le choix du régime pour la montée (Mach_constant ou CAS_constant) : ")
while regime_climb not in ("Mach_constant", "CAS_constant"):
    regime_climb = input("Entrer le choix du régime pour la montée (Mach_constant ou CAS_constant) : ")

poids = input("Entrer le poids de l'avion (lbs) : ")

flaps = input("Entrer l'angle des flaps (0, 20, 45) : ")
while flaps not in ["0", "20", "45"]:
    flaps = input("Entrer l'angle des flaps (0, 20, 45) : ")

lg = input("Entrer la position des flaps (UP, DOWN) : ")
while lg not in ("UP", "DOWN"):
    lg = input("Entrer la position des flaps (UP, DOWN) : ")
    
cg = input("Entrer la position du centre de gravité (%MAC) : ")

choix_nz_ou_f = input("Choisir entre entrer la valeur de Nz ou f : ")

while choix_nz_ou_f not in ("Nz", "f"):
    choix_nz_ou_f = input("Choisir entre entrer la valeur de Nz ou f : ")
    
if choix_nz_ou_f == "Nz":
    nz = input("Entrer la valeur de Nz (g) : ")
    
elif choix_nz_ou_f == "f":
    nz = input("Entrer la valeur de f (deg) : ")
  
values_forces = conditions_forces(float(vitesse), choix_vitesse, 
                                  regime, regime_climb, no_moteurs, float(poids), 
                                  float(isa), float(hauteur_pression), values_atm[0], 
                                  values_atm[2], values_atm[3], values_atm[4], 
                                  values_atm[5], values_atm[6], float(flaps), lg, 
                                  float(nz), choix_nz_ou_f, float(cg))

print("\nRésultats forces :")
print("-" * 60)
print("CD : " + str("{:.5g}".format(values_forces[0])))
print("CL : " + str("{:.5g}".format(values_forces[1])))
print("Portance (lbs) : " + str("{:.5g}".format(values_forces[2])))
print("Finesse (L/D) : " + str("{:.5g}".format(values_forces[3])))
print("CDp : " + str("{:.5g}".format(values_forces[4])))
print("Dp (lbs) : " + str("{:.5g}".format(values_forces[5])))
print("CDi : " + str("{:.5g}".format(values_forces[6])))
print("Di (lbs) : " + str("{:.5g}".format(values_forces[7])))
print("CDcomp : " + str("{:.5g}".format(values_forces[8])))
print("Dcomp (lbs) : " + str("{:.5g}".format(values_forces[9])))
print("CDcntl : " + str("{:.5g}".format(values_forces[10])))
print("Dcntl (lbs) : " + str("{:.5g}".format(values_forces[11])))
print("CDwm : " + str("{:.5g}".format(values_forces[12])))
print("Dwm (lbs) : " + str("{:.5g}".format(values_forces[13])))
print("Poussée totale (lbs) : " + str("{:.5g}".format(values_forces[14])))
print("Trainée (lbs) : " + str("{:.5g}".format(values_forces[15])))
print("Angle d'attaque (deg) : " + str("{:.5g}".format(values_forces[16])))
print("Nz sw : " + str("{:.5g}".format(values_forces[17])))
print("Phi sw (deg) : " + str("{:.5g}".format(values_forces[18])))
print("Nz buffet : " + str("{:.5g}".format(values_forces[19])))
print("Gamma montée : " + str("{:.5g}".format(values_forces[20])))
print("Taux de montée : " + str("{:.5g}".format(values_forces[21])))
print("Taux de montée pression : " + str("{:.5g}".format(values_forces[22])))
print("Facteur d'accéleration : " + str("{:.5g}".format(values_forces[23])))
print("Accélération selon l’axe de la trajectoire de vol : " + str("{:.5g}".format(values_forces[24])))
print("-" * 60)
print("\n")