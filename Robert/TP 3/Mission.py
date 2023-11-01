import numpy as np
import math as mt
from scipy.optimize import fsolve

try:
    from Atmosphere import *
    from Forces import *
    from Vitesses import *
    from Acceleration import *
    from Interval import *
    from Calcul_SAR import *
except:
    pass

# Fonction pour la montee ou la descente
def Montee_ou_Descente(Hpi, Hpf, VKCAS, Vwind_kts, mach, Wi_lbs, isa, 
            choix, Choice_regime, Choice_no_Engines, step, valeur_ROC = 100):
    
    ## Montée ou descente en mission (Partie A) ##
    
    # Crée le tableau Hp de telle sorte que les incréments soient égaux à l'étape
    if choix == "Descente":
        new_step = int((Hpi - Hpf) / step) + 1
        Hp = np.arange(Hpi, Hpf - 1, -step)
            
    elif choix == "Montee":
        new_step = int((Hpf - Hpi) / step) + 1
        Hp = np.linspace(Hpi, Hpf, new_step)

    # Ajoute Hpi ou Hpf au tableau, au cas où arange ne le prendrait pas en compte
    
    # Pour trouver Hp_transition, nous devons trouver Hp de telle sorte que la valeur calculée de Mach avec conditions_vitesse
    # moins la valeur du Mach donnée soit égale à 0 (0 de la fonction)
    def difference_mach(Hp):
        atm = conditions_atm(Hp, isa, "ISA")
        mach_calculé = conditions_vitesse(VKCAS, "Vc", Wi_lbs, atm[5], atm[0], atm[3], atm[2], atm[4], atm[6])[1]
        return mach_calculé - mach 
    
    Hp_transition = fsolve(lambda Hp: difference_mach(Hp), 5000)
    
    
    # Les calculs d'itération commencent ici
    # Initialise les valeurs
    W_lbs = [Wi_lbs, Wi_lbs]
    delta_t = 0
    delta_d = 0
    delta_carburant = 0
    
    # Valeurs d'accélération
    for i in range(len(Hp)-1):
        W1 = W_lbs[-1]
        
        if Hp[i+1] == 10000:
            acc = Acceleration(Hp[i+1], VKCAS, Vwind_kts, mach, W1, isa, Choice_regime, Choice_no_Engines)
            
        delta_W = W_lbs[-1] - W_lbs[-2]
        W2i = W_lbs[i] + delta_W
        W_moyen = (W2i + W1)/2
        
        j = 1

        while j < 100:
            [W2, dt, dd, df, taux_montee] = Interval(Hp[i], Hp[i+1], Hp_transition[0], VKCAS, Vwind_kts, mach, W1, W_moyen, isa, Choice_regime, Choice_no_Engines)
            delta_W2 = W2 - W2i
            W2i = W2
            j += 1
        
        if taux_montee < valeur_ROC and choix == "Montee":
            Hp = Hp[:i]
            break
        
        delta_t += dt
        delta_d += dd
        delta_carburant += df
        
        W_lbs.append(W2)
    
    # Ajoute les valeurs d'accélération aux valeurs d'intervalle
    if Hpf > Hpi and Hpf != 10000:
        delta_t = delta_t - abs(acc[2])
        delta_d = delta_d - abs(acc[3])
        delta_carburant = delta_carburant + abs(acc[4])
        
    elif Hpf < Hpi:
        delta_t = delta_t + abs(acc[2])
        delta_d = delta_d + abs(acc[3])
        delta_carburant = delta_carburant + abs(acc[4])
        
    elif Hpf == 10000:
        delta_t = delta_t + abs(acc[2])
        delta_d = delta_d + abs(acc[3])
        delta_carburant = delta_carburant + abs(acc[4])
        
    else:
        delta_t = delta_t + abs(acc[2])
        delta_d = delta_d + abs(acc[3])
        delta_carburant = delta_carburant + abs(acc[4])
    
    return Hp, taux_montee, delta_t, delta_d, delta_carburant, Hp_transition[0], acc[0], abs(acc[1]), abs(acc[2]), abs(acc[3]), abs(acc[4]), acc[5]

## Croisière en mission (Partie B) ##
def Croisiere(Hp, isa, Vwind_kts, W_lbs, Vitesse):
     
    # Donnees pour Flap = 0 du fichier Données avion
    V_MO = 330 # kts
    M_MO = 0.85
    CDp = 0.0206
    K = 0.0364
    S = 520 # ft^2
    
    # Donnees atmospheriques 
    
    valeurs_atm = conditions_atm(Hp, isa, "ISA")
    
    # Calcul de V_MD
    CL_MD = np.sqrt(CDp/K)
    V_MD = np.sqrt((2*W_lbs)/(valeurs_atm[-1]*CL_MD*S))/1.688
    
    # Calcul de V_LRC
    CL_LRC = np.sqrt(CDp/(3*K))
    V_LRC = np.sqrt((2*W_lbs)/(valeurs_atm[-1]*CL_LRC*S))/1.688
    
    # Vitesse en fonction du choix utilisateur
    if Vitesse == 'Vmd':
        VKTAS = V_MD
    elif Vitesse == 'LRC': 
        Vitesses = np.linspace(50,500,1000)
        Valeurs_SAR = []
        
        # Calcul des valeurs SAR
        for i in range(len(Vitesses)):
            SARs = SAR(Hp, isa, W_lbs, Vitesses[i])
            Valeurs_SAR.append(SARs)
            
        # Calcul de LRC
        LRC = max(Valeurs_SAR)*0.99 # Page 22 du Module 6

        # Afin de calculer VKTAS, V_LRC doit etre calcule en fonction de SAR
        # Ainsi, on cree une fonction avec la fonction Calcul_SAR et on utilise
        # la fonction fsolve
        
        def fonction_SAR(Vitesses):
            Valeurs_SAR = SAR(Hp, isa, W_lbs, Vitesses)
            return Valeurs_SAR - LRC 
        
        V_LRC = fsolve(lambda Vitesses: fonction_SAR(Vitesses), 1.1*Vitesses[Valeurs_SAR.index(max(Valeurs_SAR))])[0]
        VKTAS = V_LRC
        
    # Verification que la vitesse de croisiere ne depasse pas V_MO/V_MO et
    # Verification que la vitesse n'est pas inferieure a V_MD
    
    valeurs_atm = conditions_atm(Hp, isa, "ISA")
    
    forces = conditions_forces(VKTAS, "V", "MCR", "CAS_constant",
                      "AEO", W_lbs, isa, Hp, valeurs_atm[0], valeurs_atm[2], 
                      valeurs_atm[3], valeurs_atm[4], valeurs_atm[5], valeurs_atm[6], 0, "UP", 
                      1, "Nz", 0)
    
    # forces[-3] : VITESSE CALIBREE
    # forces[-2] : NOMBRE DE MACH
    # forces[-1] : VITESSE VRAIE
    
    if forces[-3] > V_MO or forces[-2] > M_MO:
        print("La vitesse de croisière dépasse V_MO/M_MO")
        
    if forces[-1] < 0.99*V_MD:
        print("La vitesse de croisière est inférieure à V_MD")
    
    T_req = W_lbs/(forces[14]/forces[15])
    SFC = 0.58 + (0.035 * Hp/10000)
    
    Wf = SFC*T_req

    Valeur_SAR = SAR(Hp, isa, W_lbs, forces[-1])
    Valeur_SR = (forces[-1] + Vwind_kts)/Wf

    Valeur_VG = Vwind_kts + forces[-1]

    return Valeur_SAR, Valeur_SR, Valeur_VG