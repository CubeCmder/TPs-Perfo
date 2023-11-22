import numpy as np
import math as mt
from scipy.interpolate import interp1d
import sys

try:
    from Atmosphere import *
    from Forces import *
    from Vitesses import *
    from Decollage import *
    from Transformations import *
except:
    pass

def Longueur_Piste(V1VR, W, Hp, isa, S = 520):
    
    # Données
    RAD = 0 # distance d'alignement sur la piste
    flap = 20 # degres
    LG = "UP"
    CG = 0.25
    Choice_Nzf = "Nz"
    f = 1
    Choice_regime_climb = "CAS_constant"
    
    # Calcul des conditions atmospheriques
    atm = conditions_atm(Hp, isa, "ISA")
    theta = atm[2]
    
    # Calcul des valeurs de décollage
    valeurs_decollage = Decollage(W, Hp, isa)
    Vr = valeurs_decollage[2] # fts
    
    # Calcul de V1
    V1 = V1VR * Vr # fts
    V1mcg = 95 # kts # fichier avion
    V1mcg = conditions_vitesse(V1mcg, "Vc", W, atm[5], atm[0], atm[3], 
                       atm[2], atm[4], atm[6], LG, flap) # kts
    V1mcg = trans_kts_to_fts(V1mcg) # fts
    
    # Verification que V1/VR spécifié en entrée résulte en une valeur de V1 qui est supérieure à V1MCG
    if V1mcg > V1:
        print("\n")
        print("Attention ! La vitesse V1 est trop faible. Celle-ci sera augmentée pour respecter le ratio V1/Vr fourni.")
        Vr = V1mcg/V1VR # fts
        V1 = V1mcg # fts
        
    # Decollage AEO
    Choice_no_Engines = "AEO"
    
    ############################### Portion V0 = 0 a V = V1 ###############################
    
    # MU
    rollmu = 0.020 # fichier avion
    
    # Vitesse initiale et finale
    V_initial = 0 # fts
    V_final = V1 # fts
    
    # Calcul de vitesse quadratique moyenne (RMS)
    Vrms = (np.sqrt(2)/2)*(V_initial**2 + V_final**2)**(1/2) # fts
    
    Choice_regime = "MTO"
    # Calcul de la poussee, q, CL et CD
    forces = conditions_forces(trans_fts_to_kts(Vrms), "V", Choice_regime, Choice_regime_climb,
                          Choice_no_Engines, W, isa, Hp, atm[0], atm[2], atm[3], 
                          atm[4], atm[5], atm[6], flap, LG, f, Choice_Nzf, CG)
    T_lbs = forces[14]
    q = forces[-1]
    CL = 0.8290 # fichier avion - flap = 20
    CD = 0.0750 # fichier avion - flap = 20
    
    # Calcul de l'acceleration a la vitesse RMS
    acc_moy_RMS = (cst.g/W)*((T_lbs - rollmu*W) - (CD - rollmu*CL)*q*S) # Module 10 p. 29
    
    # Calcul de delta_t et delta_s
    delta_t_V0_V1_AEO = (V_final - V_initial)/acc_moy_RMS # s
    delta_s_V0_V1_AEO = delta_t_V0_V1_AEO*(V_initial + ((V_final - V_initial)/2)) # m
    
    # ############################### Portion V = V1 a V = Vr ###############################
    
    # MU
    rollmu = 0.020 # fichier avion
    
    # Vitesse initiale et finale
    V_initial = V1 # fts
    V_final = Vr # fts
    
    # Calcul de vitesse quadratique moyenne (RMS)
    Vrms = (np.sqrt(2)/2)*(V_initial**2 + V_final**2)**(1/2) # fts

    Choice_regime = "MTO"
    # Calcul de la poussee, q, CL et CD
    forces = conditions_forces(trans_fts_to_kts(Vrms), "V", Choice_regime, Choice_regime_climb,
                          Choice_no_Engines, W, isa, Hp, atm[0], atm[2], atm[3], 
                          atm[4], atm[5], atm[6], flap, LG, f, Choice_Nzf, CG)
    T_lbs = forces[14]
    q = forces[-1]
    CL = 0.8290 # fichier avion - flap = 20
    CD = 0.0750 # fichier avion - flap = 20
    
    # Calcul de l'acceleration a la vitesse RMS
    acc_moy_RMS = (cst.g/W)*((T_lbs - rollmu*W) - (CD - rollmu*CL)*q*S) # Module 10 p. 29
    
    # Calcul de delta_t et delta_s
    delta_t_V1_Vr_AEO = (V_final - V_initial)/acc_moy_RMS # s
    delta_s_V1_Vr_AEO = delta_t_V1_Vr_AEO*(V_initial + ((V_final - V_initial)/2)) # m

    
    ############################### Portion V = V1 a V = V0 ###############################
    
    # MU
    mudry = 0.400 # fichier avion
    
    # Vitesse initiale et finale
    V_initial = V1 # fts
    V_final = 0 # fts
    
    # Calcul de vitesse quadratique moyenne (RMS)
    Vrms = (np.sqrt(2)/2)*(V_initial**2 + V_final**2)**(1/2) # fts
    
    Choice_regime = "Idle"
    # Calcul de la poussee, q, CL et CD
    forces = conditions_forces(trans_fts_to_kts(Vrms), "V", Choice_regime, Choice_regime_climb,
                          Choice_no_Engines, W, isa, Hp, atm[0], atm[2], atm[3], 
                          atm[4], atm[5], atm[6], flap, LG, f, Choice_Nzf, CG)
    T_lbs = forces[14]
    q = forces[-1]
    CL = 0.2090 # fichier avion - flap = 20
    CD = 0.1171 # fichier avion - flap = 20
    
    # Calcul de l'acceleration a la vitesse RMS
    acc_moy_RMS = (cst.g/W)*((T_lbs - mudry*W) - (CD - mudry*CL)*q*S) # Module 10 p. 29
    
    # Calcul de delta_t et delta_s
    delta_t_Vr_V0_AEO = (V_final - V_initial)/acc_moy_RMS # s
    delta_s_Vr_V0_AEO = delta_t_Vr_V0_AEO*(V_initial + ((V_final - V_initial)/2)) # m
    
    # Decollage OEI
    Choice_no_Engines = "OEI"
    
    ############################### Portion V = V1 a V = Vr ###############################
    
    # MU
    rollmu = 0.020 # fichier avion
    
    # Vitesse initiale et finale
    V_initial = V1 # fts
    V_final = Vr # fts
    
    # Calcul de vitesse quadratique moyenne (RMS)
    Vrms = (np.sqrt(2)/2)*(V_initial**2 + V_final**2)**(1/2) # fts
    Vrms = conditions_vitesse(trans_fts_to_kts(Vrms), "V", W, atm[5], atm[0], atm[3], 
                    atm[2], atm[4], atm[6], LG, flap) # kts
    
    Choice_regime = "MTO"
    # Calcul de la poussee, q, CL et CD
    forces = conditions_forces(Vrms, "V", Choice_regime, Choice_regime_climb,
                          Choice_no_Engines, W, isa, Hp, atm[0], atm[2], atm[3], 
                          atm[4], atm[5], atm[6], flap, LG, f, Choice_Nzf, CG)
    T_lbs = forces[14]
    q = forces[-1]
    CL = 0.8290 # fichier avion - flap = 20
    CD = 0.0750 # fichier avion - flap = 20
    
    # Calcul de l'acceleration a la vitesse RMS
    acc_moy_RMS = (cst.g/W)*((T_lbs - rollmu*W) - (CD - rollmu*CL)*q*S) # Module 10 p. 29
    
    # Calcul de delta_t et delta_s
    delta_t_V1_Vr_OEI = (V_final - V_initial)/acc_moy_RMS # s
    delta_s_V1_Vr_OEI = delta_t_V1_Vr_OEI*(V_initial + ((V_final - V_initial)/2)) # m
    
    ############################### Portion V = V1 a V = V0 ###############################
    
    # MU
    mudry = 0.400 # fichier avion
    
    # Vitesse initiale et finale
    V_initial = V1 # fts
    V_final = 0 # fts
    
    # Calcul de vitesse quadratique moyenne (RMS)
    Vrms = (np.sqrt(2)/2)*(V_initial**2 + V_final**2)**(1/2) # fts
    Vrms = conditions_vitesse(trans_fts_to_kts(Vrms), "V", W, atm[5], atm[0], atm[3], 
                    atm[2], atm[4], atm[6], LG, flap) # kts
    
    Choice_regime = "Idle"
    # Calcul de la poussee, q, CL et CD
    forces = conditions_forces(Vrms, "V", Choice_regime, Choice_regime_climb,
                          Choice_no_Engines, W, isa, Hp, atm[0], atm[2], atm[3], 
                          atm[4], atm[5], atm[6], flap, LG, f, Choice_Nzf, CG)
    T_lbs = forces[14]
    q = forces[-1]
    CL = 0.2090 # fichier avion - flap = 20
    CD = 0.1171 # fichier avion - flap = 20
    
    # Calcul de l'acceleration a la vitesse RMS
    acc_moy_RMS = (cst.g/W)*((T_lbs - mudry*W) - (CD - mudry*CL)*q*S) # Module 10 p. 29
    
    # Calcul de delta_t et delta_s
    delta_t_Vr_V0_OEI = (V_final - V_initial)/acc_moy_RMS # s
    delta_s_Vr_V0_OEI = delta_t_Vr_V0_OEI*(V_initial + ((V_final - V_initial)/2)) # m

    
    # Calcul de FTOD_AEO
    TOD_AEO = delta_s_V0_V1_AEO + delta_s_V1_Vr_AEO + valeurs_decollage[13] + valeurs_decollage[14]
    FTOD_AEO = TOD_AEO * 1.15
    
    # Calcul de TOD_OEI
    TOD_OEI = delta_s_V0_V1_AEO + delta_s_V1_Vr_OEI + valeurs_decollage[11] + valeurs_decollage[12]
    
    # Calcul de ASD_AEO
    ASD_AEO = delta_s_V0_V1_AEO + delta_s_Vr_V0_AEO + 2*V1
    
    # Calcul de la longueur minimum requise pour la longueur de la piste de décollage avec la valeur de V1 spécifiée
    Longueur_minimum = max(FTOD_AEO, TOD_AEO, ASD_AEO, TOD_OEI)
    
    # Calcul de l'energie totale absorbée par les freins pendant l’arrêt (en millions ft-lb)
    energie_totale = 0.4*(W - q*S*CL)*max(delta_s_Vr_V0_AEO, delta_s_Vr_V0_OEI)/1e6
    
    # return delta_s_V0_V1_AEO, delta_t_V0_V1_AEO
    return FTOD_AEO, TOD_OEI, ASD_AEO, Longueur_minimum, energie_totale