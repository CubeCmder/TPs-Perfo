import numpy as np
import math as mt
from scipy.interpolate import interp1d
import sys

try:
    from Atmosphere import *
    from Forces import *
    from Vitesses import *
    from Transformations import *
except:
    pass

def Decollage(W, Hp, isa, S = 520):
    
    # Donnees
    flap = 20 # degres
    LG = "UP"
    CLmax = 1.85 # fichier avion
    CG = 0.25
    Choice_Nzf = "Nz"
    f = 1
    Choice_regime = "MTO"
    Choice_regime_climb = "CAS_constant"
    
    # Calcul des conditions atmospheriques
    atm = conditions_atm(Hp, isa, "ISA")
    
    # Calcul des vitesses minimum permises par les regles de certification
    # Calcul de la vitesse vraie basee sur Vsr
    Vsr = np.sqrt((295.37*W)/(CLmax*S))
    Vsr = conditions_vitesse(Vsr, "Ve", W, atm[5], atm[0], atm[3], 
                        atm[2], atm[4], atm[6], LG, flap)
    
    # Calcul de la vitesse vraie basee sur Vmca
    Vmca = 95 # fichier avion
    Vmca = conditions_vitesse(Vmca, "Vc", W, atm[5], atm[0], atm[3], 
                       atm[2], atm[4], atm[6], LG, flap)
    
    # Calcul de la vitesse vraie basee sur Vmcl
    Vmcl = 92 # fichier avion
    Vmcl = conditions_vitesse(Vmcl, "Vc", W, atm[5], atm[0], atm[3], 
                       atm[2], atm[4], atm[6], LG, flap)
    
    # Calcul de la vitesse vraie basee sur Vmcg
    Vmcg = 90 # fichier avion
    Vmcg = conditions_vitesse(Vmcg, "Vc", W, atm[5], atm[0], atm[3], 
                       atm[2], atm[4], atm[6], LG, flap)
    
    # Calcul de la vitesse vraie basee sur V1mcg
    V1mcg = 95 # fichier avion
    V1mcg = conditions_vitesse(V1mcg, "Vc", W, atm[5], atm[0], atm[3], 
                       atm[2], atm[4], atm[6], LG, flap)
    
    # Calcul de V2 max
    V2 = max(1.13 * Vsr, 1.1 * Vmca)

    # Calcul du gradient lorsque les deux moteurs sont en fonction (AEO)
    Choice_no_Engines = "AEO"
    
    gradient_montee_AEO = conditions_forces(V2, "V", Choice_regime, Choice_regime_climb,
                          Choice_no_Engines, W, isa, Hp, atm[0], atm[2], atm[3], 
                          atm[4], atm[5], atm[6], flap, LG, f, Choice_Nzf, CG)[20]
    
    # Calcul du gradient lorsque un seul des moteurs est en fonction (OEI)
    Choice_no_Engines = "OEI"
    
    gradient_montee_OEI = conditions_forces(V2, "V", Choice_regime, Choice_regime_climb,
                          Choice_no_Engines, W, isa, Hp, atm[0], atm[2], atm[3], 
                          atm[4], atm[5], atm[6], flap, LG, f, Choice_Nzf, CG)[20]
    
    if gradient_montee_OEI < 2.4:
        print("\n")
        print("!!!! WARNING !!!!")
        print("On ne rencontre pas le gradient minimum requis de 2.4% spécifié par le FAR 25.121 (b).")
        sys.exit()
        
    # Calcul des vitesses selon les increments de vitesse du fichier avion

    increments_vitesse = [[0.000,  0.010,  0.020,  0.100,  0.120,  0.200,  0.400,  0.60], # (Climb gradient at 35 ft – θ) (based on V2, gear up and AF=0)
                          [1.80,   2.30,   2.75,   6.50,   7.35,  10.80,  19.80,  30.0], # Speed spread from rotation to lift-off (normal rotation)(ft/s)
                          [0.00,   0.00,   1.00,   9.00,  11.00,  16.70,  31.00,  45.0], # Speed spread from  lift-off to 35 ft  (ft/s)
                          [0.00,   0.00,   1.00,   8.05,   9.80,  13.00,  21.00,  29.0]] # Speed spread from  lift-off to 15 ft (ft/s)
    
    # Transformation de V2 a fts
    V2 = trans_kts_to_fts(V2)
    
    # OEI
    Vloei = -interp1d(increments_vitesse[0], increments_vitesse[2])(gradient_montee_OEI/100) + V2
    Vr = Vloei - interp1d(increments_vitesse[0], increments_vitesse[1])(gradient_montee_OEI/100)
    
    # Calcul de VRmin
    Vrmin = trans_kts_to_fts(max(1.05*Vmca, V1mcg))
    
    # Condition de Vr donnee dans l'enonce
    if Vr < Vrmin:
        Vr = Vrmin
        Vloei = Vr + interp1d(increments_vitesse[0], increments_vitesse[1])(gradient_montee_OEI/100)
        V2 = Vloei + interp1d(increments_vitesse[0], increments_vitesse[2])(gradient_montee_OEI/100)
    
    # AEO
    Vlaeo = Vr + interp1d(increments_vitesse[0], increments_vitesse[1])(gradient_montee_AEO/100)
    V35aeo = Vlaeo + interp1d(increments_vitesse[0], increments_vitesse[2])(gradient_montee_AEO/100)
    
    # Calcul du temps
    
    increments_temps = [[0.0,   0.02,  0.03,  0.05,  .065,  .075,  0.10,  0.18,   0.2,   0.4,   0.6], # (Climb gradient at 35 ft – θ) (based on V2, gear up and AF=0)
                          [2.35,   2.25,  2.19,  2.09,  2.01,  1.96,  1.83,  1.41,  1.30,  1.30,  1.30], # Time between Vr and Vlo
                          [17.00,  10.00,  6.60,  5.50,  5.10,  4.90,  4.50,  3.80,  3.80,  3.80,  3.80], # Time between Vlo and 35 ft
                          [12.00,   7.20,  4.70,  4.10,  3.92,  3.80,  3.50,  2.55,  2.55,  2.55,  2.55]] # Time between Vlo and 15 ft
    
    # OEI
    DtrvloOEI = interp1d(increments_temps[0], increments_temps[1])(gradient_montee_OEI/100)
    Dtlo35OEI = interp1d(increments_temps[0], increments_temps[2])(gradient_montee_OEI/100)
    
    # AEO
    DtrvloAEO = interp1d(increments_temps[0], increments_temps[1])(gradient_montee_AEO/100)
    Dtlo35AEO = interp1d(increments_temps[0], increments_temps[2])(gradient_montee_AEO/100)
    
    # Calcul des distances
    
    #OEI
    
    DistrvloOEI = DtrvloOEI*(Vr + Vloei)/2
    Distlo35OEI = Dtlo35OEI*(V2 + Vloei)/2

    #AEO
    
    DistrvloAEO = DtrvloAEO*(Vr + Vlaeo)/2
    Distlo35AEO = Dtlo35AEO*(V35aeo + Vlaeo)/2

    # Transformations
       
    V1min = V1mcg*1.68781
    V1max = Vr
    
    return V1min, V1max, Vr, V2, Vloei, Vlaeo, V35aeo, DtrvloOEI, Dtlo35OEI, DtrvloAEO, Dtlo35AEO, DistrvloOEI, Distlo35OEI, DistrvloAEO, Distlo35AEO