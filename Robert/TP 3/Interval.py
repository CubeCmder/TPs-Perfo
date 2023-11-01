import numpy as np
import math as mt
from scipy.optimize import fsolve

try:
    from Atmosphere import *
    from Forces import *
    from Vitesses import *
except:
    pass

def Interval(Hp1, Hp2, Hp_transition, VKCAS, Vwind_kts, mach, W1, W_avg, isa, Choix_regime, Choix_no_Engines):
    # Hp moyen
    Hp_moyen = (Hp1 + Hp2) / 2
    
    if Hp_moyen < 10000:
        Vitesse = 250
        Choix_vitesse = "Vc"
        Choix_regime_montee = "CAS_constant"
        
    elif Hp_moyen > Hp_transition:
        Vitesse = mach
        Choix_vitesse = "Mach"
        Choix_regime_montee = "Mach_constant"
        
    else:
        Vitesse = VKCAS
        Choix_vitesse = "Vc"
        Choix_regime_montee = "CAS_constant"
    
    # Température avec Hp_moyen calculé
    T_K = conditions_atm(Hp_moyen, isa, "ISA")[0]
    T_K_std = conditions_atm(Hp_moyen, 0, "ISA")[0]

    # Variation de Hp
    delta_Hp = Hp2 - Hp1
    delta_Hp_geo = delta_Hp * (T_K / T_K_std)
    
    valeurs_atm = conditions_atm(Hp_moyen, isa, "ISA")
    
    forces = conditions_forces(Vitesse, Choix_vitesse, Choix_regime, Choix_regime_montee,
                      Choix_no_Engines, W_avg, isa, Hp_moyen, valeurs_atm[0], valeurs_atm[2], 
                      valeurs_atm[3], valeurs_atm[4], valeurs_atm[5], valeurs_atm[6], 0, "UP", 
                      1, "Nz", 0)
    
    # Calcul de V, T, D, taux_montee et FA
    V = forces[-1] # kts
    T = forces[14] # lbs
    D = forces[15] # lbs
    taux_montee = forces[21]
    FA = forces[23] # -
    
    # Donné dans la description du TP
    if T < 0:
        T = 1200 # lbs

    # Delta t
    delta_t = (delta_Hp_geo/taux_montee)*60 # s
    
    # Delta d
    delta_d = (V + Vwind_kts)*(delta_t/3600) # NM
    
    # SFC
    SFC = 0.58 + (0.035 * Hp_moyen/10000)
    
    # Calcul de Wf_moyen
    Wf_moyen = SFC * T
    
    # delta_fuel
    delta_fuel = Wf_moyen * (delta_t/3600)
    
    # W2
    W2 = W1 - delta_fuel
    
    return W2, delta_t, delta_d, delta_fuel, taux_montee
