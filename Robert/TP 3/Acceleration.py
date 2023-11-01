import numpy as np
import math as mt
from scipy.optimize import fsolve

try:
    from Atmosphere import *
    from Forces import *
    from Vitesses import *
except:
    pass

def Acceleration(Hp, VKCAS, Vwind_kts, mach, W1, isa, Choix_regime, Choice_no_Engines):
    
    Choix_regime_montee = "CAS_constant"
    
    Vc_moyen = (VKCAS + 250)/2
    
    # Conditions atmosphériques
    valeurs_atm = conditions_atm(Hp, isa, "ISA") 
    
    # Calcul des forces lorsque Vc = Vc_moyen
    forces_Vc_moyen = conditions_forces(Vc_moyen, "Vc", Choix_regime, Choix_regime_montee,
                      Choice_no_Engines, W1, isa, Hp, valeurs_atm[0], valeurs_atm[2], 
                      valeurs_atm[3], valeurs_atm[4], valeurs_atm[5], valeurs_atm[6], 0, "UP", 
                      1, "Nz", 0)
    
    # Calcul de la vraie vitesse lorsque VKCAS ou 250 est donné
    forces_VKCAS = conditions_forces(VKCAS, "Vc", Choix_regime, Choix_regime_montee,
                      Choice_no_Engines, W1, isa, Hp, valeurs_atm[0], valeurs_atm[2], 
                      valeurs_atm[3], valeurs_atm[4], valeurs_atm[5], valeurs_atm[6], 0, "UP", 
                      1, "Nz", 0)
    
    forces_250 = conditions_forces(250, "Vc", Choix_regime, Choix_regime_montee,
                      Choice_no_Engines, W1, isa, Hp, valeurs_atm[0], valeurs_atm[2], 
                      valeurs_atm[3], valeurs_atm[4], valeurs_atm[5], valeurs_atm[6], 0, "UP", 
                      1, "Nz", 0)
    
    delta_V_kts = abs(forces_250[-1] - forces_VKCAS[-1])
    delta_V_fts = trans_kts_to_fts(delta_V_kts)
        
    # Accélération
    acc_Vc_moyen = ((forces_Vc_moyen[14] - forces_Vc_moyen[15])/W1)*cst.g # ft/s^2
    
    # Vraie vitesse
    V_Vc_moyen_kts = forces_Vc_moyen[-1] # kts
    V_Vc_moyen_fts = trans_kts_to_fts(V_Vc_moyen_kts) # fts
    
    # Delta t
    delta_t = delta_V_fts/acc_Vc_moyen # s
    
    # Delta d
    delta_d = (V_Vc_moyen_kts + Vwind_kts)*delta_t/3600 # NM
    
    # SFC
    SFC = 0.58 + (0.035 * Hp / 10000)
    
    # Calcul de Wf_moyen
    Wf_moyen = SFC * forces_Vc_moyen[14]
    
    # delta_fuel
    delta_fuel = Wf_moyen * (delta_t/3600)
    
    # W2
    W2 = W1 - delta_fuel
    
    return acc_Vc_moyen, V_Vc_moyen_fts, delta_t, delta_d, delta_fuel, W2
