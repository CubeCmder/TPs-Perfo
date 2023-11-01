import numpy as np
import math as mt
from scipy.optimize import fsolve

try:
    from Atmosphere import *
    from Forces import *
    from Vitesses import *
    
    from Acceleration import *
    from Interval import *
except:
    pass

def SAR(Hp, isa, W_lbs, Vitesse):
    
    valeurs_atm = conditions_atm(Hp, isa, "ISA")
    
    forces = conditions_forces(Vitesse, "V", "MCR", "Mach_constant",
                      "AEO", W_lbs, isa, Hp, valeurs_atm[0], valeurs_atm[2], 
                      valeurs_atm[3], valeurs_atm[4], valeurs_atm[5], valeurs_atm[6], 0, "UP", 
                      1, "Nz", 0)
    
    T_req = W_lbs/(forces[14]/forces[15])
    SFC = 0.58 + (0.035 * Hp/10000)
    
    return Vitesse/(T_req*SFC)