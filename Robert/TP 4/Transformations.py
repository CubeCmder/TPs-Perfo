import numpy as np

# Takes the temperature value in C and transforms into K
def trans_C_to_K(T):
    T = T + 273.15
    return T

# Takes the temperature value in K and transforms into C
def trans_K_to_C(T):
    T = T - 273.15
    return T

# m/s to kts
def trans_ms_to_kts(V):
    V = V*1.9438452
    return V

# kts to ft/s
def trans_kts_to_fts(V):
    V = V*1.68781
    return V

# ft/s to kts
def trans_fts_to_kts(V):
    V = V/1.68781
    return V

# degres to radians
def trans_deg_to_rad(V):
    V = V*np.pi/180
    return V