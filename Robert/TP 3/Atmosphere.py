import numpy as np
try:
    from Constantes import *
    from Transformations import *
except:
    pass

# Value defines whether ISA or T is chosen depending on
# value of Choice ("ISA", "T")
def conditions_atm(Hp, Value, Choice):
    if Choice == "T":
        T_kelvin = trans_C_to_K(Value)
        if Hp <= 36089:
            theta = trans_C_to_K(Value)/cst.T_0 
            ISA = Value - (trans_K_to_C(cst.T_0) - cst.lam*Hp)
            delta = (1-(Hp*cst.const_random))**5.2559
            sigma = delta/theta
            P = delta*cst.p_0
            rho = sigma * cst.rho_0
            
            return T_kelvin, ISA, theta, delta, sigma, P, rho

        else:
            theta = trans_C_to_K(Value)/cst.T_0 
            ISA = Value - (-56.50) # Above 36089 pi, T_ISA = -56.50C
            delta = (1/4.477)*np.exp((36089-Hp)/20806)
            sigma = delta/theta
            P = delta*cst.p_0
            rho = sigma * cst.rho_0
            
            return T_kelvin, ISA, theta, delta, sigma, P, rho
        
    elif Choice == "ISA":
        if Hp <= 36089:
            Value = Value + (trans_K_to_C(cst.T_0) - cst.lam*Hp)
            T_kelvin = trans_C_to_K(Value)
            theta = trans_C_to_K(Value)/cst.T_0 
            delta = (1-(Hp*cst.const_random))**5.2559
            sigma = delta/theta
            P = delta*cst.p_0
            rho = sigma * cst.rho_0
            
            return T_kelvin, Value, theta, delta, sigma, P, rho
        
        else:
            Value = Value + (-56.50) # Above 36089 pi, T_ISA = -56.50C
            T_kelvin = trans_C_to_K(Value)
            theta = trans_C_to_K(Value)/cst.T_0 
            delta = (1/4.477)*np.exp((36089-Hp)/20806)
            sigma = delta/theta
            P = delta*cst.p_0
            rho = sigma * cst.rho_0
            
            return T_kelvin, Value, theta, delta, sigma, P, rho