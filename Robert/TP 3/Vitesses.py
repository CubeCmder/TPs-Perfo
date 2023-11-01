import numpy as np
try:
    from Constantes import *
    from Transformations import *
except:
    pass

def conditions_vitesse(Vitesse, Choice_vitesse, W, p, T, delta, theta, sigma, rho, S = 520, MAC = 8.286, K = 1):
    if Choice_vitesse == "Vc":
        
        # Vitesse du son
        a_0 = cst.a_0 # Noeuds
        a = a_0*(theta)**0.5

        # Pressions et nombre de Mach
        q_c = cst.p_0 * (((1+0.2*(Vitesse/a_0)**2)**(3.5)) - 1)
        mach = np.sqrt(5*((((q_c/p)+1)**(0.2857))-1))
        q = 1481.3 * delta * mach**2
        
        # Coefficient de portance
        C_L = W/(q*S)
        
        # Temperature et pression totale
        T_t = T * (1 + 0.2*K*mach**2)
        p_t = q_c + p
        
        # Vitesses
        V = mach * a_0 * np.sqrt(theta)
        V_fts = trans_kts_to_fts(V)
        V_e = V * np.sqrt(sigma)
        
        # Viscosite et nombre de Reynolds
        mu = (0.3125e-7*T**(1.5))/(T+120)
        rey = (rho*V_fts*MAC/mu)
        
        return a, mach, V, V_e, p_t, q, q_c, T_t, mu, rey, C_L
    
    elif Choice_vitesse == "Ve":
        
        # Vitesse du son
        a_0 = cst.a_0 # Noeuds
        a = a_0*(theta)**0.5
        
        # Pressions et nombre de Mach
        q_c = p * ((((V_fts**2*cst.rho_0)/(7*p) + 1)**(3.5)) - 1)
        mach = np.sqrt(5*((((q_c/p)+1)**(0.2857))-1))
        q = 1481.3 * delta * mach**2
        
        # Coefficient de portance
        C_L = W/(q*S)
        
        # Temperature et pression totale
        T_t = T * (1 + 0.2*K*mach**2)
        p_t = q_c + p
        
        # Vitesses
        V = a*((5*((((q_c/p) + 1)**(0.2857))-1))**(1/2))
        V_c = a_0*((5*((((q_c/cst.p_0) + 1)**(0.2857))-1))**(1/2))
        V_fts = trans_kts_to_fts(V)
        
        # Viscosite et nombre de Reynolds
        mu = (0.3125e-7*T**(1.5))/(T+120)
        rey = (rho*V_fts*MAC/mu)
        
        return a, mach, V, V_c, p_t, q, q_c, T_t, mu, rey, C_L
    
    elif Choice_vitesse == "V":
        
        # Vitesse du son
        a_0 = cst.a_0 # Noeuds
        a = a_0*(theta)**0.5
        
        # Pressions et nombre de Mach
        q_c = p * (((1+0.2*(Vitesse/a)**2)**(3.5))- 1)
        mach = np.sqrt(5*((((q_c/p)+1)**(0.2857))-1))
        q = 1481.3 * delta * mach**2
        
        # Coefficient de portance
        C_L = W/(q*S)
        
        # Temperature et pression totale
        T_t = T * (1 + 0.2*K*mach**2)
        p_t = q_c + p
        
        # Vitesses
        V_c = a_0*((5*((((q_c/cst.p_0) + 1)**(0.2857))-1))**(1/2))
        V_e = Vitesse * np.sqrt(sigma)
        Vitesse = trans_kts_to_fts(Vitesse)

        # Viscosite et nombre de Reynolds
        mu = (0.3125e-7*T**(1.5))/(T+120)
        rey = (rho*Vitesse*MAC/mu)
        
        return a, mach, V_c, V_e, p_t, q, q_c, T_t, mu, rey, C_L
    
    elif Choice_vitesse == "Mach":
        
        # Vitesse du son
        a_0 = cst.a_0 # Noeuds
        a = a_0*(theta)**0.5
        
        #Pressions 
        q_c = p*((((Vitesse**2/5)+1)**(3.5))-1)
        q = 1481.3 * delta * Vitesse**2
           
        # Coefficient de portance
        C_L = W/(q*S)
        
        # Temperature et pression totale
        T_t = T * (1 + 0.2*K*Vitesse**2)
        p_t = q_c + p
        
        # Vitesses
        V = Vitesse * a_0 * np.sqrt(theta)
        V_c = a_0*((5*((((q_c/cst.p_0) + 1)**(0.2857))-1))**(1/2))
        V_e = V * np.sqrt(sigma)
        V_fts = trans_kts_to_fts(V)

        # Viscosite et nombre de Reynolds
        mu = (0.3125e-7*T**(1.5))/(T+120)
        rey = (rho*V_fts*MAC/mu)
        
        return a, V, V_e, V_c, p_t, q, q_c, T_t, mu, rey, C_L
        