import numpy as np
try:
    from Constantes import *
    from Transformations import *
except:
    pass

def conditions_vitesse(Vitesse, Choice_vitesse, W, p, T, delta, theta, sigma, rho, LG, flap, S = 520, MAC = 8.286, K = 1):

    if Choice_vitesse == "Vc":
        
        # Vitesse du son
        a_0 = cst.a_0 # Noeuds
        a = a_0*(theta)**0.5
    
        # Pressions et nombre de Mach
        q_c = cst.p_0 * (((1+0.2*(Vitesse/a_0)**2)**(3.5)) - 1)
        mach = np.sqrt(5*((((q_c/p)+1)**(0.2857))-1))
        q = 1481.3 * delta * mach**2
        
        # Calcul de V
        V = mach * a_0 * np.sqrt(theta)
        V_fts = trans_kts_to_fts(V)
        V_ftmin = V_fts * 60
        
        # Temperature et pression totale
        T_t = T * (1 + 0.2*K*mach**2)
        p_t = q_c + p
        
        # Vitesse_equivalente
        V_e = V * np.sqrt(sigma)
        Ve_fts = trans_kts_to_fts(V_e)
        Ve_ftmin = Ve_fts * 60
        
        # Viscosite et nombre de Reynolds
        mu = (0.3125e-7*T**(1.5))/(T+120)
        rey = (rho*V_fts*MAC/mu)
        
    elif Choice_vitesse == "Ve":
        
        # Vitesse du son
        a_0 = cst.a_0 # Noeuds
        a = a_0*(theta)**0.5
        
        # # Pressions et nombre de Mach
        V_kts = Vitesse/sigma**0.5
        mach = V_kts/a
        q_c = p*((1+mach**2/5)**(1/0.2857)-1)
        
        # Calcul de V
        V = a*((5*((((q_c/p) + 1)**(0.2857))-1))**(1/2))
        V_fts = trans_kts_to_fts(V)
        V_ftmin = V_fts * 60
        
        # Temperature et pression totale
        T_t = T * (1 + 0.2*K*mach**2)
        p_t = q_c + p
        
        # Vitesses
        V_c = a_0*((5*((((q_c/cst.p_0) + 1)**(0.2857))-1))**(1/2))
        Vc_fts = trans_kts_to_fts(V_c)
        Vc_ftmin = Vc_fts * 60
        
        # Viscosite et nombre de Reynolds
        mu = (0.3125e-7*T**(1.5))/(T+120)
        rey = (rho*V_fts*MAC/mu)
        
    elif Choice_vitesse == "V":
        
        # Vitesse du son
        a_0 = cst.a_0 # Noeuds
        a = a_0*(theta)**0.5
        
        # Pressions et nombre de Mach
        q_c = p * (((1+0.2*(Vitesse/a)**2)**(3.5))- 1)
        mach = np.sqrt(5*((((q_c/p)+1)**(0.2857))-1))
        q = 1481.3 * delta * mach**2
        
        # Calcul de V
        V = Vitesse
        V_fts = trans_kts_to_fts(V)
        V_ftmin = V_fts * 60
        
        # Temperature et pression totale
        T_t = T * (1 + 0.2*K*mach**2)
        p_t = q_c + p
        
        # Vitesses
        V_c = a_0*((5*((((q_c/cst.p_0) + 1)**(0.2857))-1))**(1/2))
        Vc_fts = trans_kts_to_fts(V_c)
        Vc_ftmin = Vc_fts * 60
        
        V_e = Vitesse * np.sqrt(sigma)
        Ve_fts = trans_kts_to_fts(V_e)
        Ve_ftmin = Ve_fts * 60

        # Viscosite et nombre de Reynolds
        mu = (0.3125e-7*T**(1.5))/(T+120)
        rey = (rho*Vitesse*MAC/mu)
        
        
    elif Choice_vitesse == "Mach":
        
        # Vitesse du son
        a_0 = cst.a_0 # Noeuds
        a = a_0*(theta)**0.5
        
        # Pressions 
        q_c = p*((((Vitesse**2/5)+1)**(3.5))-1)
        q = 1481.3 * delta * (Vitesse)**2
        
        # Nombre de mach
        mach = Vitesse
        
        # Calcul de V
        V = mach * a_0 * np.sqrt(theta)
        V_fts = trans_kts_to_fts(V)
        V_ftmin = V_fts * 60
        
        # Calcul du oefficient de portance
        if Choice_Nzf == "f":
            Nz = 1/np.cos(trans_deg_to_rad(f))
            CL = (Nz*W)/(q*S)
            
        elif Choice_Nzf == "Nz":
            CL = (f*W)/(q*S)
            
        # Calcul de la portance
        L_lbs = q*S*CL
        
        # Temperature et pression totale
        T_t = T * (1 + 0.2*K*Vitesse**2)
        p_t = q_c + p
        
        # Vitesses
        V_c = a_0*((5*((((q_c/cst.p_0) + 1)**(0.2857))-1))**(1/2))
        Vc_fts = trans_kts_to_fts(V_c)
        Vc_ftmin = Vc_fts * 60
        
        V_e = V * np.sqrt(sigma)
        Ve_fts = trans_kts_to_fts(V_e)
        Ve_ftmin = Ve_fts * 60

        # Viscosite et nombre de Reynolds
        mu = (0.3125e-7*T**(1.5))/(T+120)
        rey = (rho*V_fts*MAC/mu)
        
    elif Choice_vitesse == "Vsr":
        
        # Vitesse du son
        a_0 = cst.a_0 # Noeuds
        a = a_0*(theta)**0.5
        
        # Calcul du CL
        if LG == "UP" and flap == 0:
            CLmax = 1.65
            
        elif LG == "DOWN" and flap == 0:
            CLmax = 1.60
            
        elif LG == "UP" and flap == 20:
            CLmax = 1.85
            
        elif LG == "DOWN" and flap == 20:
            CLmax = 1.80
            
        elif LG == "DOWN" and flap == 45:
            CLmax = 2.10
    
        CL = CLmax/(Vitesse)**2
        
        # # Calcul du nombre de mach (Module 4 - page 36)
        # if Choice_Nzf == "f":
        #     Nz = 1/np.cos(trans_deg_to_rad(f))
            
        # elif Choice_Nzf == "Nz":
        #     Nz = f
            
        mach = np.sqrt((W*1)/(delta*1481.3*CL*S))
        
        # Calcul du q
        q = 1481.3 * delta * mach**2
        
        # Calcul de la portance
        L_lbs = q*S*CL
        
        # Calcul de V
        V = mach * a_0 * np.sqrt(theta)
        V_fts = trans_kts_to_fts(V)
        V_ftmin = V_fts * 60
        
    return V
        