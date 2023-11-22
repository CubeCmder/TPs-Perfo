import numpy as np
import math as mt
from scipy.interpolate import interp1d

try:
    from Constantes import *
    from Transformations import *
    from Atmosphere import *
except:
    pass

def conditions_forces(Vitesse, Choice_vitesse, Choice_regime, Choice_regime_climb, Choice_no_Engines, W, isa, hp, T, \
                       theta, delta, sigma, p, rho, flap, LG, f, Choice_Nzf, CG, K = 1, S = 520, MAC = 8.286, LT = 40.56):
    
    # Partie A
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
        
        # Calcul du oefficient de portance
        if Choice_Nzf == "f":
            Nz = 1/np.cos(trans_deg_to_rad(f))
            CL = (Nz*W)/(q*S)
            
        elif Choice_Nzf == "Nz":
            CL = (f*W)/(q*S)
            
        # Calcul de la portance
        L_lbs = q*S*CL
        
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
        
        # Pressions et nombre de Mach
        q_c = p * ((((Vitesse**2*cst.rho_0)/(7*p) + 1)**(3.5)) - 1)
        mach = np.sqrt(5*((((q_c/p)+1)**(0.2857))-1))
        q = 1481.3 * delta * mach**2
        
        # Calcul de V
        V = a*((5*((((q_c/p) + 1)**(0.2857))-1))**(1/2))
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
        
        # Calcul du oefficient de portance
        if Choice_Nzf == "f":
            Nz = 1/np.cos(trans_deg_to_rad(f))
            CL = (Nz*W)/(q*S)
            
        elif Choice_Nzf == "Nz":
            CL = (f*W)/(q*S)
            
        # Calcul de la portance
        L_lbs = q*S*CL
        
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
        
        # Calcul du nombre de mach (Module 4 - page 36)
        if Choice_Nzf == "f":
            Nz = 1/np.cos(trans_deg_to_rad(f))
            
        elif Choice_Nzf == "Nz":
            Nz = f
            
        mach = np.sqrt((W*Nz)/(delta*1481.3*CL*S))
        
        # Calcul du q
        q = 1481.3 * delta * mach**2
        
        # Calcul de la portance
        L_lbs = q*S*CL
        
        # Calcul de V
        V = mach * a_0 * np.sqrt(theta)
        V_fts = trans_kts_to_fts(V)
        V_ftmin = V_fts * 60
    
    # Calcul du Thrust
    if Choice_regime == "Poussee_totale":
        T_lbs = 4500
        
    if Choice_regime == "MTO" or Choice_regime == "GA" or Choice_regime == "MCT":
        T_lbs = 8775 - 0.1915*hp - (8505-0.195*hp)*mach
        
        if isa > 15:
            # Note: Above ISA+15, MTOFN reduces by 1 % per degree C
            Delta_ISA = isa - 15
            T_lbs = T_lbs*(100-Delta_ISA)/100
            
        elif Choice_regime == "MCT":
            T_lbs = T_lbs*0.9             
                
    elif Choice_regime == "MCL" or Choice_regime == "MCR":
        T_lbs = 5690 - 0.0968*hp - (1813-0.0333*hp)*mach
    
        if isa > 10:
            #Note: Above ISA+10, MCLFN reduces by 1 % per degree C
            Delta_ISA = isa - 10
            T_lbs = T_lbs*(100-Delta_ISA)/100
        
        if Choice_regime == "MCR":
            T_lbs = T_lbs*0.98
                
    elif Choice_regime == "Idle":
        T_lbs = 600-1000*mach
                
    # Determiner le trust selon le nombre de moteurs
    if Choice_no_Engines == "AEO" and Choice_regime != "Poussee_totale":
        T_lbs = 2*T_lbs
        
    elif Choice_no_Engines == "OEI" and Choice_regime != "Poussee_totale":
        T_lbs = T_lbs
    
    # Calcul de CDp et Cdi
    if flap == 0:
        CDp = 0.0206
        K = 0.0364
        CDi = K*CL**2
        
    elif flap == 20:
        CDp = 0.0465
        K = 0.0334
        CDi = K*CL**2
        
    elif flap == 45:
        CDp = 0.1386
        K = 0.0301
        CDi = K*CL**2
        
    # Calcul de CDcomp
    if mach >= 0 and mach <= 0.6:
        CDcomp = 0
        
    elif mach > 0.6 and mach <= 0.78:
        CDcomp = (0.0508 - 0.1748*mach + 0.1504*mach**2)*CL**2
        
    elif mach > 0.78 and mach <= 0.85:
        CDcomp = (-99.3434 + 380.888*mach - 486.8*mach**2 + 207.408*mach**3)*CL**2
        
    # Calcul de CDcntl
    if Choice_no_Engines == "AEO":
        CDcntl = 0
        CDwm = 0
        
    elif Choice_no_Engines == "OEI":   
    #DRAG INCREMENTS DURING OPERATION WITH ONE ENGINE OUT (in the air)
        CT = T_lbs/(q*S)
        CDcntl = 0.1*CT**2
        CDwm = 0.003
        
    # Calcul de CD et du D Total
    CD = CDp + CDi + CDcomp + CDcntl + CDwm
    
    if LG == "UP":
        CD = CD
        
    elif LG == "DOWN":
        #DRAG INCREMENT WITH LANDING GEAR DOWN
        CD = CD + 0.02
        
    # Calcul de la trainee
    D_lbs = q*S*CD
 
    # Calcul de l'angle d'attaque (formule fournie dans le module 4 - page 24)
    CLfwd = CL*(1 + (MAC/LT)*(0.09 - CG))
    if LG == "UP":
        if flap == 0:
            aoa = (CLfwd - 0.05)/0.1
        elif flap == 20:
            aoa = (CLfwd - 0.25)/0.1
        elif flap == 45:
            aoa = (CLfwd - 0.55)/0.1
            
    if LG == "DOWN":
        if flap == 0:
            aoa = (CLfwd - 0)/0.1
        elif flap == 20:
            aoa = (CLfwd - 0.20)/0.1
        elif flap == 45:
            aoa = (CLfwd - 0.50)/0.1

    # Calcul de CLsw (base sur AOAsw)
    if LG == "UP":
        if flap == 0:
            CLsw = 0.05 + 0.1*14.7
            CLsw = CLsw/(1 + (MAC/LT)*(0.09 - CG))
        elif flap == 20:
            CLsw = 0.25 + 0.1*14.6
            CLsw = CLsw/(1 + (MAC/LT)*(0.09 - CG))
        elif flap == 45:
            CLsw = 0.55 + 0.1*14.4
            CLsw = CLsw/(1 + (MAC/LT)*(0.09 - CG))
            
    if LG == "DOWN":
        if flap == 0:
            CLsw = 0.00 + 0.1*14.7
            CLsw = CLsw/(1 + (MAC/LT)*(0.09 - CG))
        elif flap == 20:
            CLsw = 0.20 + 0.1*14.6
            CLsw = CLsw/(1 + (MAC/LT)*(0.09 - CG))
        elif flap == 45:
            CLsw = 0.50 + 0.1*14.4
            CLsw = CLsw/(1 + (MAC/LT)*(0.09 - CG))

    # Calcul de Nz sw et de phi sw (base sur CLsw)
    Nzsw = CLsw*q*S/W
    if 1/Nzsw >= -1 and 1/Nzsw <= 1:
        Phisw = mt.acos(1/Nzsw) * (180/np.pi)
    else:
        Phisw = 0

    # Calcul de Nz buffet (base sur CLbuffet)
    machs = [0.2750, 0.3000, 0.3250, 0.3500, 0.3750, 0.4000, 0.4250, 0.4500,
              0.4750, 0.5000, 0.5250, 0.5500, 0.5750, 0.6000, 0.6250, 0.6500,
              0.6750, 0.7000, 0.7250, 0.7500, 0.7750, 0.8000, 0.8250, 0.8500,
              0.8750, 0.9000]
    CLbuffets = [1.342, 1.3199, 1.2974, 1.2667, 1.2310, 1.1930, 1.1551, 1.1191,  
                1.0863, 1.0577, 1.0337, 1.0142, 0.9989, 0.9868, 0.9764, 0.9659,
                0.9530, 0.9349, 0.9085, 0.8698, 0.8149, 0.7391, 0.6373, 0.5039,
                0.3330, 0.118]
    
    if mach >= machs[0] and mach <= machs[-1]:
        CLbuffet = interp1d(machs, CLbuffets)(mach)/(1 + (MAC/LT)*(0.09 - CG))
        Nzbuffet = CLbuffet*q*S/W
    else:
        CLbuffet = 0
        Nzbuffet = 0
    
    # Partie B
    # Calcul du facteur phi
    phi = (1/(0.7*mach**2))*((1+0.2*mach**2)**3.5-1)/((1+0.2*mach**2)**2.5)
    
    # Calcul de la temperature standard
    T_std = conditions_atm(hp, 0, "ISA")[0] # K
    
    # Calcul du facteur d'acceleration
    if Choice_regime_climb == "Mach_constant":
        if hp <= 36089:
            FA = -0.133184 * mach**2 * (T_std/T)
        elif hp > 36089:
            FA = 0
            
    if Choice_regime_climb == "CAS_constant":
        if hp <= 36089:
            FA = 0.7 * mach**2 * (phi - 0.190263*(T_std/T))
        elif hp > 36089:
            FA = 0.7 * mach**2 * phi
            
    # Calcul du gradient de montee (Module 5 - page 9)
    gamma_montee = (T_lbs/W - CD/CL)/(1 + FA)
    
    # Calcul du taux de montee - ROC (Module 5 - page 6)
    taux_montee = V_ftmin * gamma_montee
    
    # Calcul du taux de montee pression - ROCp (Module 5 - page 18)
    taux_montee_pression = taux_montee * (T_std/T)
    
    # Calcul de l'acceleration selon l'axe de la trajectoire de vol (Module 5 - page)
    axfp = (T_lbs/W - CD/CL) - gamma_montee
 
    return CD, CL, L_lbs, CL/CD, CDp, q*S*CDp, CDi, q*S*CDp, CDcomp, q*S*CDcomp, CDcntl, q*S*CDcntl, CDwm, q*S*CDwm, T_lbs, D_lbs, aoa, Nzsw, Phisw, Nzbuffet, gamma_montee*100, taux_montee, taux_montee_pression, FA, axfp, V, q