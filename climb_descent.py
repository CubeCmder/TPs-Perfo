from typing import Any
import numpy as np
import warnings
from velocities import *
from atmos import *
from units import knots2fps
def climb_descent(aircraft, R, alt, T, W, RM, n_engines, aoa, flap, landing_gear, cg, CAS=None, Mach=None):
    """
    Function that returns all information about climb or descent phase of flight

    :param R: Régime de vol ["CAS","Mach"]
    :param alt: altitude pression [ft]
    :param T: temperature [k]
    :param W: Weight [lb]
    :param RM: Param moteur [lb]
    :param n_engines: Number of engines operative
    :param aoa: Angle of attack
    :param flap: Flap angle
    :param landing_gear: landing gear position
    :param CAS: Calibrated air speed (constant if given with regime CAS)
    :param Mach: Mach (constant if given with regime Mach)
    :return: gradient, ROC, ROCp, AF, a_xfp
    """

    p = pressure_from_alt(alt)
    deltaISA = get_delta_ISA(alt, T)
    T_std = T - deltaISA

    if CAS is None and Mach is not None:
        pass
    elif CAS is not None and Mach is None:
        Mach = get_mach_from_calibrated_airspeed(p, CAS)
    else:
        raise Exception('Cas or Mach must be given')

    if abs(n_engines - aircraft.n_engines) < 1:
        OEI = 0
    else:
        OEI = 1

    V = get_true_airspeed(p, Mach, temp=T)
    V = knots2fps(V)

    q = get_dynamic_pressure(p, T, mach=Mach)

    Thrust = aircraft.get_thrust(RM, alt, Mach, T, n_engines=n_engines)
    CL = aircraft.get_lift_coefficient(Nz=1, aoa=aoa, flap_angle=flap, gear_up=landing_gear, cg=cg)
    CD = aircraft.get_drag_coefficient(CL, flap, Mach, LDG=not landing_gear, NZ=1, OEI=OEI, q=q, thrust=Thrust)['CDtot']
    D = q * CD * aircraft.S

    AF = get_AF(R, alt, Mach, T, T_std)
    ROC = get_ROC(V, Thrust, D, W, AF)
    ROCp = get_ROCp(ROC, T, T_std)
    gradient = get_gradient(AF,Thrust,W,CD,CL)

    a_xfp = get_accel_xfp(Thrust, W, CD, CL, gradient)

    return gradient, ROC, ROCp, AF, a_xfp

def get_gradient(AF,T,W,CD,CL):
    gradient = (T/W-CD/CL)/(1+AF)
    return gradient

def get_AF(R, alt, Mach, t , t_std):
    phi = (1/0.7/Mach**2)*((1+0.2*Mach**2)**3.5-1)/((1+0.2*Mach**2)**2.5)

    if alt<36089:
        if R == "CAS" :
            AF = 0.7 * Mach**2 * (phi-0.190263*(t_std/t))
        elif R == "Mach":
            AF = -0.133184*Mach**2 * (t_std/t)
    else:
        if R == "CAS":
            AF = 0.7*Mach**2 * phi
        elif R == "Mach":
            AF = 0

    return AF

def get_ROC(V, T, D, W, AF):
    ROC = ((V*(T-D))/W)/(1+AF)
    return ROC*60

def get_ROCp(ROC,t ,t_std):
    ROCp = ROC*(t_std/t)
    return ROCp

def get_accel_xfp(T, W, CD, CL, gamma):
    """
    The acceleration of the aircraft in the direction of the flight trajectory.

    :param T: Total Thrust
    :param W: Aircraft Weight
    :param CD: Drag coefficient
    :param CL: Lift Coefficient
    :param gamma: Climb gradient
    :return: accel
    """
    accel = (T/W - CD/CL) - gamma
    return accel
