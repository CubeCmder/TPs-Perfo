import numpy as np
from atmos import *
from units import *

# Speed of sound at sea level
a0 = 661.48  # knots
P_0 = 2116.22  # Standard SL Pressure [psf]
T_0 = 288.15  # Standard SL Temperature  [Kelvin]
RHO_0 = 0.002377  # Standard SL Density [Slugs/pi^3]
Gamma = 0.14


def get_impact_pressure(p, mach):
    """
    Get compressible impact pressure (difference between total pressure and static pressure) at the given conditions.

    :param p: Local static pressure in [psf] - required
    :param mach: Aircraft mach number - required

    :return: compressible impact pressure (q_c) in [psf]
    """
    #qc = p * (((1 + mach**2/5)**(1/0.2857))-1)
    return p * ((1 + ((gamma - 1) / 2) * mach ** 2) ** (gamma / (gamma - 1)) - 1)


def get_total_pressure(p, mach):
    """
    Get total pressure at given conditions.

    :param p: Local static pressure in [psf] - required
    :param mach: Aircraft mach number - required

    :return: Total pressure in [psf]
    """

    qc = get_impact_pressure(p, mach)
    return qc + p


def get_dynamic_pressure(p, T=None, v=None, mach=None):
    """
    Get total pressure at given conditions.

    :param p: Local static pressure in [psf] - required
    :param T: Temperature [K]
    :param v: Aircraft true airspeed (TAS) - required if 'mach' missing
    :param mach: Aircraft mach number - required if 'v' missing

    :return: Dynamic pressure in [psf]
    """

    if v is None and mach is not None:
        return 1481.3 * p / P_0 * mach ** 2

    elif v is not None and mach is None and T is None:
        hp = get_pressure_altitude(p)
        rho = density_from_alt(hp)
        return 1 / 2 * rho * v ** 2
    elif v is not None and T is not None and mach is None:
        hp = get_pressure_altitude(p)
        dISA = get_delta_ISA(hp,T)
        p, rho, t = get_atmos_from_dISA(hp,dISA, False)
        return 1 / 2 * rho * v ** 2
    else:
        raise Exception('Ambiguous or erroneous function arguments.')


def get_total_temperature(T, mach, K=1):
    """

    :param T: Local atmospheric static temperature in [K] or [degC] - required
    :param mach: Aircraft mach number - required
    :param K: Recuperation factor of temperature sensor - default value of 1.0

    :return: Free-stream total temperature in [K] or [degC]
    """
    return T*(1+0.2*K*mach**2)


def get_SOS(temp=None, theta=None, knots=True):
    """
    Get the speed of sound at the given temperature or temperature ratio. Only one parameter required.

    :param temp: Absolute Temperature in [K]
    :param theta: Temperature Ratio
    :return: SOS in [knots]
    """
    if theta is None and temp is not None:
        if knots:
            return a0 * (temp / T_0) ** 0.5
        else:
            return knots2fps(a0 * (temp / T_0) ** 0.5)
    elif theta is not None and temp is None:
        if knots:
            return a0 * theta ** 0.5
        else:
            return knots2fps(a0 * theta ** 0.5)
    else:
        raise Exception('Ambiguous or erroneous function arguments.')


def get_mach(v=None, a=None, temp=None, qc=None, p=None):
    """
    Get aircraft mach number. Based on aircraft true speed and local speed of sound (sos).
    Local temperature can be given instead of local sos.

    :param v: Aircraft True Airspeed in [knots] - required if qc missing
    :param a: Local SOS in [knots] - required if qc and temp are missing
    :param temp: Local temperature in [K] - required if qc and a are missing
    :param qc: Compressible impact pressure in [psf] - required otherwise
    :param p: Static pressure in [psf] - required if qc is used
    :return: Aircraft mach number
    """
    if (v is not None and temp is None and a is not None) and qc is None:
        return v / a
    elif (v is not None and temp is not None and a is None) and qc is None:
        return v / get_SOS(temp=temp)
    elif qc is not None and p is not None:
        return (5 * ((qc / p + 1) ** 0.2857 - 1)) ** 0.5
    else:
        raise Exception('Ambiguous or erroneous function arguments.')


def get_calibrated_airspeed(p, mach, knots=True):
    """
    Get the calibrated airspeed (CAS or V_c) of an aircraft at the given conditions.

    :param p: Local static pressure in [psf] - required
    :param mach: Aircraft mach number - required

    :return: Calibrated airspeed in [knots]
    """

    qc = get_impact_pressure(p, mach)

    if knots:
        return a0 * (5 * ((qc / P_0 + 1) ** 0.2857 - 1)) ** 0.5
    else:
        return knots2fps(a0 * (5 * ((qc / P_0 + 1) ** 0.2857 - 1)) ** 0.5)

def get_mach_from_calibrated_airspeed(p, Vc):
    """
    Get Mach (M) from calibrated airspeed (V_c).

    :param p: Local static pressure in [psf] - required
    :param Vc: Equivalent airspeed in [knots] - required

    :return: Mach
    """
    qc = P_0 * ((1 + 0.2 * (Vc / a0) ** 2) ** 3.5 - 1)
    mach = ( (2 / (gamma - 1)) * ( (1 + qc/p) ** ( (gamma - 1) / gamma) -1))**0.5

    return mach

def get_pressure_from_mach_and_CAS(mach,Vc):
     p  = (P_0 * ((1 + 0.2 * (Vc / a0) ** 2) ** 3.5 - 1))/(((mach**2)/(2 / (gamma - 1))+1)**(1/((gamma - 1) / gamma))-1)
     #p = (1 + (P_0 * (1 + 0.2 * (Vc / a0) ** 2) ** 3.5 - 1)) / (((mach**2)/(2 / (gamma - 1))+1)**(1/((gamma - 1) / gamma)))

     return p

def get_equivalent_airspeed(p, mach, knots=True):
    """
    Get equivalent airspeed (V_e) of an aircraft at the given conditions.

    :param p: Local static pressure in [psf] - required
    :param mach: Aircraft mach number - required

    :return: Equivalent airspeed in [knots]
    """

    qc = get_impact_pressure(p, mach)

    if knots:
        return knots2fps((7 * p / RHO_0 * ((qc / p + 1) ** 0.2857 - 1)) ** 0.5, True)
    else:
        return (7 * p / RHO_0 * ((qc / p + 1) ** 0.2857 - 1)) ** 0.5

def get_mach_from_equivalent_airspeed(p, Ve, rho, a):
    """
    Get Mach (M) from equivalent airspeed (Ve).

    :param p: Local static pressure in [psf] - required
    :param Ve: Equivalent airspeed in [knots] - required

    :return: Mach
    """
    sigma = rho/RHO_0
    TAS = Ve/(sigma**0.5)
    mach = TAS/a

    return mach

def get_true_airspeed(p, mach, a=None, temp=None, knots=True):
    """
    Get the true airspeed (TAS or V) of an aircraft at the given conditions.

    :param p: Local static pressure in [psf] - required
    :param mach: Aircraft mach number - required
    :param a: Local SOS in [knots] or [fps] - required if 'temp' is missing
    :param temp: Local atmospheric temperature in [K] - required if 'a' is missing

    :return: True airspeed in [knots] if 'temp' is used, otherwise same unit as 'a'
    """

    qc = get_impact_pressure(p, mach)

    if a is None and temp is not None:
        a = get_SOS(temp)
    elif a is not None and temp is None:
        pass
    else:
        raise Exception('Ambiguous or erroneous function arguments.')

    if knots:
        return a * (5 * ((qc / p + 1) ** 0.2857 - 1)) ** 0.5
    else:
        return knots2fps(a * (5 * ((qc / p + 1) ** 0.2857 - 1)) ** 0.5)


def get_viscosity(p,T):
    """
    Get the dynamic viscosity (mu) of air at the given conditions.

    :param p: Local static pressure in [psf] - required

    :return: dynamic viscosity (lb*sec/pi2)
    """

    return 0.3125E-7*T**1.5/(T+120)

def get_reynolds(p,V,T,L=8.286):
    """
    Get the number of Reynolds (RN) at the given conditions.

    :param p: Local static pressure in [psf] - required
    :param V: True air speed
    :param L: MAC (mean aerodynamic chord)

    :return: Reynold Number
    """
    mu = get_viscosity(p,T)
    V = knots2fps(V)
    h = get_pressure_altitude(p)
    dISA = get_delta_ISA(h, T)
    p, rho, T = get_atmos_from_dISA(h, dISA, False)

    return rho*V*L/mu

def get_lift_coefficient(p,V,W,T,S=520,N_z = 1):
    """
    Get the number of Reynolds (RN) at the given conditions.

    :param p: Local static pressure in [psf] - required
    :param V: True airspeed [knots]
    :param S: Wing surface
    :param W: Aircraft weight
    :param N_z: Coefficient

    :return : Lift coefficient
    """
    L = N_z*W
    V = knots2fps(V)
    q = get_dynamic_pressure(p,T,V)

    return L/(q*S)


