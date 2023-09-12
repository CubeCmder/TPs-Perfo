import numpy as np
from atmos import *

# Speed of sound at sea level
a0 = 661.48  # knots


def get_impact_pressure(p, mach):
    """
    Get compressible impact pressure (difference between total pressure and static pressure) at the given conditions.

    :param p: Local static pressure in [psf] - required
    :param mach: Aircraft mach number - required

    :return: compressible impact pressure (q_c) in [psf]
    """

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


def get_dynamic_pressure(p, v=None, mach=None):
    """
    Get total pressure at given conditions.

    :param p: Local static pressure in [psf] - required
    :param v: Aircraft true airspeed (TAS) - required if 'mach' missing
    :param mach: Aircraft mach number - required if 'v' missing

    :return: Dynamic pressure in [psf]
    """

    if v is None and mach is not None:
        return 1481.3 * p / P_0 * mach ** 2

    elif v is not None and mach is None:
        hp = get_pressure_altitude(p)
        rho = density_from_alt(hp)
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


def get_SOS(temp=None, theta=None):
    """
    Get the speed of sound at the given temperature or temperature ratio. Only one parameter required.

    :param temp: Absolute Temperature in [K]
    :param theta: Temperature Ratio
    :return: SOS in [knots]
    """
    if theta is None and temp is not None:
        return a0 * (temp / T_0) ** 0.5
    elif theta is not None and temp is None:
        return a0 * theta ** 0.5
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
        return (5 * ((qc / p + 1) * 0.2857 - 1)) ** 0.5
    else:
        raise Exception('Ambiguous or erroneous function arguments.')


def get_calibrated_airspeed(p, mach):
    """
    Get the calibrated airspeed (CAS or V_c) of an aircraft at the given conditions.

    :param p: Local static pressure in [psf] - required
    :param mach: Aircraft mach number - required

    :return: Calibrated airspeed in [knots]
    """

    qc = get_impact_pressure(p, mach)

    return a0 * (5 * ((qc / P_0 + 1) * 0.2857 - 1)) ** 0.5

def get_mach_from_calibrated_airspeed(p, V_e):
    h = get_pressure_altitude(p)
    sigma = density_from_alt(h, True)
    V = V_e/sigma**(0.5)
    temp = temp_from_alt(h)
    a = get_SOS(temp)
    qc = get_dynamic_pressure(p)
    Mach = get_mach(V, a, temp, qc, p)

    return Mach

def get_equivalent_airspeed(p, mach):
    """
    Get equivalent airspeed (V_e) of an aircraft at the given conditions.

    :param p: Local static pressure in [psf] - required
    :param mach: Aircraft mach number - required

    :return: Equivalent airspeed in [knots]
    """

    qc = get_impact_pressure(p, mach)

    return (7 * p / P_0 * ((qc / p + 1) ** 0.2857 - 1)) ** 0.5


def get_true_airspeed(p, mach, a=None, temp=None):
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

    return a * (5 * ((qc / P_0 + 1) * 0.2857 - 1)) ** 0.5
