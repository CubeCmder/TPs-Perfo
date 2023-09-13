import numpy as np

P_0 = 2116.22  # Standard SL Pressure [psf]
T_0 = 288.15  # Standard SL Temperature  [Kelvin]
RHO_0 = 0.002377  # Standard SL Density [Slugs/pi^3]
dT_dh = 0.0019812  # Tropospheric Temperature Gradient - Absolute Value [K/pi]

R = 96.0  # Gas constant [ft/K]
g = 32.174  # Gravitational acceleration [ft/s^2]

h_tr = 36089  # Altitude of the tropopause [ft]

delta_tr = 0.22336  # Pressure ratio at the tropopause

gamma = 1.4  # Atmospheric heat capacity ratio

# delta -> pressure ratio
# sigma -> density ratio
# theta -> temperature ratio

def temp_from_alt(h, ratio=False):
    """
    Get the atmospheric temperature at a given altitude.

    :param h: Altitude in [ft]
    :param ratio: if True returns temperature expressed as a ratio to SL temperature.
    :return: Temperature in [K]
    """
    if h < h_tr:
        theta = (T_0 - dT_dh * h) / T_0
    elif h_tr <= h <= 65617:
        theta = 216.65 / T_0
    else:
        raise Exception('Invalid altitude.')

    if not ratio:
        return theta * T_0
    else:
        return theta


def pressure_from_alt(h, ratio=False):
    """
    Get the atmospheric pressure at a given altitude.

    :param h: Altitude in [ft]
    :param ratio: if True returns pressure expressed as a ratio to SL pressure.
    :return: Pressure in [psf]
    """

    if -2001 <= h < h_tr:
        delta = (1 - dT_dh / T_0 * h) ** (1 / dT_dh / R)

    elif h_tr <= h <= 65617:
        delta = 0.22336 * np.exp(-(h - h_tr) / (R * temp_from_alt(h_tr)))

    else:
        raise Exception('Invalid altitude.')

    if not ratio:
        return delta * P_0
    else:
        return delta


def density_from_alt(h, ratio=False):
    """
    Get the atmospheric density at a given altitude.

    :param h: Altitude in [ft]
    :param ratio: if True returns density expressed as a ratio to SL density.
    :return: Density in [slugs/ft^3]
    """

    if h < h_tr:
        sigma = temp_from_alt(h, ratio=True) ** 4.2559
    elif h_tr <= h <= 65617:
        sigma = 0.29707 * np.exp(-(h - h_tr) / (R * temp_from_alt(h_tr)))
    else:
        raise Exception('Invalid altitude.')

    if not ratio:
        return sigma * RHO_0
    else:
        return sigma


def get_pressure_altitude(P):
    '''
    Get altitude equivalent to given ambient pressure.

    :param P: True ambient pressure [psf]
    :param h: True geometric altitude [ft]

    :return hp: Equivalent pressure altitude [ft]
    '''

    delta = P / P_0

    if delta > delta_tr:
        hp = (1 - delta ** (1 / 5.2559)) / (6.87535 * 10 ** -6)

    elif delta_tr >= delta >= 0.0540041:
        hp = h_tr - 20806 * np.log(4.477 * delta)

    else:
        raise Exception('Invalid altitude.')

    return hp


def get_delta_ISA(hp, T):
    '''
    Return temperature deviation from standard temperature at given altitude.

    :param hp: Pressure altitude [ft]
    :param T: Measured Temperature in [K] or in [°C]
    :return :
    '''

    T_ISA = temp_from_alt(hp)

    return T - T_ISA

def get_atmos_from_dISA(hp, dISA, ratio=True):
    """
    Get atmospheric properties given pressure altitude and
    temperature deviation from standard (dISA)

    :param hp: Pressure Altitude [ft]
    :param dISA: Temperature deviation in [K] or in [°C]
    :param ratio: If True returns properties expressed as ratios to SL values.

    :return: Atmospheric properties
    """

    T_std = temp_from_alt(hp)
    T_dISA = T_std + dISA

    delta = pressure_from_alt(hp, True)
    theta = T_dISA / T_0
    sigma = delta / theta

    if ratio:
        return delta, sigma, theta
    else:
        return delta*P_0, sigma*RHO_0, theta*T_0


