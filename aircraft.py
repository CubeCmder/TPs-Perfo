from typing import Any
import numpy as np
import warnings

import atmos
from velocities import *
from atmos import *
from units import *
from climb_descent import *


class Aircraft():
    CLMAX_F45_GD: Any

    def __init__(self):
        # WEIGHTS - [lb]
        self.MWE = 28000.0
        self.OWE = 31500.0
        self.MZFW = 44000.0
        self.MLW = 47000.0
        self.MTOW = 53000.0
        self.MRW = 53250.0
        self.MXFUEL = 15000.0

        # CG
        self.FWDCG = 9.0 / 100

        # AIRCRAFT GEOMETRY
        self.S = 520.0  # Wing area (ft^2)
        self.span = 67.85  # Wing span (ft)
        self.mac = 8.286  # Mean aerodynamic chord (ft)
        self.zmac = 5.016  # Height of the MAC above ground (ft)
        self.lt = 40.56  # Tail arm (ft)
        self.yeng = 7.00  # Thrust line moment arm relative to C.G. (ft)
        self.zeng = 121.35  # Thrust line water line
        self.xnw = 228.40  # Fuselage station at nose landing gear
        self.xuc = 676.41  # Fuselage station at main landing gear
        self.xcg = 630.70  # Fuselage station at center of gravity at 9% MAC
        self.xcp = 646.65  # Fuselage station at aerodynamic center
        self.zcp = 52.25  # Water line station at aerodynamic center
        self.zcg = 94.50  # Center of gravity water line
        self.zgl = 9.55  # Static ground line water line

        # AIRCRAFT CONFIGURATION - notation => "key" = [flap angle, LG state (0=up, 1=deployed)]
        self.config = {"cruise": [0, 0],
                       "takeoff": [20, 0],
                       "approach": [20, 0],
                       "landing": [45, 1]}

        # AIRCRAFT PROPULSION
        self.K = 1  # Temperature probe recovery factor

        # AIRCRAFT ENGINE DATA
        self.n_engines = 2  # number of engines

        #                **********************
        #                     Speed limits
        #                **********************

        self.V_MO = 330  # Maximum CAS with no flaps deployed [kts]
        self.M_MO = 0.85  # Maximum Mach number with no flaps deployed
        self.V_FE_F20 = 210  # Maximum CAS with flaps deployed at 20° [kts]
        self.V_FE_F45 = 180  # Maximum CAS with flaps deployed at 20° [kts]
        self.V_LO = self.V_LE = 220  # Maximum CAS with LG deployed [kts]
        self.V_mca = 95  # minimum control speed in the air (take-off) (KCAS)
        self.V_mcl = 92  # minimum control speed in the air (landing) (KCAS)
        self.V_mcg = 90  # minimum control speed on the ground (take-off) (KCAS)
        self.V_1mcg = 95  # minimum v1 based on engine failure at vmcg (KCAS)

        #  *********************************************************
        #                 MAXIMUM LIFT COEFFICIENT
        #  *********************************************************

        # FXX : XX defines flap setting
        # GU / GD : landing gear up / landing gear down

        self.CLMAX_F00_GU = 1.65
        self.CLMAX_F00_GD = 1.60
        self.CLMAX_F20_GU = 1.85
        self.CLMAX_F20_GD = 1.80
        self.CLMAX_F45_GD = 2.10

        #  *********************************************************
        #          FUSELAGE ANGLE-OF-ATTACK AT STALL WARNING
        #  *********************************************************

        self.fuse_AOA_SW = {0: 14.7, 20: 14.6, 45: 14.4}  # Keys: Flap angle [°]; Values: AoA @ Stall warning [°]

        #  *********************************************************
        #                   DRAG COEFFICIENTS
        #  *********************************************************

        self.d_CD_LG = 0.02  # DRAG INCREMENT WITH LANDING GEAR DOWN
        self.d_CD_WM = 0.0030  # (windmilling drag coefficient) - OEI

        # *********************************************************
        #                   Speed and time increments
        #  *******************************************************
        # self.twdv = #(Climb gradient at 35 ft – θ)(based on V2, gear up and AF = 0)
        # self.dvrvl = #Speed spread from rotation to lift - off(normal rotation)(ft / s)
        # self.dvlo35 = #Speed spread from lift-off to 35 ft(ft / s)
        # self.dvlo15 = #Speed spread from lift-off to 15 ft(ft / s)

        # self.twdt = #(Climb gradient at 35 ft – θ)(based on V2, gear up and AF = 0)
        # self.dtvrvl = #Time between Vr and Vlo
        # self.dtvlo35 = #Time between Vlo and 35 ft
        # self.dtvlo15 = #Time between Vlo and 15 ft

        self.ndv = 8
        self.twdv = [0.000, 0.010, 0.020, 0.100, 0.120, 0.200, 0.400, 0.60]
        self.dvrvl = [1.80, 2.30, 2.75, 6.50, 7.35, 10.80, 19.80, 30.0]
        self.dvlo35 = [0.00, 0.00, 1.00, 9.00, 11.00, 16.70, 31.00, 45.0]
        self.dvlo15 = [0.00, 0.00, 1.00, 8.05, 9.80, 13.00, 21.00, 29.0]

        self.ndt = 11
        self.twdt = [0.0, 0.02, 0.03, 0.05, .065, .075, 0.10, 0.18, 0.2, 0.4, 0.6]
        self.dtvrvl = [2.35, 2.25, 2.19, 2.09, 2.01, 1.96, 1.83, 1.41, 1.30, 1.30, 1.30]
        self.dtvlo35 = [17.00, 10.00, 6.60, 5.50, 5.10, 4.90, 4.50, 3.80, 3.80, 3.80, 3.80]
        self.dtvlo15 = [12.00, 7.20, 4.70, 4.10, 3.92, 3.80, 3.50, 2.55, 2.55, 2.55, 2.55]

        #  *********************************************************
        #                       Brake Coefficient
        #  *********************************************************
        self.kemax = 16.7  # max demonstrated brake energy per brake (million ft-lb)
        self.fuse_plug_limit = 11.5  # fuse plug limit brake energy per brake (million ft-lb)
        self.mu_dry = 0.400  # Airplane braking coefficient on dry runway
        self.rollmu = 0.0200  # Rolling coefficient of friction
        self.mu_wet = 0.225  # Brake mu on wet runways
        self.mu_snow = 0.200  # Brake mu on compacted snow covered runway
        self.mu_ice = 0.050  # Brake mu on ice covered runway
        self.Pm = 165  # Main tire pressure (psi)
        self.tiremax = 210  # Maximum allowable tire speed (mph)
        self.nb_brake = 4  # Number of brakes

        #  *********************************************************
        #                       T.O / R.T.O. Time Delays
        #  *********************************************************
        self.dtdwm = 0  # Time to reach windmilling drag level following an engine cut (sec)
        self.dtapr = 0  # Time to reach apr thrust (sec)
        self.dtrec = 1  # Engine failure recognition time (sec)
        self.dtidle = 0  # Time to reach idle thrust after throttles chopped (sec)
        self.tbrake = 0  # Time from V1 to brake application (sec)
        self.tpower = 0  # Time from V1 to throttle chop to idle( sec)
        self.tdump = 0  # Time from V1 to GLD fully deployed (sec)

        #  *********************************************************
        #        Landing Delay Time/Decel and parametric param
        #  *********************************************************
        self.dtdelay = 1.2
        self.adelay = 4
        self.tair_min = 3.48
        self.at = 2.6337
        self.bt = 0.4443
        self.ct = 0.3584
        self.av = 1.1097
        self.bv = -0.0043
        self.cv = -0.0027



    #  ===============================================================================
    def _induced_drag_efficiency_factor(self, flap_angle):
        """
        Induced drag K-factor

        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :return:
        """

        if flap_angle == 0:
            K = 0.0364
        elif flap_angle == 20:
            K = 0.0334
        elif flap_angle == 45:
            K = 0.0301
        else:
            raise Exception('arg <flap_angle> can only be 0, 20 or 45 degrees.')

        return K

    def CL_at_buffet_vs_mach(self, mach):
        """
        LIFT COEFFICIENT AT BUFFET (CL_BUFFET) AS A FUNCTION ON MACH
        DATA BASED ON 9 % MAC AND FLAP 0

        :param mach: Mach number
        :return:
        """

        x_pts = [0.2750, 0.3000, 0.3250, 0.3500, 0.3750, 0.4000, 0.4250, 0.4500,
                 0.4750, 0.5000, 0.5250, 0.5500, 0.5750, 0.6000, 0.6250, 0.6500,
                 0.6750, 0.7000, 0.7250, 0.7500, 0.7750, 0.8000, 0.8250, 0.8500,
                 0.8750, 0.9000]

        y_pts = [1.3424, 1.3199, 1.2974, 1.2667, 1.2310, 1.1930, 1.1551, 1.1191, 1.0863,
                 1.0577, 1.0337, 1.0142, 0.9989, 0.9868, 0.9764, 0.9659, 0.9530, 0.9349,
                 0.9085, 0.8698, 0.8149, 0.7391, 0.6373, 0.5039, 0.3330, 0.118]

        return np.interp(mach, x_pts, y_pts)

    def speed_spread_vs_climb_gradient(self, climb_gradient_35, type ):
        """

        :param climb_gradient_35: The expected climb gradient at 35 ft after TO
        :param type: Indicates which speed spread is of interest (1 or 2 or 3)
            1: Speed spread from rotation to lift - off(normal rotation)(ft / s)
            2: Speed spread from lift-off to 35 ft(ft / s)
            3: Speed spread from lift-off to 15 ft(ft / s)

        :return: Speed spread
        """

        # self.twdv = #(Climb gradient at 35 ft – θ)(based on V2, gear up and AF = 0)
        # self.dvrvl = #Speed spread from rotation to lift - off(normal rotation)(ft / s)
        # self.dvlo35 = #Speed spread from lift-off to 35 ft(ft / s)
        # self.dvlo15 = #Speed spread from lift-off to 15 ft(ft / s)

        if type == 1:
            return np.interp(climb_gradient_35, self.twdv, self.dvrvl)
        elif type == 2:
            return np.interp(climb_gradient_35, self.twdv, self.dvlo35)
        elif type == 3:
            return np.interp(climb_gradient_35, self.twdv, self.dvlo15)

    def time_spread_vs_climb_gradient(self, climb_gradient_35, type):
        """

        :param climb_gradient_35: The expected climb gradient at 35 ft after TO
        :param type: Indicates which speed spread is of interest (1 or 2 or 3)
            1: Time between Vr and Vlo
            2: Time between Vlo and 35 ft
            3: Time between Vlo and 15 ft

        :return: Speed spread
        """

        # self.twdt = #(Climb gradient at 35 ft – θ)(based on V2, gear up and AF = 0)
        # self.dtvrvl = #Time between Vr and Vlo
        # self.dtvlo35 = #Time between Vlo and 35 ft
        # self.dtvlo15 = #Time between Vlo and 15 ft

        if type == 1:
            return np.interp(climb_gradient_35, self.twdt, self.dtvrvl)
        elif type == 2:
            return np.interp(climb_gradient_35, self.twdt, self.dtvlo35)
        elif type == 3:
            return np.interp(climb_gradient_35, self.twdt, self.dtvlo15)

    def NZ_buffet(self, mach, CL, cg):
        """
        Load factor at buffet start. CL buffet is given for CG at 9% MAC, needs to be adjusted!

        :param mach: Mach number
        :param CL: Current Actual Lift coefficient
        :param cg: Current cg position (% MAC)

        :return NZ_buffet: Load factor buffet start
        """

        CL_buffet = self.CL_at_buffet_vs_mach(mach) / (1 + self.mac / self.lt * (self.FWDCG - cg))
        return CL_buffet / CL

    def lift_curve_aoa(self, aoa, flap_angle, cg_act=None, gear_up=True):
        """
        BASED ON A CG POSITION OF 9% MAC.

        :param aoa: Angle of attack [°]
        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :param cg_act: cg en %
        :param gear_up: Whether the landing gear is up or not. True means up (not deployed)
        :return: Lift coefficient (CL) at given config.
        """

        if gear_up not in [False, True, 0, 1]:
            raise Exception('arg <gear_up> must be a boolean value.')

        if flap_angle == 0:
            CL_0 = 0.05
        elif flap_angle == 20:
            CL_0 = 0.25
        elif flap_angle == 45:
            CL_0 = 0.55
        else:
            raise Exception('arg <flap_angle> can only be 0, 20 or 45 degrees.')

        if not gear_up:
            CL_0 -= 0.05

        CL_fwd = CL_0 + 0.10 * aoa

        if cg_act is not None:
            return CL_fwd / (1 + self.mac / self.lt * (self.FWDCG - cg_act))
        else:
            return CL_fwd

    def get_aoa(self, CL, flap_angle, cg_act, gear_up: bool):
        """
        BASED ON A CG POSITION OF 9% MAC.

        :param CL: Angle of attack [°]
        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :param cg_act: cg en %
        :param gear_up: Whether the landing gear is up or not. True means up (not deployed)
        :return: Lift coefficient (CL) at given config.
        """

        if gear_up not in [False, True, 0, 1]:
            raise Exception('arg <gear_up> must be a boolean value.')

        if flap_angle == 0:
            CL_0 = 0.05
        elif flap_angle == 20:
            CL_0 = 0.25
        elif flap_angle == 45:
            CL_0 = 0.55
        else:
            raise Exception('arg <flap_angle> can only be 0, 20 or 45 degrees.')

        if not gear_up:
            CL_0 -= 0.05

        CL_fwd = CL * (1 + self.mac / self.lt * (self.FWDCG - cg_act))

        return (CL_fwd - CL_0) / 0.1

    def drag_curve_aoa(self, aoa, flap_angle, gear_up: bool):
        """
        DRAG DATA IS VALID FOR ALL CG LOCATIONS AND FOR ALL REYNOLDS NUMBERS
        LOW SPEED DRAG POLARS   – LANDING GEAR UP, ALL ENGINES OPERATING

        :param aoa: Angle of attack [°]
        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :param gear_up: Weither the landing gear is up or not. True means up (not deployed)
        :return: Drag coefficient (CL) at given config.
        """

        CL = self.lift_curve_aoa(aoa, flap_angle, gear_up)

        CD = self.CD_profile(flap_angle) + self.CD_induced(CL, flap_angle)

        return CD

    def CD_profile(self, flap_angle):
        """
        PROFILE DRAG COEFFICIENT

        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :return:
        """

        if flap_angle == 0:
            Cdp = 0.0206
        elif flap_angle == 20:
            Cdp = 0.0465
        elif flap_angle == 45:
            Cdp = 0.1386
        else:
            raise Exception('arg <flap_angle> can only be 0, 20 or 45 degrees.')

        return Cdp

    def CD_induced(self, CL, flap_angle):
        """
        INDUCED DRAG COEFFICIENT

        :param CL: Lift coefficient
        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :return:
        """

        K = self._induced_drag_efficiency_factor(flap_angle)

        return K * CL ** 2

    def d_CD_CTL_OEI(self, q, thrust, air=True):
        """
        ENGINE-OUT DRAG INCREMENT (in the air)

        :param q: Dynamic pressure [lbs/ft^2]
        :param thrust: Thrust from operating engine [lbs]
        :param air: Whether the aircraft is in the air or on the ground (True means in the air)
        :return:
        """

        # Engine-Out Control Drag
        if air:
            CT = thrust / (q * self.S)
            d_CD_CTL = 0.10 * CT ** 2
        else:
            d_CD_CTL = 0.0020

        return d_CD_CTL

    def CD_comp(self, mach, CL):
        """
        DRAG INCREMENT DUE TO COMPRESSIBILITY ΔCD_COMP –  APPLIES TO FLAP 0 ONLY

        :param mach: Mach number
        :param CL: Lift coefficient
        :return: ΔCD_COMP
        """
        if 0 <= mach <= 0.6:
            return 0
        elif 0.6 < mach <= 0.78:
            return (0.0508 - 0.1748 * mach + 0.1504 * mach ** 2) * CL ** 2
        elif 0.78 < mach <= 0.85:
            return (-99.3434 + 380.888 * mach - 486.8 * mach ** 2 + 207.408 * mach ** 3) * CL ** 2
        else:
            raise Exception("Mach number is beyond defined limits [0.00, 0.85]")

    def d_CD_turn(self, CL, flap_angle=0, phi=None, Nz=None):
        """
        Drag due to banking.

        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :param phi: Banking angle [°]
        :param Nz: Load Factor
        :return:
        """
        if phi == 0 or Nz == 1:
            return 0

        if phi is None and Nz is not None:
            pass
        elif phi is not None and Nz is None:
            Nz = 1 / np.cos(np.radians(phi))
        else:
            raise Exception('Ambiguous or erroneous function arguments.')

        K = self._induced_drag_efficiency_factor(flap_angle)

        return K * CL ** 2 * (Nz - 1)

    def aero_coefficient_taxi(self, flap_angle=0, spoilers=False):
        """
        Drag due to banking.

        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :param spoilers: Spoilers extension, False = not extended
        :return: cl, cd
        """

        if flap_angle == 0:
            if spoilers:
                cl = -0.0870
                cd = 0.0905
            else:
                cl = 0.2630
                cd = 0.0413
        if flap_angle == 20:
            if spoilers:
                cl = 0.2090
                cd = 0.1171
            else:
                cl = 0.8290
                cd = 0.0750
        if flap_angle == 45:
            if spoilers:
                cl = 0.4110
                cd = 0.1753
            else:
                cl = 1.361
                cd = 0.1593

        return cl, cd

    def get_Nz_stall_warning(self, p, mach, flap_angle, weight, cg, gear_up, return_phi=False):

        aoa_SW = self.fuse_AOA_SW[flap_angle]

        CL_sw = self.get_lift_coefficient(Nz=1, aoa=aoa_SW, flap_angle=flap_angle, cg=cg, gear_up=gear_up)
        q = get_dynamic_pressure(p, mach=mach)
        L = q * self.S * CL_sw
        NZ_sw = L / weight

        if not return_phi:
            return NZ_sw
        else:
            return self.get_phi(NZ_sw)

    def get_ground_coefficients(self, flaps, spoilers=False):
        """
        LIFT AND DRAG COEFFICIENTS IN TAXI ATTITUDE (WITH LANDING GEAR DOWN, IN GROUND EFFECT, VALID FOR ALL CG LOCATIONS)
        :param flaps:
        :param spoilers:
        :return:
        """
        if flaps == 0:
            if spoilers:
                CL = -0.0870
                CD = 0.0905
            else:
                CL = 0.2630
                CD = 0.0413
        elif flaps == 20:
            if spoilers:
                CL = 0.2090
                CD = 0.1171
            else:
                CL = 0.8290
                CD = 0.0750
        elif flaps == 40:
            if spoilers:
                CL = 0.4110
                CD = 0.1753
            else:
                CL = 1.361
                CD = 0.1593
        else:
            CL = None
            CD = None

        return CL, CD
    def get_lift_coefficient(self, **kwargs):
        """

        :param kwargs:

            Several parameters possible:

            ->  Keywords ['Nz', 'weight', 'q', 'S_ref'] Can be used to calculate aircraft lift coefficient
                based on aircraft weight.

            ->  Keywords ['Nz', 'aoa', 'flap_angle', 'cg', 'gear_up'] Can be used to calculate aircraft lift
                coefficient based on aircraft AoA, flap, CG and LG configuration.

            Keyword 'Nz' is optional (assumed to be equal to 1). Keyword 'phi' can also be given to calculate the
            load factor.

        :Keyword Arguments:
        * *Nz* --
          Extra stuff
        * *weight* --
          Additional content

        :return:
        """

        if 'Nz' in kwargs:
            Nz = kwargs['Nz']
        elif 'phi' in kwargs:
            Nz = 1 / np.cos(np.radians(kwargs['phi']))
        else:
            Nz = 1

        if all(key in kwargs for key in ['aoa', 'flap_angle', 'gear_up']):
            CL = Nz * self.lift_curve_aoa(aoa=kwargs['aoa'],
                                          flap_angle=kwargs['flap_angle'],
                                          cg_act=kwargs['cg'],
                                          gear_up=kwargs['gear_up'])

        elif all(key in kwargs for key in ['weight', 'q', 'S_ref']):
            CL = Nz * kwargs['weight'] / (kwargs['q'] * kwargs['S_ref'])
            # CL = CL/(1+self.mac/self.lt *(self.FWDCG-kwargs['cg']))
        else:
            raise Exception('Ambiguous or erroneous function arguments.')

        # if all(key in kwargs for key in ['cg']):
        #    CL = CL#*(1+self.mac/self.lt *(self.FWDCG-kwargs['cg']))

        return CL

    def get_lift_force(self, P, T, Weight, mach, **kwargs):
        # if 'Nz' in kwargs:
        #     Nz = kwargs['Nz']
        # elif 'phi' in kwargs:
        #     Nz = 1 / np.cos(np.radians(kwargs['phi']))
        # else:
        #     Nz = 1

        if 'S_ref' in kwargs:
            Sref = kwargs['S_ref']
        else:
            Sref = self.S

        if 'q' in kwargs:
            q = kwargs['q']
        else:
            q = get_dynamic_pressure(P, T, mach=mach)
            kwargs['q'] = q

        CL = self.get_lift_coefficient(**kwargs)
        Lift = q * CL * Sref

        return Lift

    def get_CL_max(self, flap_angle, gear_up):

        if gear_up:
            if flap_angle == 0:
                return self.CLMAX_F00_GU
            elif flap_angle == 20:
                return self.CLMAX_F20_GU
            elif flap_angle == 45:
                raise Exception('Configuration not defined')
        else:
            if flap_angle == 0:
                return self.CLMAX_F00_GD
            elif flap_angle == 20:
                return self.CLMAX_F20_GD
            elif flap_angle == 45:
                return self.CLMAX_F45_GD

    def get_drag_coefficient(self, CL, flap_angle, mach, **kwargs):
        """

        :param CL:
        :param flap_angle:
        :param mach:
        :param kwargs:

            Required Keywords:
                -> CL
                -> flap_angle
                -> mach

            Optional Keywords:
                -> OEI: Account for OEI drag contributions if True, else False
                -> LDG: Set True if landing gear deployed
                -> Nz: Load factor
                -> phi: Banking angle to calculate load factor, in degrees
                -> q
                -> thrust

        :return:
        """

        if 'OEI' in kwargs:
            OEI = kwargs['OEI']
        else:
            OEI = 0

        if 'LDG' in kwargs:
            LDG = kwargs['LDG']
        else:
            LDG = 0

        if 'Nz' in kwargs:
            Nz = kwargs['Nz']
        elif 'phi' in kwargs:
            Nz = 1 / np.cos(np.radians(kwargs['phi']))
        else:
            Nz = 1

        CDp = self.CD_profile(flap_angle)
        CDi = self.CD_induced(CL, flap_angle)
        CDi += self.d_CD_turn(CL, flap_angle, Nz)
        CDcomp = self.CD_comp(mach, CL)

        d_CD_CTL = 0
        d_CD_WM = 0

        if OEI:
            if 'q' in kwargs and 'thrust' in kwargs:
                d_CD_CTL = self.d_CD_CTL_OEI(kwargs['q'], kwargs['thrust'], air=True)
                d_CD_WM = self.d_CD_WM
            else:
                warnings.warn("Missing 'q' and 'thrust' keywords in function arguments! Defaulting to d_CD_OEI = 0.0")

        if LDG:
            CDp += self.d_CD_LG

        CDtotal = CDp + CDi + CDcomp + d_CD_CTL + d_CD_WM

        drag = {'CDp': CDp,
                'CDi': CDi,
                'CDcomp': CDcomp,
                'CDctl': d_CD_CTL,
                'CDwm': d_CD_WM,
                'CDtot': CDtotal}

        return drag

    def get_drag_force(self, P, T, Weight, flap_angle, mach, Sref=None, **kwargs):
        if 'Nz' in kwargs:
            Nz = kwargs['Nz']
        elif 'phi' in kwargs:
            Nz = 1 / np.cos(np.radians(kwargs['phi']))
        else:
            Nz = 1

        if Sref is None:
            Sref = self.S

        q = get_dynamic_pressure(P, T, mach=mach)
        kwargs['q'] = q
        CL = self.get_lift_coefficient(Nz=Nz, weight=Weight, q=q, S_ref=Sref)
        CD = self.get_drag_coefficient(CL, flap_angle, mach, **kwargs)['CDtot']
        Drag = q * CD * self.S

        return Drag

    def get_minimum_drag_speed(self, W, rho):
        """
        :param:
        :return:
        """
        CDp = self.CD_profile(0)
        K = self._induced_drag_efficiency_factor(0)

        CL = np.sqrt(CDp / K)

        V_md = np.sqrt(W / (0.5 * rho * CL * self.S)) / 1.68781

        return V_md

    def get_max_range_cruise_speed(self, W, rho):
        """
        :param:
        :return:
        """
        CDp = self.CD_profile(0)
        K = self._induced_drag_efficiency_factor(0)

        CL = np.sqrt(CDp / (3 * K))

        V_MRC = np.sqrt(W / (0.5 * rho * CL * self.S)) / 1.68781

        return V_MRC

    def get_long_range_cruise_mach(self, W, T, rho, P, hp):
        """
        :param:
        :return:
        """
        Mach_v = np.linspace(0.8, 0.45, 100)
        SAR_v = np.zeros(len(Mach_v))
        for i in range(0, len(Mach_v)):
            q = get_dynamic_pressure(P, T=T, mach=Mach_v[i])
            L = W
            CL = self.get_lift_coefficient(Nz=1, weight=W, q=q, S_ref=self.S)
            CD = self.get_drag_coefficient(CL=CL, flap_angle=0, mach=Mach_v[i])['CDtot']
            D = CD * q * self.S
            SFC = self.get_SFC(hp)

            SAR_v[i] = ((a0 * np.sqrt(T / T_0)) / SFC) * (Mach_v[i] * L / D) * (1 / W)

        SAR_MRC = max(SAR_v)
        SAR_LRC = 0.99 * SAR_MRC
        SAR_v_cut = SAR_v[:np.where(SAR_v == SAR_MRC)[0][0]]
        Mach_v_cut = Mach_v[:np.where(SAR_v == SAR_MRC)[0][0]]

        Mach_LRC = np.interp(SAR_LRC, SAR_v_cut, Mach_v_cut)

        return Mach_LRC

    def get_thrust(self, RM, PALT, M, T, n_engines=None):

        """
        :param RM:
        :param PALT:
        :param n_engines:

        :return: Thrust
        """

        if n_engines is None:
            n_engines = self.n_engines

        dISA = get_delta_ISA(PALT, T)

        RM = RM.upper()
        if RM == 'MTO' or RM == 'GA' or RM == 'MCT':
            # MTOFN = Maximum Take-Off (MTO) thrust per engine (lb) (flat rated to ISA+15) (valid up to ISA+15)
            # Note: Above ISA+15, MTOFN reduces by 1 % per degree C
            T_OE = 8775 - 0.1915 * PALT - (8505 - 0.195 * PALT) * M

            if dISA > 15:
                T_OE *= (1 - 1 / 100 * (dISA - 15))

            if RM == 'MCT':
                T_OE *= 0.90

        elif RM == 'MCL' or RM == 'MCR':
            # MCLFN = Maximum Climb (MCL) thrust per engine (lb) (flat rated to ISA+10) (valid up to ISA+10)
            # Note: Above ISA+10, MCLFN reduces by 1 % per degree C
            T_OE = 5690 - 0.0968 * PALT - (1813 - 0.0333 * PALT) * M

            if dISA > 10:
                T_OE *= (1 - 1 / 100 * (dISA - 10))

            if RM == 'MCR':
                T_OE *= 0.98

        elif RM == 'ID' or RM == 'IDLE':
            T_OE = 600 - 1000 * M  # Idle thrust per engine (lb) (independent of altitude & temperature)
        elif RM == 'MXR':
            T_OE = -1300 - 12000 * M  # Max. Reverse thrust per engine (lb) (indep. of altitude & temperature)
        elif RM == 'IDR':
            T_OE = -160 - 3700 * M  # Idle Reverse thrust per engine (lb) (indep. of altitude & temperature)
        else:
            raise Exception('Unexpected engine rating')

        return T_OE * n_engines

    def get_SFC(self, PALT):
        """
        :param PALT:

        :return: SFC
        """
        # Specific Fuel Consumption (lb/hr of fuel per lb of thrust per engine)
        # Note: If thrust is below 0 lb, use T = 600 lb/engine for fuel flow calculation.
        SFC = 0.58 + (0.035 * PALT / 10000)

        return SFC

    def get_NZ(self, L, W):
        """
        :param L:
        :param W:

        :return: Nz
        """

        return L / W

    def get_NZ_sw(self, W, flap_angle, acg_act, gear_up, p, M=None):
        aoa_SW = self.fuse_AOA_SW[flap_angle]
        CL_sw = self.lift_curve_aoa(aoa_SW, flap_angle, acg_act, gear_up)
        qp = get_dynamic_pressure(p, mach=M)
        L_stall = CL_sw * qp * self.S
        NZ_sw = L_stall / W
        if NZ_sw:
            return NZ_sw
        else:
            return self.get_phi(NZ_sw)

    def get_phi(self, NZ):
        """
        :param NZ:

        :return: phi: Bank angle
        """

        return np.degrees(np.arccos(1 / NZ))

    def get_cg_from_CL(self, CL_act, AoA, flap_angle, gear_up):

        CL_fwd = self.lift_curve_aoa(AoA, flap_angle, gear_up=gear_up)

        # cg_act = (1-CL_fwd/CL_act)*lt/mac

        return (1 - CL_fwd / CL_act) * self.lt / self.mac + self.FWDCG

    def get_mach_from_VSR(self, p, weight, vsr_mult, Nz, flap_angle, gear_up):

        delta = p / P_0
        CL = self.get_CL_max(flap_angle, gear_up=gear_up) / vsr_mult ** 2
        mach = np.sqrt((weight * Nz / delta) / 1481.3 / CL / self.S)
        return mach

    def get_fuel_burn_rate(self, hp, thrust):
        """
        Fuel burn rate in lbs/s.

        :param p:
        :param thrust:
        :return:
        """

        SFC = self.get_SFC(hp)
        if thrust <= 0:
            thrust = 1200

        return SFC * thrust / 60 / 60

    def mission_performance_climb_descent(self, hpi, hpf, dISA, Vwind, V_CAS_cst, Mach_cst, Wi, RM, int_hp_resol=50,
                                          OEI_tag=False, ROC_min=300):

        cg = 0.25

        if OEI_tag:
            n_engines = 1
        else:
            n_engines = 2

        # identify the situation
        isAscent = True if hpi < hpf else False

        # in case of descent, integration interval is negative
        if not isAscent:
            int_hp_resol *= -1

        # Calculate the integration intervals
        int_steps = np.arange(start=hpi, stop=hpf + int_hp_resol, step=int_hp_resol)

        # Last altitude needs to be adjusted
        if int_steps[-1] != hpf:
            int_steps[-1] = hpf

        # This accounts for discontinuities in the ROC based on regulatory restrictions under 10 000 ft
        if hpi <= 10000 <= hpf:  # Account for a discontinuity at 10000 ft
            if 10000 not in int_steps:
                int_steps = np.insert(int_steps, min(np.where(10000 < int_steps))[0], [10000, 10000])
            else:
                int_steps = np.insert(int_steps, min(np.where(10000 == int_steps))[0], [10000])
        elif hpi >= 10000 >= hpf and 10000 not in int_steps:  # Account for a discontinuity at 10000 ft
            if 10000 not in int_steps:
                int_steps = np.insert(int_steps, max(np.where(10000 < int_steps))[0], [10000, 10000])
            else:
                int_steps = np.insert(int_steps, max(np.where(10000 == int_steps))[0], [10000])

        # Calculate transition altitude
        Ps_tr = get_impact_pressure(P_0, V_CAS_cst / a0) / ((1 + 0.2 * Mach_cst ** 2) ** 3.5 - 1)
        hp_tr = get_pressure_altitude(Ps_tr)
        hp_tr = round(hp_tr / 100) * 100

        # This accounts for discontinuities in the ROC based on transition from cst CAS climb to cst Mach climb
        if hpi < hp_tr < hpf and hp_tr not in int_steps:  # Account for a discontinuity at 10000 ft
            int_steps = np.insert(int_steps, min(np.where(hp_tr < int_steps))[0], hp_tr)
        elif hpi > hp_tr > hpf and hp_tr not in int_steps:  # Account for a discontinuity at 10000 ft
            int_steps = np.insert(int_steps, max(np.where(hp_tr < int_steps))[0], hp_tr)

        intermediate_data = {}

        # Acc, AccTASAvg, AccTime, AccDist, AccFuel, AccWfi = [0, 0, 0, 0, 0, 0]
        W1 = Wi
        W2 = W1
        fuel_burned_idx = 1 / 100000
        Wfuel = 0
        t_total = 0
        dist_total = 0

        for idx, hp1 in enumerate(int_steps[:-1]):

            while 1 - abs((W1 - W2) / fuel_burned_idx) > 1 / 100:

                W2 = W1 - fuel_burned_idx
                hp2 = int_steps[idx + 1]
                hp_moy = hp2 / 2 + hp1 / 2

                p = pressure_from_alt(hp_moy)
                T_ISA = temp_from_alt(hp_moy)
                T_case = T_ISA + dISA

                delta_hp = abs(hp2 - hp1)
                delta_hg = delta_hp * (T_case / T_ISA)
                W_avg = (W2 + W1) / 2

                if hp1 == 10000:  # Acceleration Segment
                    V1 = 250
                    V2 = V_CAS_cst
                    V_CAS_AVG = (V1 + V2) / 2

                    Mach1 = get_mach_from_calibrated_airspeed(p, V1)
                    Mach2 = get_mach_from_calibrated_airspeed(p, V2)
                    Mach = get_mach_from_calibrated_airspeed(p, V_CAS_AVG)
                    q = get_dynamic_pressure(p, T_case, mach=Mach)

                    Thrust = self.get_thrust(RM, hp_moy, Mach, T_case, n_engines=n_engines)
                    CL = self.get_lift_coefficient(Nz=1, weight=W_avg,
                                                   q=get_dynamic_pressure(p, T_case, mach=Mach),
                                                   S_ref=self.S)
                    CD = self.get_drag_coefficient(CL, 0, Mach, LDG=0, NZ=1, OEI=OEI_tag, q=q,
                                                   thrust=Thrust)['CDtot']
                    D = q * CD * self.S

                    Acc = (Thrust - D) / W_avg * g
                    AccTASAvg = get_true_airspeed(p, Mach, temp=T_case, knots=False)
                    AccTime = abs((get_true_airspeed(p, Mach2, temp=T_case, knots=False) - get_true_airspeed(p, Mach1,
                                                                                                             temp=T_case,
                                                                                                             knots=False)) / Acc)
                    AccDist = AccTime * (AccTASAvg + knots2fps(Vwind))
                    fuel_burn_rate = self.get_fuel_burn_rate(hp_moy, Thrust)
                    AccFuel = fuel_burn_rate * AccTime
                    AccWfi = W1

                    fuel_burned_idx = AccFuel
                    d_time_idx = AccTime
                    d_dist_idx = AccDist

                else:
                    if hp1 < hp_tr:
                        if hp1 < 10000:
                            V_CAS = 250
                        else:
                            V_CAS = V_CAS_cst

                        Mach = get_mach_from_calibrated_airspeed(p, V_CAS)
                        RV = 'CAS'

                    else:
                        Mach = Mach_cst
                        RV = 'Mach'

                    V_TAS = get_true_airspeed(p, Mach, temp=T_case, knots=False)
                    q = get_dynamic_pressure(p, T_case, mach=Mach)
                    CL = self.get_lift_coefficient(Nz=1, weight=W_avg, q=q, S_ref=self.S)
                    Thrust = self.get_thrust(RM, hp_moy, Mach, T=T_case, n_engines=n_engines)
                    D = q * self.get_drag_coefficient(CL, 0, Mach, LDG=0, NZ=1, OEI=OEI_tag, q=q,
                                                      thrust=Thrust)['CDtot'] * self.S
                    AF = get_AF(RV, hp_moy, Mach, T_case, T_ISA)

                    if isAscent:
                        ROC = get_ROC(V_TAS, Thrust, D, W_avg, AF)
                    else:
                        ROC = get_ROC(V_TAS, -Thrust, -D, W_avg, AF)

                    d_time_idx = delta_hg / (ROC / 60)
                    d_dist_idx = (V_TAS + knots2fps(Vwind)) * d_time_idx
                    fuel_burn_rate = self.get_fuel_burn_rate(hp_moy, Thrust)
                    fuel_burned_idx = fuel_burn_rate * d_time_idx

            W1 = W2

            if isAscent and ROC / (T_case / T_ISA) < ROC_min:
                break

            Wfuel += fuel_burned_idx
            t_total += d_time_idx
            dist_total += d_dist_idx

        Wf = Wi - Wfuel
        hpf = hp2

        return t_total, dist_total / 6076, Wfuel, hp_tr, hpf, Acc, AccTASAvg, AccTime, AccDist / 6076, AccFuel, AccWfi

    def cruise(self, hp, dISA, Vwind, V, V_type, W_cri, dW_fuel_cruise, OEI_tag=False):
        """
        Paramêtre de cruise.

        :param hp: Altitude(ft)
        :param dISA: Delta ISA (C)
        :param Vwind: Vitesse vend (kts)
        :param V: Aircraft speed [-, mach number, -, -]
        :param V_type: ["Vmd", "Mach", "MRC", "LRC"] [Min drag, mach number cste, long range cruise, max cruise speed]
        :param W: weight (lb)
        :return KCAS: Vitesse CAS en kts
        :return V_g: vitesse sol kts
        :return SAR: Distance aerienne fanchissable specifique (nm/lb)
        :return SR: Distance fanchissable specifique (nm/lb)
        :return Wf: Debit massique total de carburant (lb/hr)
        """
        P, rho, T = get_atmos_from_dISA(hp, dISA)
        W_avg = W_cri - dW_fuel_cruise/2
        if V_type == "Vmd":
            V_md = self.get_minimum_drag_speed(W_avg, rho)
            V = V_md
            Mach = get_mach(v=V_md, temp=T)
        elif V_type == "Mach":
            Mach = V
            V = get_true_airspeed(p=P, mach=Mach, temp=T)
        elif V_type == "MRC":
            V_MRC = self.get_max_range_cruise_speed(W_avg, rho)
            V = V_MRC
            Mach = get_mach(v=V_MRC, temp=T)
        elif V_type == "LRC":
            Mach = self.get_long_range_cruise_mach(W_avg, T, rho, P, hp)
            V = get_true_airspeed(p=P, mach=Mach, temp=T)


        q = get_dynamic_pressure(P, T=T, mach=Mach)

        L = W_avg
        CL = self.get_lift_coefficient(Nz=1, weight=W_avg, q=q, S_ref=self.S)
        CD = self.get_drag_coefficient(CL=CL, flap_angle=0, mach=Mach)['CDtot']
        D = CD * q * self.S
        thrust = D
        Wf = self.get_fuel_burn_rate(hp, thrust)
        SFC = self.get_SFC(hp)

        SAR = ((a0 * np.sqrt(T / T_0)) / SFC) * (Mach * L / D) * (1 / W_avg)
        V_g = V + Vwind
        SR = SAR * (V_g / V)

        KCAS = get_calibrated_airspeed(p=P, mach=Mach)

        return KCAS, V_g, Mach, SAR, SR, Wf

    def mission(self, dW_cargo, hpi, hpf, RM_cl, V_CAS_cst, Mach_cst, hp_cruise, RM_d, dISA, Vwind, Wi_fuel=15000,
                dW_fuel_to_taxi=200, dW_fuel_to=250, dW_fuel_ldg=200, dW_fuel_ldg_taxi=100, W_fuel_reserve=2000):
        '''

        :param Wi_cargo:
        :param hpi:
        :param hpf:
        :param RM_cl:
        :param V_CAS_cst:
        :param Mach_cst:
        :param hp_cruise:
        :param dISA:
        :param Vwind:
        :param Wi_fuel:
        :param dW_fuel_to_taxi:
        :param dW_fuel_to:
        :param dW_fuel_ldg:
        :param dW_fuel_ldg_taxi:
        :param W_fuel_reserve:
        :return:
        '''

        ZFW = self.OWE + dW_cargo  # Zero Fuel Weight - constant unless passengers are lost (or found?) along the way
        LW = ZFW + W_fuel_reserve
        RW = ZFW + Wi_fuel  # Ramp weight - weight before taxi

        # Taxi
        TOW = RW - dW_fuel_to_taxi  # Takeoff Weight (after taxi)

        if LW > self.MLW:
            raise Exception('LW depasse MLW.')
        elif RW > self.MRW:
            raise Exception('RW depasse MRW.')
        elif TOW > self.MTOW:
            raise Exception('TOW depasse MTOW.')

        # Takeoff
        W_cli = TOW - dW_fuel_to  # Initial climb weight (after initial climb to 1500ft)
        hpi_cl = hpi + 1500  # Initial climb altitude

        # Climb
        Climb = self.mission_performance_climb_descent(hpi_cl, hp_cruise, dISA, Vwind, V_CAS_cst, Mach_cst, W_cli,
                                                       RM_cl)
        dist_climb = Climb[1]
        dW_fuel_climb = Climb[2]
        hp_cruise = round(Climb[4]/1000)*1000

        hpi_d = hp_cruise  # First estimate of initial descent altitude (same as cruise altitude)
        hpf_d = hpf + 1500  # Altitude at final approach (end of descent)

        W_cri = W_cli - dW_fuel_climb  # Initial cruise weight
        W_di = 0  # Declare Variable
        W_crf = W_cri  # First estimate of initial descent altitude (same as cruise altitude)
        dW_fuel_cruise = 0  # First estimate of initial descent Weight
        while abs(W_di - W_crf) / W_crf > 1 / 10000:
            # Descent
            W_di = W_cri - dW_fuel_cruise  # Estimate of initial descent Weight
            Descent = self.mission_performance_climb_descent(hpi_d, hpf_d, dISA, Vwind, V_CAS_cst, Mach_cst, W_di, RM_d)
            dist_descent = Descent[1]
            dW_fuel_descent = Descent[2]

            # Cruise
            W_crf = W_fuel_reserve + dW_fuel_ldg_taxi + dW_fuel_ldg + dW_fuel_descent + ZFW  # Estimate of Weight after cruise
            dW_fuel_cruise = W_cri - W_crf
            Cruise = self.cruise(hp_cruise, dISA, Vwind, Mach_cst, 'Mach', W_cri, dW_fuel_cruise)

            dist_cruise = Cruise[3] * dW_fuel_cruise

            # hpi_d = get_pressure_altitude(P_0*W_crf/(W_cri/pressure_from_alt(hp_cruise, ratio=True)))
            # if hpi_d > 41000:
            #     hpi_d = 41000

        dist_tot = Climb[1] + Descent[1] + dist_cruise
        dWf_tot = RW - (W_crf - dW_fuel_descent - dW_fuel_ldg - dW_fuel_ldg_taxi)

        return dist_climb, hp_cruise, dist_cruise, dist_descent, dist_tot, dWf_tot


    def takeoff_run_velocities(self, W, Hp, dISA, theta=0, OEI_tag=False):
        """
        Calcul des vitesses pertinentes lors du décollage.

        :param W: Poids de l'avion [lbs]
        :param Hp: altitude pression [ft]
        :param dISA: déviation ISA [²C ou K]

        :return: V1Min, V1Max, VR, V2, V_LO_OEI, V_LO_AEO, V_35_AEO
        """
        # Hypotheses
        flap_angle = 20
        RM = 'MTO'

        # Conditions Atmospheriques
        P = pressure_from_alt(Hp)
        T = temp_from_alt(Hp) + dISA
        delta = P / P_0

        # V_SR
        CL = self.get_CL_max(flap_angle, gear_up=True) / 1 ** 2
        mach_Vsr = np.sqrt((W * 1 / delta) / 1481.3 / CL / self.S)
        V_SR = get_calibrated_airspeed(P,mach_Vsr)

        # Autres Vitesses
        V1_MCG = self.V_1mcg # CAS
        V_MCA = self.V_mca # CAS

        # V1, VR et V2 Minimums
        V1Min = knots2fps(CAS2TAS(Hp,dISA, V1_MCG))
        VRMin = knots2fps(CAS2TAS(Hp,dISA, max(V1_MCG, 1.05*V_MCA)))
        V2Min = knots2fps(CAS2TAS(Hp,dISA, max(1.13*V_SR,1.1*V_MCA)))

        # Vitesses OEI et VR
        n_engines = 1
        OEI_tag = True
        V2 = V2Min
        Mach = V2/get_SOS(T,knots=False)
        q = get_dynamic_pressure(P, T, mach=Mach)
        Thrust = self.get_thrust(RM, Hp, Mach, T, n_engines=n_engines)
        Drag = self.get_drag_force(P, T, W, flap_angle, Mach, LDG=0, NZ=1, thrust=Thrust, OEI=OEI_tag)
        grad35ft = (Thrust-Drag)/W-theta
        if grad35ft < 2.4/100:
            print('WARNING! Climb gradient at 35ft is less than 2.4%. Calculations aborted.\n')
            return None
        dV_LO_35_OEI = self.speed_spread_vs_climb_gradient(grad35ft, 2)
        V_LO_OEI = V2 - dV_LO_35_OEI
        dV_ROT_LO_OEI = self.speed_spread_vs_climb_gradient(grad35ft, 1)
        VR = V_LO_OEI - dV_ROT_LO_OEI

        #Time OEI
        dt_VLO_VR_OEI = self.time_spread_vs_climb_gradient(grad35ft, 1)
        dt_V35_VLO_OEI = self.time_spread_vs_climb_gradient(grad35ft, 2)
        # Dist AEO
        D_VLO_VR_OEI = (VR + dV_ROT_LO_OEI / 2) * dt_VLO_VR_OEI
        D_V35_VLO_OEI = (V_LO_OEI + dV_LO_35_OEI / 2) * dt_V35_VLO_OEI


        # Correction VR, V2 et V_LO_OEI si VR < VRmin
        if VR < VRMin:
            VR = VRMin
            V_LO_OEI = VR + dV_ROT_LO_OEI
            V2 = V_LO_OEI + dV_LO_35_OEI

        # Vitesses AEO
        n_engines = 2
        OEI_tag = False
        Mach = V2 / get_SOS(T, knots=False)
        q = get_dynamic_pressure(P, T, mach=Mach)
        Thrust = self.get_thrust(RM, Hp, Mach, T, n_engines=n_engines)
        Drag = self.get_drag_force(P, T, W, flap_angle, Mach, LDG=0, NZ=1, thrust=Thrust, OEI=OEI_tag)
        grad35ft = (Thrust - Drag) / W - theta
        dV_ROT_LO_AEO = self.speed_spread_vs_climb_gradient(grad35ft, 1)
        V_LO_AEO = VR + dV_ROT_LO_AEO
        dV_LO_35_AEO = self.speed_spread_vs_climb_gradient(grad35ft, 2)
        V_35_AEO = V_LO_AEO + dV_LO_35_AEO

        # Time AEO
        dt_VLO_VR_AEO = self.time_spread_vs_climb_gradient(grad35ft, 1)
        dt_V35_VLO_AEO = self.time_spread_vs_climb_gradient(grad35ft, 2)

        # Dist AEO
        D_VLO_VR_AEO = (VR + V_LO_AEO) / 2 * dt_VLO_VR_AEO
        D_V35_VLO_AEO = (V_LO_AEO + V_35_AEO) * dt_V35_VLO_AEO/2

        # V1_Max
        V1Max = VR


        return V1Min, V1Max, VR, V2, V_LO_OEI, V_LO_AEO, V_35_AEO, dt_VLO_VR_OEI, dt_V35_VLO_OEI, dt_VLO_VR_AEO, dt_V35_VLO_AEO, D_VLO_VR_OEI, D_V35_VLO_OEI, D_VLO_VR_AEO, D_V35_VLO_AEO

    def takeoff_run_distances(self, Hp, W, V1VR, V_wind,T=None,dISA = None, flap_angle=20, cg=20, theta=0):
        """

        :return:
            - FTOD_AEO: AEO Factored takeoff distance [ft]
            - TOD_OEI: OEI Takeoff Distance [ft]
            - ASD_AEO: AEO Acceleration-stop distance [ft]
            - BFL(?): Balanced field length [ft]
            - E_br: Braking energy [in millions of ft-lbs]
        """
        def segment(Vi, dVg, V_wind, OEI_tag, braking=False):
            if braking:
                RM = 'ID'
                mu = self.mu_dry
                spoilers = True
            else:
                RM = 'MTO'
                mu = self.rollmu
                spoilers = False

            if OEI_tag:
                n_engines = 1
            else:
                n_engines = 2

            V_IA = Vi + V_wind  # V_wind is positive when wind is facing the aircraft
            V_FA = V_IA + dVg
            V_RMS = 2 ** .5 / 2 * (V_IA ** 2 + V_FA ** 2) ** 0.5
            Mach = V_RMS / get_SOS(T, knots=False)
            q = get_dynamic_pressure(P, T, mach=Mach)
            Thrust = self.get_thrust(RM, Hp, Mach, T, n_engines=n_engines)
            CL_fwd, CD = self.get_ground_coefficients(flap_angle, spoilers=spoilers)
            CL = CL_fwd / (1 + self.mac / self.lt * (self.FWDCG - cg))
            Lift = q * self.S * CL
            Drag = q * self.S * CD
            a = g / W * (Thrust - Drag - mu * (W - Lift) - W * theta)
            dt = dVg / a
            ds = dt * (Vi+dVg/2)

            if braking:
                braking_energy = 0.4 * (W - q * self.S * CL) * ds / 1e6
            else:
                braking_energy = None

            return ds, braking_energy

        # Conditions Atmospheriques
        P = pressure_from_alt(Hp)
        if dISA == None:
            dISA = get_delta_ISA(Hp, T)
        elif T==None:
            P,Rho,T = get_atmos_from_dISA(Hp, dISA, ratio=False)
        else:
            print(f'Soit T ou dISA doivent etre fournis')

        # Hypotheses
        mu = self.rollmu

        # Velocities
        velocities = self.takeoff_run_velocities(W, Hp, dISA)
        if velocities is None:
            return [None]*5
        else:
            V1Min, V1Max, VR, V2, V_LO_OEI, V_LO_AEO, V_35_AEO = velocities [0:7]
            dt_VLO_VR_OEI, dt_V35_VLO_OEI, dt_VLO_VR_AEO, dt_V35_VLO_AEO = velocities[7:11]
            D_VLO_VR_OEI, D_V35_VLO_OEI, D_VLO_VR_AEO, D_V35_VLO_AEO = velocities[11:16]

        V1 = V1VR * VR
        if V1<V1Min:
            print('WARNING! Given V1/VR value results in V1 value lower than V1mcg. Corrections applied.\n')
            # Velocity increments
            dVR_VLO_AEO = V_LO_AEO - VR
            dVR_VLO_OEI = V_LO_OEI - VR
            dVLO_V35 = V_35_AEO - V_LO_AEO

            VR = V1Min/V1VR
            V1 = V1Min

            # New velocities at LO and 35ft
            V_LO_AEO = VR + dVR_VLO_AEO
            V_LO_OEI = VR + dVR_VLO_OEI
            V_35 = V_LO_AEO + dVLO_V35
            # Distances correction
            D_VLO_VR_AEO = (VR + V_LO_AEO) / 2 * dt_VLO_VR_AEO
            D_V35_VLO_AEO = (V_LO_AEO + V_35) / 2 * dt_V35_VLO_AEO
            D_VLO_VR_OEI = (VR + V_LO_OEI) / 2 * dt_VLO_VR_OEI
            D_V35_VLO_OEI = (V_LO_OEI + V_35) / 2 * dt_V35_VLO_OEI

        if V_wind < 0:
            V_wind *= 1.5
        else:
            V_wind *= 0.5

        # CALCULATE FTOD_AEO - AEO Factored takeoff distance [ft]
        OEI_tag = False
        d_V0_V1_AEO, __ = segment(0, V1, V_wind,OEI_tag)
        d_V1_VR_AEO, __ = segment(V1, VR-V1, V_wind, OEI_tag)
        d_VR_VLO_AEO = D_VLO_VR_AEO
        d_VL0_V35_AEO = D_V35_VLO_AEO
        TOD_AEO = d_V0_V1_AEO+d_V1_VR_AEO+d_VR_VLO_AEO+d_VL0_V35_AEO
        FTOD_AEO = 1.15 * TOD_AEO

        # CALCULATE TOD_OEI - OEI Takeoff Distance [ft]
        OEI_tag = True
        d_V1_VR_OEI, __ = segment(V1, VR - V1, V_wind, OEI_tag)
        d_VR_VLO_OEI = D_VLO_VR_OEI
        d_VL0_V35_OEI = D_V35_VLO_OEI
        TOD_OEI = d_V0_V1_AEO + d_V1_VR_OEI + d_VR_VLO_OEI + d_VL0_V35_OEI

        # CALCULATE ASD_AEO - AEO Acceleration-stop distance [ft]
        OEI_tag = False
        d_V1_V0_AEO, braking_energy = segment(V1, -V1, V_wind, OEI_tag, braking=True)
        ASD_AEO = d_V0_V1_AEO + d_V1_V0_AEO + 2*V1

        Longueur_min = max(FTOD_AEO, TOD_OEI, ASD_AEO)

        return FTOD_AEO, TOD_OEI, ASD_AEO, Longueur_min, braking_energy
