from typing import Any
import numpy as np


class aircraft():
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
        self.FWDCG = 9.0

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

        # Note: Above ISA+15, MTOFN reduces by 1 % per degree C
        self.MTOFN = 8775 - 0.1915 * PALT - (8505 - 0.195 * PALT) * M  # Maximum Take-Off (MTO) thrust per engine (lb) (flat rated to ISA+15) (valid up to ISA+15)
        self.GAFN = self.MTOFN

        # Note: Above ISA+10, MCLFN reduces by 1 % per degree C
        self.MCLFN = 5690 - 0.0968 * PALT - (1813 - 0.0333 * PALT) * M  # Maximum Climb (MCL) thrust per engine (lb) (flat rated to ISA+10) (valid up to ISA+10)

        self.MCRFN = MCLFN * 0.98  # Maximum Cruise thrust (MCR) per engine (lb) (flat rated to ISA+10)
        self.MCTFN = MTOFN * 0.90  # Maximum Continuous Thrust (MCT) per engine (lb) (flat rated to ISA+15)
        self.IDLEFN = 600 - 1000*M  # Idle thrust per engine (lb) (independent of altitude & temperature)
        self.MXRVFN = -1300 -12000*M  # Max. Reverse thrust per engine (lb) (indep. of altitude & temperature)
        self.IDRVFN = -160 -3700*M  # Idle Reverse thrust per engine (lb) (indep. of altitude & temperature)

        # Note: If thrust is below 0 lb, use T = 600 lb/engine for fuel flow calculation.
        self.SFC = 0.58 + (0.035 * PALT/ 10000)  # Specific Fuel Consumption (lb/hr of fuel per lb of thrust per engine)

        #                **********************
        #                     Speed limits
        #                **********************

        self.V_MO = 330  # Maximum CAS with no flaps deployed [kts]
        self.M_MO = 0.85  # Maximum Mach number with no flaps deployed
        self.V_FE_F20 = 210  # Maximum CAS with flaps deployed at 20° [kts]
        self.V_FE_F45 = 180  # Maximum CAS with flaps deployed at 20° [kts]
        self.V_LO = self.V_LE = 220  # Maximum CAS with LG deployed [kts]
        self.V_mca = 95  #minimum control speed in the air (take-off) (KCAS)
        self.V_mcl = 92  #minimum control speed in the air (landing) (KCAS)
        self.V_mcg = 90  #minimum control speed on the ground (take-off) (KCAS)
        self.V_1mcg = 95  #minimum v1 based on engine failure at vmcg (KCAS)

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
        self.d_CD_WM = 0.0030 # (windmilling drag coefficient) - OEI

        # *********************************************************
        #                   Speed and time increments
        #  *******************************************************
        #self.twdv = #(Climb gradient at 35 ft – θ)(based on V2, gear up and AF = 0)
        #self.dvrvl = #Speed spread from rotation to lift - off(normal rotation)(ft / s)
        #self.dvlo35 = #Speed spread from lift-off to 35 ft(ft / s)
        #self.dvlo15 = #Speed spread from lift-off to 15 ft(ft / s)

        #self.twdt = #(Climb gradient at 35 ft – θ)(based on V2, gear up and AF = 0)
        #self.dtvrvl = #Time between Vr and Vlo
        #self.dtvlo35 = #Time between Vlo and 35 ft
        #self.dtvlo15 = #Time between Vlo and 15 ft

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
        self.kemax = 16.7 #max demonstrated brake energy per brake (million ft-lb)
        self.fuse_plug_limit = 11.5 #fuse plug limit brake energy per brake (million ft-lb)
        self.mu_dry = 0.400  #Airplane braking coefficient on dry runway
        self.rollmu = 0.200  #Rolling coefficient of friction
        self.mu_wet = 0.225  #Brake mu on wet runways
        self.mu_snow = 0.200  #Brake mu on compacted snow covered runway
        self.mu_ice = 0.050  #Brake mu on ice covered runway
        self.Pm = 165  #Main tire pressure (psi)
        self.tiremax = 210  #Maximum allowable tire speed (mph)
        self.nb_brake = 4  #Number of brakes

        #  *********************************************************
        #                       T.O / R.T.O. Time Delays
        #  *********************************************************
        self.dtdwm = 0  #Time to reach windmilling drag level following an engine cut (sec)
        self.dtapr = 0  #Time to reach apr thrust (sec)
        self.dtrec = 1  #Engine failure recognition time (sec)
        self.dtidle = 0  #Time to reach idle thrust after throttles chopped (sec)
        self.tbrake = 0  #Time from V1 to brake application (sec)
        self.tpower = 0  #Time from V1 to throttle chop to idle( sec)
        self.tdump = 0  #Time from V1 to GLD fully deployed (sec)


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

    def lift_curve_aoa(self, aoa, flap_angle, gear_up: bool):
        """
        BASED ON A CG POSITION OF 9% MAC.

        :param aoa: Angle of attack [°]
        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :param gear_up: Whether the landing gear is up or not. True means up (not deployed)
        :return: Lift coefficient (CL) at given config.
        """

        if gear_up not in [False, True, 0 , 1]:
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

        return CL_0 + 0.10*aoa

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

        if flap_angle == 0:
            Cdp = 0.0206
            K = 0.0364
        elif flap_angle == 20:
            Cdp = 0.0465
            K = 0.0334
        elif flap_angle == 45:
            Cdp = 0.1386
            K = 0.0301
        else:
            raise Exception('arg <flap_angle> can only be 0, 20 or 45 degrees.')

        return Cdp + K*CL^2

    def CL_at_buffet_vs_mach(self, mach):
        """
        LIFT COEFFICIENT AT BUFFET (CL_BUFFET) AS A FUNCTION ON MACH
        DATA BASED ON 9 % MAC AND FLAP 0

        :param mach: Mach number
        :return:
        """

        x_pts = [0.2750,  0.3000,  0.3250,  0.3500,  0.3750,  0.4000,  0.4250,  0.4500,
                 0.4750,  0.5000,  0.5250,  0.5500,  0.5750,  0.6000,  0.6250,  0.6500,
                 0.6750,  0.7000,  0.7250,  0.7500,  0.7750,  0.8000,  0.8250,  0.8500,
                 0.8750,  0.9000]

        y_pts = [1.3424, 1.3199, 1.2974, 1.2667, 1.2310, 1.1930, 1.1551, 1.1191, 1.0863,
                 1.0577, 1.0337, 1.0142, 0.9989, 0.9868, 0.9764, 0.9659, 0.9530, 0.9349,
                 0.9085, 0.8698, 0.8149, 0.7391, 0.6373, 0.5039, 0.3330, 0.118]

        return np.interp(mach, x_pts, y_pts)

    def d_CD_OEI(self, q, thrust, air=True):
        """
        ENGINE-OUT CONTROL DRAG (in the air)

        :param q: Dynamic pressure [lbs/ft^2]
        :param thrust: Thrust from operating engine [lbs]
        :param air: Whether the aircraft is in the air or on the ground (True means in the air)
        :return:
        """
        #Windmilling Drag
        d_CD_WM = self.d_CD_WM

        # Engine-Out Control Drag
        if air:
            CT = thrust / (q * self.S)
            d_CD_NTL = 0.10*CT^2
        else:
            d_CD_NTL = 0.0020

        return d_CD_WM + d_CD_NTL

    def d_CD_comp(self, mach, CL):
        """
        DRAG INCREMENT DUE TO COMPRESSIBILITY ΔCD_COMP –  APPLIES TO FLAP 0 ONLY

        :param mach: Mach number
        :param CL: Lift coefficient
        :return: ΔCD_COMP
        """
        if 0 <= mach <= 0.6:
            return 0
        elif 0.6 < mach <= 0.78:
            return (0.0508 - 0.1748*mach + 0.1504*mach^2)*CL^2
        elif 0.78 < mach <= 0.85:
            return (-99.3434 + 380.888*mach - 486.8*mach^2 + 207.408*mach^3)*CL^2
        else:
            raise Exception("Mach number is beyond defined limits [0.00, 0.85]")



    def d_CD_turn(self, CL, flap_angle=0, phi=None, Nz = None):
        """
        Drag due to banking.

        :param flap_angle: Flap deployment angle [0°, 20° or 45°]
        :param phi: Banking angle [°]
        :param Nz: Load Factor
        :return:
        """

        if phi is None and Nz is not None:
            pass
        elif phi is not None and Nz is None:
            Nz = 1/np.cos(np.radians(phi))
        else:
            raise Exception('Ambiguous or erroneous function arguments.')

        K = self._induced_drag_efficiency_factor(flap_angle)

        return K*CL**2*(Nz-1)

    def aero_coefficient_taxi(self, flap_angle=0, spoilers = False):
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

        return cl , cd














