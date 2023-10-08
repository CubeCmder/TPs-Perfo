import numpy as np
import openpyxl
import re

from tabulate import tabulate

from velocities import *
from atmos import *
from aircraft import Aircraft

def extract_float_from_string(s):
    # Define a regular expression pattern to match floating-point numbers
    # This pattern will match numbers with optional positive/negative signs, integer part, optional decimal point, and fractional part.
    pattern = r"[-+]?\d*\.\d+|\d+\.\d*|\d+"

    # Use the findall() function to extract all matching floats from the input string
    floats = re.findall(pattern, s)

    # Convert the extracted strings to actual float values and store them in a list
    float_values = [float(num) for num in floats]
    return float_values
    l = []
    for t in s.split():
        try:
            l.append(float(t))
        except ValueError:
            pass
    return l

if __name__ == '__main__':

    aircraft = Aircraft()

    sol = openpyxl.load_workbook(r'C:\Users\Admin\OneDrive - polymtl.ca\Poly & School\POLY\2023 Automne\Performance\TPs\solutions\AER8375_TP2A_sol.xlsx')['Sheet1']
    cases = {}

    for row in sol.iter_rows(min_row=15, max_row=22, min_col=1, max_col=9, values_only=True):
        case = {}

        if isinstance(row[1], str):
            hp = extract_float_from_string(row[1])[0] # pressure altitude
        else:
            hp = float(row[1])

        p = pressure_from_alt(hp)

        dISA = extract_float_from_string(row[2])
        if len(dISA):
            dISA = dISA[0]
        else:
            dISA = 0


        temp = temp_from_alt(hp) + dISA

        weight = float(row[4])

        cg = float(row[5])

        flap_angle = float(row[6].split('/')[0].strip())

        g_u = row[6].split('/')[1].strip()

        if g_u == 'up':
            g_u = True
        else:
            g_u = False


        if isinstance(row[7], str):
            Nz = extract_float_from_string(row[7])[0]
            Nz = 1 / np.cos(np.radians(Nz))
        else:
            Nz = float(row[7])

        if 'Mach' in row[3]:
            mach = extract_float_from_string(row[3])[0]
        elif 'Vsr' in row[3]:
            vsr_mult = extract_float_from_string(row[3])[0]

            mach = aircraft.get_mach_from_VSR(p, weight, vsr_mult, Nz, flap_angle, g_u)
            #Vs = get_stall_speed(weight, Nz,hp,aircraft.S,CL)
            #mach = get_mach_from_calibrated_airspeed(p, Vs * 0.592484)

            #CLMAX = aircraft.get_CL_max(flap_angle,gear_up=g_u)*(1 + aircraft.mac / aircraft.lt * (aircraft.FWDCG - cg/100))
            #Vs = get_stall_speed(weight, Nz, hp, aircraft.S, CLMAX)
            #mach = get_mach_from_calibrated_airspeed(p, vsr_mult * Vs * 0.592484)
        elif 'CAS' in row[3]:
            CAS = extract_float_from_string(row[3])[0]
            mach = get_mach_from_calibrated_airspeed(p, CAS)

        RM = row[8]

        case['hp'] = hp
        case['T'] = temp
        case['mach'] = mach
        case['weight'] = weight
        case['cg'] = cg/100
        case['flap'] = flap_angle
        case['gear_up'] = g_u
        case['Nz'] = Nz
        case['RM'] = RM

        cases[row[0]] = case

    results = {}
    for id in cases:
        case = cases[id]

        hp = case['hp']
        T = case['T']
        mach = case['mach']
        weight = case['weight']
        cg = case['cg']
        flap_angle = case['flap']
        g_u = case['gear_up']
        Nz = case['Nz']
        RM = case['RM']

        n_engines = 2
        if 'AEO' in RM:
            OEI_tag = False
        elif 'OEI' in RM:
            OEI_tag = True
            n_engines = 1
        else:
            OEI_tag = False

        p = pressure_from_alt(hp)

        q = get_dynamic_pressure(p, mach=mach)

        # OUTPUTS

        if len(extract_float_from_string(RM)):
            thrust = extract_float_from_string(RM)[0]
        else:
            if 'idle' in RM:
                thrust = aircraft.get_thrust('ID', hp, mach, T, n_engines=n_engines)
            else:
                thrust = aircraft.get_thrust(RM.split()[0], hp, mach, T, n_engines=n_engines)

        CL = aircraft.get_lift_coefficient(Nz=Nz, weight=weight, q=q, cg=cg, S_ref=aircraft.S) # Not corrected for CG

        CD = aircraft.get_drag_coefficient(CL, flap_angle, mach, Nz=Nz, OEI=OEI_tag, LDG=not g_u, q=q, thrust=thrust)
        CDp = CD['CDp']
        CDi = CD['CDi']
        CDcomp = CD['CDcomp']
        CDctl = CD['CDctl']
        CDwm = CD['CDwm']
        CDtot = CD['CDtot']

        lift_to_drag = CL / CDtot

        drag_force = get_dynamic_pressure(p,T,mach=mach) * aircraft.S * CDtot
        lift_force = get_dynamic_pressure(p, T, mach=mach) * aircraft.S * CL

        # load factor at stall warning
        Nz_sw = aircraft.get_Nz_stall_warning(p, mach, flap_angle, weight, cg, g_u)
        # bank angle at stall warning
        phi_sw = aircraft.get_Nz_stall_warning(p, mach, flap_angle, weight, cg, g_u, True)
        Nz_buffet = aircraft.NZ_buffet(mach, CL, cg)

        AoA = aircraft.get_aoa(CL, flap_angle, cg, g_u)
        results[id] = {}
        results[id]['CD'] = CDtot
        results[id]['CL'] = CL
        results[id]['L/D'] = CL/CDtot
        results[id]['CDp'] = CDp
        results[id]['CDi'] = CDi
        results[id]['CDcomp'] = CDcomp
        results[id]['CDctl'] = CDctl
        results[id]['CDwm'] = CDwm
        results[id]['Thrust'] = thrust
        results[id]['Drag'] = drag_force
        results[id]['AOA'] = AoA
        results[id]['Nz_sw'] = Nz_sw
        results[id]['Phi_sw'] = phi_sw
        results[id]['Nz_buffet'] = Nz_buffet


    headers = ['Cas', 'Cd', 'Cl', 'L/D', 'Cdp', 'Cdi', 'Cdcomp', 'Cdcntl', 'Cdwm', 'Pousée Totale', 'Trainée', 'AOA', 'Nz_sw', 'Phi_sw', 'Nz_buffet']

    print(tabulate([[name, *inner.values()] for name, inner in results.items()],
                   headers = headers,
                   tablefmt="github",
                   floatfmt=(".0f",".6f", ".5f", ".5f", ".5f", ".5f", ".5f", ".5f", ".5f", ".1f", ".4f")))



