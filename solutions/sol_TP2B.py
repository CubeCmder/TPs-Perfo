import numpy as np
import openpyxl
import re

from tabulate import tabulate

from velocities import *
from atmos import *
from aircraft import Aircraft
from climb_descent import climb_descent


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

    sol = openpyxl.load_workbook(
        r'C:\Users\Admin\OneDrive - polymtl.ca\Poly & School\POLY\2023 Automne\Performance\TPs\solutions\AER8375_TP2B_sol.xlsx')[
        'Sheet1']
    cases = {}

    for row in sol.iter_rows(min_row=15, max_row=17, min_col=1, max_col=10, values_only=True):
        case = {}

        if isinstance(row[1], str):
            hp = extract_float_from_string(row[1])[0]  # pressure altitude
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

            delta = p / P_0
            CL = aircraft.get_CL_max(flap_angle, gear_up=g_u) / vsr_mult ** 2
            mach = np.sqrt((weight * Nz / delta) / 1481.3 / CL / aircraft.S)

        elif 'CAS' in row[3]:
            CAS = extract_float_from_string(row[3])[0]
            mach = get_mach_from_calibrated_airspeed(p, CAS)

        RV = row[8]
        RM = row[9]

        case['hp'] = hp
        case['T'] = temp
        case['mach'] = mach
        case['weight'] = weight
        case['cg'] = cg / 100
        case['flap'] = flap_angle
        case['gear_up'] = g_u
        case['Nz'] = Nz
        case['RM'] = RM
        case['RV'] = RV

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
        RV = case['RV']


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



        CL = aircraft.get_lift_coefficient(Nz=Nz, weight=weight, q=q, cg=cg, S_ref=aircraft.S)  # Not corrected for CG




        AoA = aircraft.get_aoa(CL, flap_angle, cg, g_u)

        if RV == 'Cst M':
            RV = 'Mach'
            CAS_cd = None
            mach_cd = mach
        else:
            RV = 'CAS'
            CAS_cd = get_calibrated_airspeed(p, mach)
            mach_cd = None

        gradient, ROC, ROCp, AF, a_xfp = climb_descent(aircraft, RV, hp, T, weight, RM.split()[0], n_engines,
                                                       AoA, flap_angle, g_u,cg, CAS_cd, mach_cd)

        results[id] = {}
        results[id]['gradient'] = gradient*100
        results[id]['ROC'] = ROC
        results[id]['ROCP'] = ROCp
        results[id]['AF'] = AF
        results[id]['a'] = a_xfp


    headers = ['case', 'gradient', 'ROC', 'ROCp', 'AF', 'a_xfp']

    print(tabulate([[name, *inner.values()] for name, inner in results.items()],
                   headers=headers,
                   tablefmt="github",
                   floatfmt=(".0f", ".4f", ".1f", ".1f", ".6f", ".5f")))
