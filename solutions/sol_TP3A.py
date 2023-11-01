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
        r'C:\Users\Admin\OneDrive - polymtl.ca\Poly & School\POLY\2023 Automne\Performance\TPs\solutions\AER8375_TP3A_sol.xlsx')[
        'Sheet1']
    cases = {}

    for row in sol.iter_rows(min_row=15, max_row=19, min_col=1, max_col=9, values_only=True):
        case = {}

        if row[0] is None:
            continue
        # Initial pressure altitude
        if isinstance(row[1], str):
            hpi = extract_float_from_string(row[1])[0]
        else:
            hpi = float(row[1])

        # Final pressure altitude
        if isinstance(row[2], str):
            hpf = extract_float_from_string(row[2])[0]
        else:
            hpf = float(row[2])

        # Temperature deviation
        dISA = extract_float_from_string(row[3])
        if len(dISA):
            dISA = dISA[0]
        else:
            dISA = 0

        # Wind Speed (kts) - Note: Rear wind positive
        if isinstance(row[4], str):
            Vwind = extract_float_from_string(row[4])[0]
        else:
            Vwind = float(row[4])


        # Aircraft CAS at 10000 < hp < hp_tr
        if isinstance(row[5], str):
            V_CAS = extract_float_from_string(row[5])[0]
        else:
            V_CAS = float(row[5])

        # Aircraft Mach at  hp > hp_tr
        if isinstance(row[6], str):
            mach = extract_float_from_string(row[6])[0]
        else:
            mach = float(row[6])

        # Aircraft Initial Weight at hpi
        if isinstance(row[7], str):
            Wi = extract_float_from_string(row[7])[0]
        else:
            Wi = float(row[7])

        # Aircraft Engine Rating
        RM = row[8]

        case['hpi'] = hpi
        case['hpf'] = hpf
        case['dISA'] = dISA
        case['V_wind'] = Vwind
        case['V_CAS'] = V_CAS
        case['Mach'] = mach
        case['Wi'] = Wi
        case['RM'] = RM

        cases[row[0]] = case

    results = {}
    for id in cases:
        case = cases[id]

        hpi = case['hpi']
        hpf = case['hpf']
        dISA = case['dISA']
        Vwind = case['V_wind']
        V_CAS_cst = case['V_CAS']
        Mach_cst = case['Mach']
        Wi = case['Wi']
        Wf = Wi
        RM = case['RM']
        cg = 0.25

        n_engines = 2
        if 'AEO' in RM:
            OEI_tag = False
        elif 'OEI' in RM:
            OEI_tag = True
            n_engines = 1
        else:
            OEI_tag = False

        RM = RM.split()[0]



        int_hp_resol = 1000 # The altitude résolution used during integration [ft]

        if hpi < hpf: # Climb
            int_steps = np.arange(start=hpi, stop=hpf + int_hp_resol, step=int_hp_resol)

            if int_steps[-1] != hpf:
                int_steps[-1] = hpf

            if hpi < 10000 < hpf and 10000 not in int_steps: # Account for a discontinuity at 10000 ft
                int_steps = np.insert(int_steps, min(np.where(10000 < int_steps))[0], 10000)

            Ps_tr = get_impact_pressure(P_0, V_CAS_cst / a0) / ((1 + 0.2 * (Mach_cst) ** 2) ** 3.5 - 1)
            hp_tr = get_pressure_altitude(Ps_tr)
            hp_tr = round(hp_tr/100)*100

            if hpi < hp_tr < hpf and hp_tr not in int_steps: # Account for a discontinuity at 10000 ft
                int_steps = np.insert(int_steps, min(np.where(hp_tr < int_steps))[0], hp_tr)

            while True:
                Wfuel = 0
                t_total = 0
                dist_total = 0
                W1 = Wi
                W2 = Wi
                for idx, hp1 in enumerate(int_steps[:-2]):
                    hp2 = int_steps[idx+1]
                    hp_moy = hp2/2+hp1/2

                    p = pressure_from_alt(hp_moy)
                    T_ISA = temp_from_alt(hp_moy)
                    T_case = T_ISA + dISA

                    delta_hp = hp2-hp1
                    delta_hg = delta_hp * T_case/T_ISA
                    W_avg = (W2+W1)/2

                    if hp1 < hp_tr:
                        if hp1 < 10000:
                            V_CAS = 250
                        else:
                            V_CAS = V_CAS_cst
                        Mach = get_mach_from_calibrated_airspeed(p, V_CAS)
                        V_TAS = get_true_airspeed(p, mach, temp = T_case, knots=False)
                        CL = aircraft.get_lift_coefficient(Nz=1, weight=W_avg,
                                                           q=get_dynamic_pressure(p, T_case, mach=Mach),
                                                           S_ref=aircraft.S)  # Not corrected for CG
                        AoA = aircraft.get_aoa(CL, 0, cg, True)
                        gradient, ROC, ROCp, AF, a_xfp = climb_descent(aircraft, 'CAS', hp_moy, T_case, W_avg,
                                                                       RM, n_engines, AoA, 0, True, cg,
                                                                       CAS=V_CAS)

                    else:
                        Mach = Mach_cst
                        V_TAS = get_true_airspeed(p, mach, temp=T_case, knots=False)
                        CL = aircraft.get_lift_coefficient(Nz=1, weight=W_avg,
                                                           q=get_dynamic_pressure(p, T_case, mach=Mach),
                                                           S_ref=aircraft.S)  # Not corrected for CG
                        AoA = aircraft.get_aoa(CL, 0, cg, True)
                        gradient, ROC, ROCp, AF, a_xfp = climb_descent(aircraft, 'Mach', hp_moy, T_case, W_avg,
                                                                       RM, n_engines, AoA, 0, True, cg,
                                                                       Mach=Mach)

                    if ROC < 100:
                        break

                    d_time_idx = delta_hg/(ROC/60)
                    d_dist_idx = (V_TAS-Vwind) * d_time_idx
                    thrust = aircraft.get_thrust(RM.split()[0], hp_moy, Mach, T_case, n_engines=n_engines)
                    fuel_burn_rate = aircraft.get_fuel_burn_rate(p, thrust)
                    fuel_burned_idx = fuel_burn_rate*d_time_idx
                    Wfuel += fuel_burned_idx

                    t_total += d_time_idx
                    dist_total += d_dist_idx

                    W1 = W2
                    W2 = W1-fuel_burned_idx




                if abs((Wi - Wf)-Wfuel) <= 20:
                    Wf = Wi-Wfuel
                    break
                else:
                    Wf = Wi-Wfuel

        else: # Descent
            pass


        results[id] = {}
        results[id]['Time'] = t_total
        results[id]['Distance'] = dist_total/6076
        results[id]['Wf'] = Wf
        results[id]['Hp_tr'] = HP_tr
        results[id]['Acc'] = Acc
        results[id]['AccTASAvg'] = AccTASAvg
        results[id]['AccTime'] = AccTime
        results[id]['AccDist'] = AccDist
        results[id]['AccFuel'] = AccFuel
        results[id]['AccWfi'] = AccWfi


    headers = ['case', 'Time', 'Distance', 'Wf', 'Hp_tr', 'Acc', 'AccTASAvg', 'AccTime', 'AccDist','AccFuel','AccWfi']

    print(tabulate([[name, *inner.values()] for name, inner in results.items()],
                   headers=headers,
                   tablefmt="github",
                   floatfmt=(".0f", ".4f", ".1f", ".1f", ".6f", ".5f")))
