import numpy as np
import openpyxl
import re

from tabulate import tabulate

from velocities import *
from atmos import *
from aircraft import Aircraft
from climb_descent import climb_descent

aircraft = Aircraft()

#CAS 1
hp = 15000
dISA = 0
Vwind = 0
W = 40000
V = []
V_type = "Vmd"

KCAS, V_g, Mach, SAR, SR, Wf = aircraft.cruise(hp=hp, dISA=dISA, Vwind=Vwind, V=V, V_type=V_type, W=W, OEI_tag=False)
Cas1 = [KCAS,V_g,Mach,SAR,SR,Wf]
print(Cas1)

#CAS 2
hp = 31000
dISA = 10
Vwind = 20
W = 45000
V = 0.74
V_type = "Mach"

KCAS, V_g, Mach, SAR, SR, Wf = aircraft.cruise(hp, dISA, Vwind, V, V_type, W, OEI_tag=False)
Cas2 = [KCAS,V_g,Mach,SAR,SR,Wf]
print(Cas2)

#CAS 3
hp = 35000
dISA = 10
Vwind = -50
W = 40000
V = []
V_type = "LRC"

KCAS, V_g, Mach, SAR, SR, Wf = aircraft.cruise(hp, dISA, Vwind, V, V_type, W, OEI_tag=False)
Cas3 = [KCAS,V_g,Mach,SAR,SR,Wf]
print(Cas3)
#
# def extract_float_from_string(s):
#     # Define a regular expression pattern to match floating-point numbers
#     # This pattern will match numbers with optional positive/negative signs, integer part, optional decimal point, and fractional part.
#     pattern = r"[-+]?\d*\.\d+|\d+\.\d*|\d+"
#
#     # Use the findall() function to extract all matching floats from the input string
#     floats = re.findall(pattern, s)
#
#     # Convert the extracted strings to actual float values and store them in a list
#     float_values = [float(num) for num in floats]
#     return float_values
#     l = []
#     for t in s.split():
#         try:
#             l.append(float(t))
#         except ValueError:
#             pass
#     return l
#
#
# if __name__ == '__main__':
#
#     aircraft = Aircraft()
#
#     sol = openpyxl.load_workbook(
#         r'C:\Users\hugoa\OneDrive - polymtl.ca\Ecole\A2023\AER8375\TPs-Perfo\solutions\AER8375_TP3B_sol.xlsx')[
#         'Feuil1']
#     cases = {}
#
#     for row in sol.iter_rows(min_row=14, max_row=16, min_col=1, max_col=6, values_only=True):
#         case = {}
#
#         if row[0] is None:
#             continue
#         # Pressure altitude
#         if isinstance(row[1], str):
#             hp = extract_float_from_string(row[1])[0]
#         else:
#             hp = float(row[1])
#
#
#         # Temperature deviation
#         dISA = extract_float_from_string(row[2])
#         if len(dISA):
#             dISA = dISA[0]
#         else:
#             dISA = 0
#
#         # Wind Speed (kts) - Note: Rear wind positive
#         if isinstance(row[3], str):
#             Vwind = extract_float_from_string(row[3])[0]
#         else:
#             Vwind = float(row[3])
#
#         # Aircraft Initial Weight at hpi
#         if isinstance(row[4], str):
#             Wi = extract_float_from_string(row[4])[0]
#         else:
#             Wi = float(row[4])
#         # Aircraft speed
#         if isinstance(row[6], str):
#             mach = extract_float_from_string(row[6])[0]
#         else:
#             mach = float(row[6])
#
#
#         case['hp'] = hp
#         case['dISA'] = dISA
#         case['V_wind'] = Vwind
#         case['W'] = W
#         case['V'] = V
#
#         cases[row[0]] = case
#
#     results = {}
#     for id in cases:
#         case = cases[id]
#
#         hp = case['hp']
#         dISA = case['dISA']
#         Vwind = case['V_wind']
#         W = case['W']
#         V = case['V']
#
#         n_engines = 2
#         OEI_tag = False
#
#
#         V, V_g, SAR, SR, Wf = aircraft.cruise(hp, dISA, Vwind, V, V_type, W, OEI_tag=False)
#
#         results[id] = {}
#         results[id]['Speed'] = V
#         results[id]['Ground Speed'] = V_g
#         results[id]['Mach'] = Mach
#         results[id]['SAR'] = SAR
#         results[id]['SR'] = SR
#         results[id]['Wf'] = Wf
#
#
#     headers = ['case', 'Speed', 'Ground Speed', 'Mach', 'SAR', 'SR', 'Wf']
#
#     print(tabulate([[name, *inner.values()] for name, inner in results.items()],
#                    headers=headers,
#                    tablefmt="github",
#                    floatfmt=(".0f", ".4f", ".1f", ".1f", ".6f", ".5f")))
#
#     #
#     #
#     # while True:
#     #     Wfuel = 0
#     #     t_total = 0
#     #     dist_total = 0
#     #     W1 = Wi
#     #     W2 = Wi
#     #
#     #     for idx, hp1 in enumerate(int_steps[:-2]):
#     #         hp2 = int_steps[idx + 1]
#     #         hp_moy = hp2 / 2 + hp1 / 2
#     #
#     #         p = pressure_from_alt(hp_moy)
#     #         T_ISA = temp_from_alt(hp_moy)
#     #         T_case = T_ISA + dISA
#     #
#     #         delta_hp = hp2 - hp1
#     #         delta_hg = delta_hp * (T_case / T_ISA)
#     #         W_avg = (W2 + W1) / 2
#     #
#     #         if hp1 == 10000:  # Acceleration Segment
#     #             V1 = 250
#     #             V2 = V_CAS_cst
#     #             V_CAS_AVG = (V1 + V2) / 2
#     #
#     #             Mach1 = get_mach_from_calibrated_airspeed(p, V1)
#     #             Mach2 = get_mach_from_calibrated_airspeed(p, V2)
#     #             Mach = get_mach_from_calibrated_airspeed(p, V_CAS_AVG)
#     #             q = get_dynamic_pressure(p, T_case, mach=Mach)
#     #
#     #             Thrust = self.get_thrust(RM, hp_moy, Mach, T_case, n_engines=n_engines)
#     #             CL = self.get_lift_coefficient(Nz=1, weight=W_avg,
#     #                                            q=get_dynamic_pressure(p, T_case, mach=Mach),
#     #                                            S_ref=self.S)
#     #             CD = self.get_drag_coefficient(CL, 0, Mach, LDG=0, NZ=1, OEI=OEI_tag, q=q,
#     #                                            thrust=Thrust)['CDtot']
#     #             D = q * CD * self.S
#     #
#     #             Acc = (Thrust - D) / W_avg * g
#     #             AccTASAvg = get_true_airspeed(p, Mach, temp=T_case, knots=False)
#     #             AccTime = (get_true_airspeed(p, Mach2, temp=T_case, knots=False) - get_true_airspeed(p, Mach1,
#     #                                                                                                  temp=T_case,
#     #                                                                                                  knots=False)) / Acc
#     #             AccDist = AccTime * (AccTASAvg - knots2fps(Vwind)) / 6076
#     #             fuel_burn_rate = self.get_fuel_burn_rate(hp_moy, Thrust)
#     #             AccFuel = fuel_burn_rate * AccTime
#     #             AccWfi = Wi
#     #
#     #             fuel_burned_idx = AccFuel
#     #             d_time_idx = AccTime
#     #             d_dist_idx = AccDist
#     #
#     #         elif hp1 < hp_tr:
#     #             if hp1 < 10000:
#     #                 V_CAS = 250
#     #             else:
#     #                 V_CAS = V_CAS_cst
#     #             Mach = get_mach_from_calibrated_airspeed(p, V_CAS)
#     #             CL = self.get_lift_coefficient(Nz=1, weight=W_avg,
#     #                                            q=get_dynamic_pressure(p, T_case, mach=Mach),
#     #                                            S_ref=self.S)  # Not corrected for CG
#     #             AoA = self.get_aoa(CL, 0, cg, True)
#     #             gradient, ROC, ROCp, AF, a_xfp = climb_descent(self, 'CAS', hp_moy, T_case, W_avg,
#     #                                                            RM, n_engines, AoA, 0, True, cg,
#     #                                                            CAS=V_CAS)
#     #             V_TAS = get_true_airspeed(p, Mach, temp=T_case, knots=False)
#     #
#     #             d_time_idx = delta_hg / (ROC / 60)
#     #             d_dist_idx = (V_TAS - knots2fps(Vwind)) * d_time_idx
#     #             thrust = self.get_thrust(RM.split()[0], hp_moy, Mach, T_case, n_engines=n_engines)
#     #             fuel_burn_rate = self.get_fuel_burn_rate(hp_moy, thrust)
#     #             fuel_burned_idx = fuel_burn_rate * d_time_idx
#     #
#     #         else:
#     #             Mach = Mach_cst
#     #             CL = self.get_lift_coefficient(Nz=1, weight=W_avg,
#     #                                            q=get_dynamic_pressure(p, T_case, mach=Mach),
#     #                                            S_ref=self.S)  # Not corrected for CG
#     #             AoA = self.get_aoa(CL, 0, cg, True)
#     #             gradient, ROC, ROCp, AF, a_xfp = climb_descent(self, 'Mach', hp_moy, T_case, W_avg,
#     #                                                            RM, n_engines, AoA, 0, True, cg,
#     #                                                            Mach=Mach)
#     #
#     #             V_TAS = get_true_airspeed(p, Mach, temp=T_case, knots=False)
#     #
#     #             d_time_idx = delta_hg / (ROC / 60)
#     #             d_dist_idx = (V_TAS - knots2fps(Vwind)) * d_time_idx
#     #             thrust = self.get_thrust(RM.split()[0], hp_moy, Mach, T_case, n_engines=n_engines)
#     #             fuel_burn_rate = self.get_fuel_burn_rate(hp_moy, thrust)
#     #             fuel_burned_idx = fuel_burn_rate * d_time_idx
#     #
#     #         if ROC < 100:
#     #             break
#     #
#     #         Wfuel += fuel_burned_idx
#     #         t_total += d_time_idx
#     #         dist_total += d_dist_idx
#     #
#     #         W1 = W2
#     #         W2 = W1 - fuel_burned_idx
#     #
#     #     if abs((Wi - Wf) - Wfuel) <= 20:
#     #         Wf = Wi - Wfuel
#     #         break
#     #     else:
#     #         Wf = Wi - Wfuel
