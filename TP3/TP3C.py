import numpy as np

from velocities import *
from atmos import *
from aircraft import Aircraft
from tabulate import tabulate

aircraft = Aircraft()

# Donnees probleme
OEI_tag = False  # True if OEI
V_CAS_cst = 275  # constant CAS velocity value
Mach_cst = 0.78  # constant Mach velocity value
ROC_Min = 300  # ft/min

W0 = aircraft.OWE + + aircraft.MXFUEL  # Initial weight before taxi and climb
Vwind = 0
dISA = 0

RM_cr = 'MCR'

RM_cl = 'MCL'
Wi_cl = W0 - 200 - 250  # Carburant au début du climb (enlève Taxi et montée jusqu'à 1500ft)
hpi_cl = 1500
hpf_cl = 41000

int_hp_resol = 10  # The altitude résolution used during integration [ft]

hpi = 0
hpf = 0
RM_cl = 'MCL'
V_CAS_cst = 275
Mach_cst = 0.78
hp_cruise = 41000
RM_d = 'ID'
dISA = 0
Vwind = 0
dW_cargo = 20 * 225
Wi_fuel = aircraft.MXFUEL
dW_fuel_to_taxi = 200
dW_fuel_to = 250
dW_fuel_ldg = 200
dW_fuel_ldg_taxi = 100
W_fuel_reserve = 2000 - dW_fuel_ldg_taxi

# t_total, dist_total, dW, __, __, __, __, __, __, __ = aircraft.mission_performance_climb_descent(
#             hpi_cl, hpf_cl, dISA, Vwind, V_CAS_cst, Mach_cst, Wi_cl, RM_cl, int_hp_resol=int_hp_resol, OEI_tag=OEI_tag, ROC_min = ROC_Min)

dist_climb, hp_cruise, dist_cruise, dist_descent, dist_tot, dWf_tot = aircraft.mission(dW_cargo, hpi, hpf, RM_cl,
                                                                                       V_CAS_cst, Mach_cst, hp_cruise,
                                                                                       RM_d, dISA, Vwind,
                                                                                       Wi_fuel=Wi_fuel,
                                                                                       dW_fuel_to_taxi=dW_fuel_to_taxi,
                                                                                       dW_fuel_to=dW_fuel_to, dW_fuel_ldg=dW_fuel_ldg,
                                                                                       dW_fuel_ldg_taxi=dW_fuel_ldg_taxi,
                                                                                       W_fuel_reserve=W_fuel_reserve)

data = [dist_climb, hp_cruise, dist_cruise, dist_descent, dist_tot, dWf_tot]
headers = ['Climb Distance', 'Cruise Altitude', 'Cruise Distance', 'Descent Distance', 'Total Distance',
           'Total Fuel consumption']

print(tabulate([[headers[i], data[i]] for i, v in enumerate(data)],
               tablefmt="github", numalign="right"))
