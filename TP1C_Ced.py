import matplotlib.pyplot as plt
import numpy as np

from atmos import *
from velocities import *

print('_________ Question 2 _________')
print('Un avion vola à 240 kts CAS a une altitude pression de 3400 ft.')
print("Calculer la pression d'impact. \n")

Vc = 240
h = 3400
p = pressure_from_alt(h)
mach = get_mach_from_calibrated_airspeed(p, Vc)
qc = get_impact_pressure(p, mach)
print(f'Données:')
print(f'Vc : {Vc:0.2f} kts')
print(f'Altitude : {h:0.2f} ft')

print(f'\nRecherché:')
print(f'qc : {qc:0.4f} psf')

print('\n_________ Question 4 _________')
print('Un avion vole à 270 kts CAS à une altitude pression de 3000 ft.')
print(
    "Est-il possible que la vitesse vraie (TAS) soit égale à sa vitesse CAS ? Si oui, sous quelle(s) condition(s)? \n")

CAS = 270  # kts
Hp = 3000  # ft

TAS = []
conditions_dISA = []

dISA_min = -30
dISA_max = 15

for dISA in range(dISA_min, dISA_max):
    P, Rho, T = get_atmos_from_dISA(Hp, dISA)
    mach = get_mach_from_calibrated_airspeed(P, CAS)

    TAS.append(get_true_airspeed(P, mach, temp=T))
    conditions_dISA.append([P, Rho, T])

TAS = np.array(TAS)
idx = np.argmin(np.abs((CAS - TAS)))
TAS_i = TAS[idx]
dISA_i = list(range(dISA_min, dISA_max))[idx]

# print(f'Température (K): {conditions_dISA[idx][2]:0.2f}')
# print(f'Pression (psf): {conditions_dISA[idx][0]:0.2f}')

if abs((CAS - TAS_i) / TAS_i) < 0.01:
    print(f'Si dISA est égal à environ {dISA_i:0.0f}°C, la valeur de TAS sera égal à la valeur de CAS.')
    print(f'TAS: {TAS_i:0.2f}')
else:
    print("Il n'y a pas de valeur de dISA où la valeur de TAS est égal à la valeur de CAS.")

plt.plot(list(range(dISA_min, dISA_max)), TAS, color='orange')
#plt.title("TAS en fonction de delta ISA")
plt.xlabel('Delta ISA (°C)')
plt.grid()
if idx != -1:
    plt.axhline(y=CAS, color='green', linestyle='--', label='CAS')
    plt.axvline(x=dISA_i, color='green', linestyle='--', label='dISA')
plt.ylabel('TAS (kts)')
plt.show()

print('\n_________ Question 5 _________')
print('Un avion vole à 275 kts CAS.')
print("Est-il possible que le nombre de Mach soit égal à 0.76 ? Si oui, sous quelle(s) condition(s)?\n")

CAS = 275  # kts

idx = -1
altitudes = np.linspace(20000, 45000, 32500)
conditions_dISA = []
Mach = []

for h in altitudes:
    P, Rho, T = get_atmos_from_dISA(h, 0)
    mach = get_mach_from_calibrated_airspeed(P, CAS)

    Mach.append(mach)
    conditions_dISA.append([P, Rho, T])

Mach = np.array(Mach)
idx = np.argmin(np.abs((Mach - 0.76)))
Mach_i = Mach[idx]
alt_i = altitudes[idx]

if abs(Mach_i - 0.76) / Mach_i < 0.01:
    print(f'À une altitude de {alt_i:0.0f} ft, le nombre de mach sera de 0.76.')
else:
    print("Il n'y a pas de conditions où la valeur de le nombre de mach est de 0.76.")

plt.plot(altitudes, Mach, color='blue')
#plt.title("Mach en fonction de l'altitude")
plt.xlabel('Altitude (ft)')
plt.grid()
if idx != -1:
    plt.axhline(y=0.76, color='green', linestyle='--', label='Mach Objective')
    plt.axvline(x=alt_i, color='green', linestyle='--', label='Altitude')
plt.ylabel('Mach')
plt.show()


print('\n_________ Question 6 _________')
print("Un avion vole à une altitude pression de 15 000 ft à une température de 0 °C en conditions givrantes et de \n"
      "la glace s'accumule sur les surfaces non protégées de l'avion. Le pilote désire voler l'avion de façon à ce \n"
      "que la température totale soit d'au moins 10 degrés C aux points de stagnation sur l'avion, ce qui permettra \n"
      "de faire fondre l'accumulation de glace.")
print("À quelle vitesse calibrée minimum l'avion devra-t-il voler ? \n")

Hp = 15000  # ft
T = 273.15  # K

Mach = np.linspace(0.3, 0.5, 1000)
temps = []
CAS = []
idx = -1
for i in Mach:
    dISA = get_delta_ISA(Hp, T)
    P, Rho, T = get_atmos_from_dISA(Hp, dISA)
    temps.append(get_total_temperature(T, i))
    CAS.append(get_calibrated_airspeed(P, i))

temps = np.array(temps) - 273.15
idx = np.argmin(np.abs((temps - 10)))
CAS_i = CAS[idx]
T0_i = temps[idx]

if abs(T0_i - 10) / T0_i < 0.01:
    print(f'À une vitesse calibrée minimale de {CAS_i:0.1f} ft, la température totale aux points de stagnation sera de 10°C.')

else:
    print("Il n'y a pas de conditions où la température totale aux points de stagnation sera de 10°C.")

plt.plot(CAS, temps, color='red')
#plt.title("La température totale (°C) en fonction de la vitesse calibrée (kts).")
plt.xlabel('CAS (kts)')
plt.grid()

if idx != -1:
    plt.axvline(x=CAS_i, color='green', linestyle='--', label='CAS Objective')
    plt.axhline(y=T0_i, color='green', linestyle='--', label='Altitude')
plt.ylabel('Température Totale (°C)')
plt.show()
