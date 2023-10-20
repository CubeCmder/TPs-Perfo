import numpy as np
import matplotlib.pyplot as plt

from aircraft import Aircraft
from atmos import pressure_from_alt, temp_from_alt
from climb_descent import get_gradient, get_AF
from velocities import get_dynamic_pressure

#### QUESTION 5
"""
Considérez une montée sous les conditions suivantes :

 Flap 20 / LG up
 Hp = 5,000 ft
 CG=25% MAC#
 Régime moteur MTO OEI
 1.15 Vsr (vol à CAS constant)
 Vent nul
 Température = 38 degrés C
 Nz=1

Calculez le poids maximum de l’avion auquel il est possible de monter avec un gradient 3.0%.
"""

aircraft = Aircraft()
flap_angle = 20
lg_up = True
hp = 5000
cg=0.25
RM = 'MTO'
n_engines = 1
RV = 'CAS'
T_k = 38+273.15
Nz = 1
vsr_mult = 1.15

gradient = 3/100
# grad = (T/W-CD/CL)/(1+AF)
p = pressure_from_alt(hp)
T_std = temp_from_alt(hp)
dISA = T_k - T_std

W_range = np.linspace(30000, 60000, 10000)
gradients = []
for idx, weight in enumerate(W_range):
    mach = aircraft.get_mach_from_VSR(p, weight, vsr_mult, Nz, flap_angle, lg_up)
    q = get_dynamic_pressure(p, T_k, mach=mach)

    AF = get_AF(RV, hp, mach, T_k, T_std)
    Thrust = aircraft.get_thrust(RM, hp, mach, T_k, n_engines)
    CL = aircraft.get_lift_coefficient(Nz=Nz, weight=weight, q=q, cg=cg, S_ref=aircraft.S)
    CD = aircraft.get_drag_coefficient(CL, flap_angle, mach, LDG=not lg_up, NZ=1, OEI=1, q=q, thrust=Thrust)['CDtot']
    gradient = get_gradient(AF, Thrust, weight, CD, CL)

    gradients.append(gradient*100)

gradients = np.array(gradients)
idx = np.argmin(np.abs((gradients - 3)))
gamma_i = gradients[idx]
weight_i = W_range[idx]

plt.plot(W_range, gradients, color='red', label='Gradient = fct(Weight)')
plt.xlabel('Aircraft Weight (lbs)')
plt.grid()
if idx != -1:
    plt.axvline(x=weight_i, color='green', linestyle='--', label='CAS Objective')
    plt.axhline(y=gamma_i, color='green', linestyle='--', label='Altitude')
plt.ylabel('Gradient (%)')
#plt.legend()
plt.title('Gradient as a function of weight')
plt.show()
print('Question 5')
print(f'Weight : {weight_i:0.2f} lb')
print(f'Gradient : {gamma_i:0.2f} deg')

#### QUESTION 6
"""
Considérez une montée à vitesse calibrée constante sous les conditions suivantes :

 Flap 0 / LG up
 vol à CAS constant
 Hp = 12,000 ft
 W = 35,000 lb
 CG=25% MAC
 ISA+10
 Vent nul
 Régime moteur MCT OEI
 Nz=1

Quelle est la vitesse de montée (Mach) qui maximise le gradient de montée, 
et quel est le gradient de montée correspondant?
"""

aircraft = Aircraft()
flap_angle = 0
lg_up = True
hp = 12000
cg=0.25
RM = 'MCT'
n_engines = 1
RV = 'CAS'
dISA = 10
Nz = 1
weight = 35000

p = pressure_from_alt(hp)
T_std = temp_from_alt(hp)
T_k = T_std + dISA

mach_range = np.linspace(0.1, 0.7, 1000)
gradients = []
for idx, mach in enumerate(mach_range):
    q = get_dynamic_pressure(p, T_k, mach=mach)

    AF = get_AF(RV, hp, mach, T_k, T_std)
    Thrust = aircraft.get_thrust(RM, hp, mach, T_k, n_engines)
    CL = aircraft.get_lift_coefficient(Nz=Nz, weight=weight, q=q, cg=cg, S_ref=aircraft.S)
    CD = aircraft.get_drag_coefficient(CL, flap_angle, mach, LDG=not lg_up, NZ=1, OEI=1, q=q, thrust=Thrust)['CDtot']
    gradient = get_gradient(AF, Thrust, weight, CD, CL)

    gradients.append(gradient*100)

gradients = np.array(gradients)
idx = np.argmax(gradients)
gamma_i = gradients[idx]
mach_i = mach_range[idx]

plt.plot(mach_range, gradients, color='red')
plt.xlabel('Mach ( )')
plt.grid()
if idx != -1:
    plt.axvline(x=mach_i, color='green', linestyle='--', label='CAS Objective')
    plt.axhline(y=gamma_i, color='green', linestyle='--', label='Altitude')
plt.ylabel('Gradient (%)')
plt.title('Gradient as a function of Mach number')
plt.show()
print('Question 6')
print(f'Mach Number : {mach_i:0.2f} ')
print(f'Gradient : {gamma_i:0.2f} deg')

#### QUESTION 7
"""
Considérez un vol en croisière sous les conditions suivantes :

 Flap 0 / LG up
 W = 47,000 lb
 Hp = 30,000 ft
 CG=25% MAC
 ISA+25
 Régime moteur MCR AEO
 Nz=1

À quelle vitesse maximale (Mach demandé) l’avion peut-il voler en palier? (Ne pas considérer les requis opérationnels)
"""

aircraft = Aircraft()
flap_angle = 0
lg_up = True
weight = 47000
hp = 30000
cg=0.25
dISA = 25
RM = 'MCR'
n_engines = 2
Nz = 1

p = pressure_from_alt(hp)
T_std = temp_from_alt(hp)
T_k = T_std + dISA

mach_range = np.linspace(0.1, 0.83, 40000)
thrusts = []
drags = []
for idx, mach in enumerate(mach_range):
    q = get_dynamic_pressure(p, T_k, mach=mach)

    Thrust = aircraft.get_thrust(RM, hp, mach, T_k, n_engines)
    CL = aircraft.get_lift_coefficient(Nz=Nz, weight=weight, q=q, cg=cg, S_ref=aircraft.S)
    CD = aircraft.get_drag_coefficient(CL, flap_angle, mach, LDG=not lg_up, NZ=1, OEI=0, q=q, thrust=Thrust)['CDtot']
    Drag = q*CD*aircraft.S

    drags.append(Drag)
    thrusts.append(Thrust)

drags = np.array(drags)
thrusts = np.array(thrusts)

#idx = max(np.where(np.abs(thrusts-drags)<0.1))[0]
idx = max(np.abs(thrusts-drags).argsort()[:2])
drag_i = drags[idx]
mach_i = mach_range[idx]

fig, ax = plt.subplots()
ax.plot(mach_range, drags, color='red', label='Drag')
ax.plot(mach_range, thrusts, color='orange', label='Thrust')
ax.set(xlabel='Mach', ylabel='Force (lbs)')
ax.grid()
if idx != -1:
    ax.axvline(x=mach_i, color='green', linestyle='--')
    ax.axhline(y=drag_i, color='green', linestyle='--')
plt.title('Longitudinal forces as a function of Mach number')
plt.legend()
plt.show()
print('Question 7')
print(f'Mach Number : {mach_i:0.2f} ')