import matplotlib.pyplot as plt
import numpy as np

from atmos import *
from velocities import *

print('_________ Question 2 _________')
print('Un avion vola à 240 kts CAS a une altitude pression de 3400 ft.')
print("Calculer la pression d'impact. \n\n")

Vc = 240
h = 3400
p = pressure_from_alt(h)
mach = get_mach_from_calibrated_airspeed(p, Vc)
qc = get_impact_pressure(p, mach)
print(f'Données:')
print(f'Vc : {Vc:0.2f} kts')
print(f'Altitude : {h:0.2f} ft')

print(f'Recherché:')
print(f'qc : {qc:0.4f} psf')

print('_________ Question 4 _________')
print('Un avion vola à 270 kts CAS a une altitude pression de 3000 ft.')
print("Est-il possible que TAS == CAS? Si oui, sous quelle(s) conditions(s)? \n")

Vc = 270
h = 3000

dISA = range(-55,55)
V = np.zeros(len(dISA))

test = False

for i in range(len(dISA)):
    p,rho,t = get_atmos_from_dISA(hp=h,dISA=dISA[i])
    mach = get_mach_from_calibrated_airspeed(p,Vc)
    V[i] = get_true_airspeed(p,mach,temp=t)

    if V[i] == Vc:
        print(f'Données:')
        print(f'Vc : {Vc:0.2f} kts')
        print(f'Altitude : {h:0.2f} ft')

        print(f'Recherché:')
        print(f'Température : {t:0.4f} k')
        print(f'Delta ISA : {dISA[i]:0.4f} k')
        test = True


if test:
    print(f'Les conditions ci-dessus montre les conditions pour Vc=V')
else:
    print(f'FAUX\nEn variant le delta ISA (la seul variable non fixé), acunne conditions à été trouvé pour que CAS == TAS')

plt.plot(dISA,V)
plt.title("True air speed en fonction de delta ISA ")
plt.xlabel('Delta ISA')
plt.ylabel('True air speed (knts)')
plt.show()

print('_________ Question 5 _________')
print('Un avion vola à 275 kts CAS.')
print("Est-il possible que mach == 0.76? Si oui, sous quelle(s) conditions(s)? \n")

Vc = 275
mach = 0.76

p = get_pressure_from_mach_and_CAS(Vc,mach)

print(f'Données:')
print(f'Vc : {Vc:0.2f} kts')
print(f'Altitude : {h:0.2f} ft')

print(f'Recherché:')
print(f'Température : {t:0.4f} k')
print(f'Delta ISA : {dISA[i]:0.4f} k')

