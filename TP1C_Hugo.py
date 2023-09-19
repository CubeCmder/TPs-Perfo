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
print(f'Le graphique montre que 270 knts CAS serait atteint a un valure irrealiste de delta ISA\n')

print('_________ Question 5 _________')
print('Un avion vola à 275 kts CAS.')
print("Est-il possible que mach == 0.76? Si oui, sous quelle(s) conditions(s)? \n")

Vc = 275
mach = 0.76

p = get_pressure_from_mach_and_CAS(Vc,mach)

print(f'Données:')
print(f'Vc : {Vc:0.2f} kts')
print(f'mach : {mach:0.2f}')

print(f'Recherché:')
print(f'Pression : {p:0.4f} psf')
print(f"Faux\n La valeur E-18 montre que l'operation est en fait impossible.")

dISA = [-30,0,15]
h = range(0, 40000,5000)
mach = np.zeros([len(h),len(dISA)])

for i in range(len(h)):
    for j in range(len(dISA)):
        p,rho,t = get_atmos_from_dISA(h[i],dISA[j])
        mach[i,j] = get_mach_from_calibrated_airspeed(p,Vc)

plt.plot(h,mach[:,0],"bo",h,mach[:,1],"k",h,mach[:,2],"b")
plt.title("Mach number in function of altitude for different delta ISA ")
plt.xlabel('Altitude (ft)')
plt.ylabel('mach')
plt.ylim(0, 8)
plt.xlim(0, 40000)
plt.legend(['dISA = -30','dISA = 0','dISA = 15'])
plt.show()

print(f'Le graphique montre que 0.76 mach ne sera jamais atteint de 275 knts CAS puisque les valeurs de mach')
print(f'pour cette vitesse tourne autour de 5 pour des valeurs atmospherique std.')



print('_________ Question 6 _________')
print('Un avion vole à une altitude de 15 000ft a une temperature de 0degC.')
print('Le pilote veut voler a 10degC au point de stagnation.')
print("A quelle vitesse CAS minimu l'avion devra-t-il voler? \n")

h = 15000
t_0 = 0+273.15
t_stag = 10+273.15

dISA = get_delta_ISA(h,t_stag)
qc,rho,t = get_atmos_from_dISA(hp=h,dISA=dISA)

dISA = get_delta_ISA(h,t_0)
p,rho,t = get_atmos_from_dISA(hp=h,dISA=dISA)

mach = get_mach(p=p,qc=qc)

Vc = get_calibrated_airspeed(p,mach)

print(f'Données:')
print(f'altitude : {h:0.2f} ft')
print(f'Temperature au point de stagnation : {t_stag:0.2f}')

print(f'Recherché:')
print(f'CAS : {Vc:0.4f} knts')

