import matplotlib.pyplot as plt
import numpy as np

from atmos import *
from velocities import *


print('_________ Question 4 _________')
print('Un avion vola à 270 kts CAS a une altitude pression de 3000 ft.')
print("Est-il possible que TAS == CAS? Si oui, sous quelle(s) conditions(s)? \n")

Vc = 270
h = 3000

dISA = range(-30,15)
V = np.zeros(len(dISA))


test = False
idx=-1
for i in range(len(dISA)):
    p, rho, t = get_atmos_from_dISA(hp=h, dISA=dISA[i])
    mach = get_mach_from_calibrated_airspeed(p, Vc)
    V[i] = get_true_airspeed(p, mach, temp=t)

    if abs((Vc-V[i])/V[i]) < 0.001:
        print(f'Données:')
        print(f'Vc : {Vc:0.2f} kts')
        print(f'Altitude : {h:0.0f} ft')

        print(f'Recherché:')
        print(f'Température : {t:0.2f} k')
        print(f'Delta ISA : {dISA[i]:0.2f} k')
        dISA_V_Vc = dISA[i]
        idx=i
        test = True

TAS = np.array(V)
idx = np.argmin(np.abs((Vc-TAS)))


if test:
    print(f'Les conditions ci-dessus montre les conditions pour Vc=V, soit delta ISA de {dISA[i]:0.2f} k')
else:
    print(f'FAUX\nEn variant le delta ISA (la seul variable non fixé), acunne conditions à été trouvé pour que CAS == TAS')

plt.plot(list(dISA),V, color='blue',zorder=0)
plt.plot([dISA[0],dISA[-1]],[Vc,Vc], color='orange',zorder=5)
plt.title("True air speed en fonction de delta ISA ")
plt.xlabel('Delta ISA')
plt.scatter(dISA_V_Vc,Vc, color='orange',zorder=10)
plt.text(dISA_V_Vc+1, Vc-1, f'({dISA_V_Vc:0.0f}, {Vc:0.0f})',zorder=15)
plt.ylabel('True air speed (knts)')
plt.show()


print('_________ Question 5 _________')
print('Un avion vole à 275 kts CAS.')
print("Est-il possible que mach == 0.76? Si oui, sous quelle(s) conditions(s)? \n")

CAS = 275  # kts

i = -1
alt = np.linspace(20000, 45000, 32500)
conditions_dISA = []
Mach = []

for h in alt:
    P, Rho, T = get_atmos_from_dISA(h, 0)
    mach = get_mach_from_calibrated_airspeed(P, CAS)

    Mach.append(mach)
    conditions_dISA.append([P, Rho, T])

Mach = np.array(Mach)
i = np.argmin(np.abs((Mach - 0.76)))
Mach_i = Mach[i]
alt_i = alt[i]
P, Rho, T = get_atmos_from_dISA(alt_i, 0)
mach = get_mach_from_calibrated_airspeed(P, CAS)

if abs(Mach_i - 0.76) / Mach_i < 0.01:
    print(f'À une altitude de {alt_i:0.0f} ft, le nombre de mach sera de 0.76.')
else:
    print("Il n'y a pas de conditions où la valeur de le nombre de mach est de 0.76.")


print(f'Données:')
print(f'Vc : {Vc:0.2f} kts')
print(f'mach : 0.76')

print(f'Recherché:')
print(f'Altitude : {alt_i:0.0f} ft')
print(f"Vrai\n La condition est obtenues à une altitude de {alt_i:0.0f}.")

dISA = [-30,0,15]
h = range(0, 40000,1000)
mach = np.zeros([len(h),len(dISA)])

plt.plot(alt,Mach,color='blue',zorder=0)
plt.plot([0,40000],[0.76,0.76], color='orange',zorder=5)
plt.title("Mach number in function of altitude (delta ISA of 0) ")
plt.ylabel('mach')
plt.xlabel('Altitude (ft)')
plt.xlim(20000, 40000)
plt.ylim(0.6, 1)
plt.scatter(alt_i,Mach_i, color='orange',zorder=10)
plt.text(alt_i, Mach_i-0.02, f'({alt_i:0.0f}, {Mach_i:0.2f})',zorder=15)
plt.legend(['mach fct(altitude)','Target'])
plt.show()



print('_________ Question 6 _________')
print('Un avion vole à une altitude de 15 000ft a une temperature de 0degC.')
print('Le pilote veut voler a 10degC au point de stagnation.')
print("A quelle vitesse CAS minimu l'avion devra-t-il voler? \n")

h = 15000
t_0 = 0+273.15
t_total = 10+273.15

mach = ((t_total/t_0 - 1)/0.2)**0.5
dISA = get_delta_ISA(hp=h,T=t_0)
p,rho,t = get_atmos_from_dISA(hp=h,dISA=dISA)
Vc = get_calibrated_airspeed(p,mach)

print(f'Données:')
print(f'Altitude : {h:0.0f} ft')
print(f'Temperature total : {t_total:0.2f}')

print(f'Recherché:')
print(f'CAS : {Vc:0.2f} knts')
print(f"Pour obtenir la temperature totale de 10degC, le pilote devra manoeuvrer l'avion à {Vc:0.2f} knts CAS.")

