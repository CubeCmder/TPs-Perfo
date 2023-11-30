import numpy as np
import matplotlib.pyplot as plt

from velocities import *
from atmos import *
from aircraft import Aircraft
from tabulate import tabulate

aircraft = Aircraft()

# Donnees probleme
Hp = 0              #ft
dISA = 20           #C
Flaps = 20          #deg
CG = 20             #%
PISTE = 5700        #ft
Dalignement = 200   #ft
TOD = PISTE-Dalignement #ft

W_max = int(aircraft.MTOW)
W_min = int(W_max - 10000)
W = range(W_min, W_max, 10)
FTOD_AEO = np.zeros(len(W))
TOD_OEI =np.zeros(len(W))
ASD_AEO = np.zeros(len(W))
Longueur_min = np.zeros(len(W))
braking_energy = np.zeros(len(W))
V1VR = 1
Wi = 0

for i in range(len(W)):
    FTOD_AEO[i], TOD_OEI[i], ASD_AEO[i], Longueur_min[i], braking_energy[i] = aircraft.takeoff_run_distances(Hp,W[i],V1VR,0,dISA=dISA,flap_angle=Flaps, cg=CG)

    if abs((TOD - Longueur_min[i]) / TOD) < 0.0005:
        Wi = W[i]
        Longueur_mini = Longueur_min[i]
        print(f'Si le poids est égal à environ {Wi:0.0f}lb, la distance de décollage est la même que la longueur de la piste.')


plt.plot(W,Longueur_min, color='orange')
plt.title("Distance de décollage en fonction du poids")
plt.xlabel('Poids (lb)')
plt.grid()
plt.axhline(y=TOD, color='green', linestyle='--', label='CAS')
if Wi != 0:
    plt.axvline(x=Wi, color='green', linestyle='--', label='dISA')
    plt.ylabel('FTOD (ft)')

else:
    print(f"Il n'y a pas de valeur de poids qui satisfait ce cas")

plt.show()