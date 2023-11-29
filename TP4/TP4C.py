import numpy as np

from velocities import *
from atmos import *
from aircraft import Aircraft
from tabulate import tabulate

aircraft = Aircraft()

# Donnees probleme

TO1 = aircraft.takeoff_run_velocities(45000, 0, 20, theta=0, OEI_tag=False)

print(TO1)

TO2 = aircraft.takeoff_run_velocities(30000, 2000, 0, theta=0, OEI_tag=False)

print(TO2)

TO3 = aircraft.takeoff_run_velocities(51000, 10000, 35, theta=0, OEI_tag=False)

print(TO3)
