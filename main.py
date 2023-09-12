from atmos import *
from velocities import *
# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

hp = -2000
T = 35
W = 40000
p = pressure_from_alt(hp)
Vc = 150
M = get_mach_from_calibrated_airspeed(p, Vc)
V = get_true_airspeed(p,M)
Ve = get_equivalent_airspeed(p,M)
Re = get_reynolds(p,V)
mu = get_viscosity(p)

print('Vc= 150 kts')
print('mach : ' + str(M))
print('true : ' + str(V))
print('e : ' + str(Ve))
print('Re : ' + str(Re))
print('mu : ' + str(mu))

hp = 20000
T = -40
W = 40000
p = pressure_from_alt(hp)
Ve = 250
M = get_mach_from_equivalent_airspeed(p,Ve)
V = get_true_airspeed(p,M)
Vc = get_calibrated_airspeed(p,M)
Re = get_reynolds(p,V)

print('Ve= 250 kts')
print('mach : ' + str(M))
print('true : ' + str(V))
print('c : ' + str(Vc))
print('Re : ' + str(Re))
print('mu : ' + str(mu))






