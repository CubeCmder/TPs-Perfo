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

    hp = -2000
    T = 35+273.15
    W = 40000

    p = pressure_from_alt(hp)
    dISA = get_delta_ISA(hp,T)
    print(dISA)
    p,rho,t = get_atmos_from_dISA(hp, dISA, False)
    print(p)
    print(rho)
    print(T)
    Vc = 150
    mach = get_mach_from_calibrated_airspeed(p, Vc)
    V = get_true_airspeed(p, mach, temp=T)
    Ve = get_equivalent_airspeed(p, mach)
    Re = get_reynolds(p, V, T)
    mu = get_viscosity(p,T)
    Cl = get_lift_coefficient(p, V, W, T)
    a = get_SOS(T, knots=True)

    print(f'Vc = 150 kts')
    print(f'mach : {mach:0.5f}')
    print(f'a : {a:0.2f} kts')
    print(f'TAS : {V:0.2f} kts')
    print(f'Ve : {Ve:0.2f} kts')
    print(f'Re : {Re:0.0f}')
    print(f'mu : {mu:0.9f}')
    print(f'Cl : {Cl:0.2f}')





    hp = 20000
    T = -40+273.15
    W = 40000
    dISA = get_delta_ISA(hp,T)
    p,rho,T = get_atmos_from_dISA(hp,dISA,False)
    Ve = 250
    a = get_SOS(temp=T)
    mach = get_mach_from_equivalent_airspeed(p,Ve,rho, a)
    V = get_true_airspeed(p,mach,a)
    Vc = get_calibrated_airspeed(p,mach)
    Re = get_reynolds(p,V,T)
    Cl = get_lift_coefficient(p,V,W,T)
    mu = get_viscosity(p,T)
    print('\n\nEquivalent')
    print(f'Ve : {Ve:0.2f} kts')
    print(f'Vc : {Vc:0.2f} kts')
    print(f'mach : {mach:0.5f}')
    print(f'a : {a:0.2f} kts')
    print(f'TAS : {V:0.2f} kts')
    print(f'Re : {Re:0.0f}')
    print(f'mu : {mu:0.9f}')
    print(f'Cl : {Cl:0.2f}')

    hp = 36089
    T = -60+273.15
    W = 40000
    dISA = get_delta_ISA(hp,T)
    p,rho,T = get_atmos_from_dISA(hp,dISA,False)
    V = 450
    a = get_SOS(temp=T)
    mach = get_mach(V,a)
    Ve = get_equivalent_airspeed(p,mach,a)
    Vc = get_calibrated_airspeed(p,mach)
    Re = get_reynolds(p,V,T)
    Cl = get_lift_coefficient(p,V,W,T)
    mu = get_viscosity(p,T)
    print('\n\nTrue')
    print(f'TAS : {V:0.2f} kts')
    print(f'Ve : {Ve:0.2f} kts')
    print(f'Vc : {Vc:0.2f} kts')
    print(f'mach : {mach:0.5f}')
    print(f'a : {a:0.2f} kts')
    print(f'Re : {Re:0.0f}')
    print(f'mu : {mu:0.9f}')
    print(f'Cl : {Cl:0.2f}')

    hp = 40000
    T = -50+273.15
    W = 40000
    dISA = get_delta_ISA(hp,T)
    p,rho,T = get_atmos_from_dISA(hp,dISA,False)
    mach = 0.74
    a = get_SOS(temp=T)
    V = get_true_airspeed(p,mach,a)
    Ve = get_equivalent_airspeed(p,mach,a)
    Vc = get_calibrated_airspeed(p,mach)
    Re = get_reynolds(p,V,T)
    Cl = get_lift_coefficient(p,V,W,T)
    mu = get_viscosity(p,T)
    print('\n\nTrue')
    print(f'TAS : {V:0.2f} kts')
    print(f'Ve : {Ve:0.2f} kts')
    print(f'Vc : {Vc:0.2f} kts')
    print(f'mach : {mach:0.5f}')
    print(f'a : {a:0.2f} kts')
    print(f'Re : {Re:0.0f}')
    print(f'mu : {mu:0.9f}')
    print(f'Cl : {Cl:0.2f}')





