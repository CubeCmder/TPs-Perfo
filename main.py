from atmos2 import *
# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
deltaISA = 0
h = 0
hp = 20000
T_C = -40
P = 0
ratio = False
atm = atm(deltaISA, h, hp, T_C, P, ratio)
print(atm.T_K)
print(atm.P)
print(atm.rho)
print(atm.deltaISA)