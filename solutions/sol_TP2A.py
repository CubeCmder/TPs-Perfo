
from velocities import *
from atmos import *
from aircraft import Aircraft


if __name__ == '__main__':

    aircraft = Aircraft()

    # CASE 1
    CD = aircraft.get_drag_coefficient()
    CDp = CD['CDp']
    CDi = CD['CDi']
    CDcomp = CD['CDcomp']
    CDctl = CD['CDctl']
    CDwm = CD['CDwm']
    CDtot = CD['CDtot']

    CL = aircraft.get_lift_coefficient()

    lift_to_drag = CL/CDtot

    thrust
    drag_force
    angle_of_attack
    







