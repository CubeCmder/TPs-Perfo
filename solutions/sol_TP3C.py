import numpy as np
import openpyxl
import re

from tabulate import tabulate

from velocities import *
from atmos import *
from aircraft import Aircraft
from climb_descent import climb_descent

aircraft = Aircraft()

