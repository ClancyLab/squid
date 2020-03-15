'''
Constants useful for calculations

ENERGY: *dict*
    Various units of energy in terms of Joules.
    Includes: 'Ha', 'eV', 'J', 'kcal', 'kcal/mol', 'kJ/mol', 'kT_300', 'Ry'

PRESSURE: *dict*
    Various units of pressure in terms of atmospheres.
    Includes: 'atm', 'bar', 'Pa', 'GPa'

DISTANCE: *dict*
    Various units of idstance in terms of Angstroms.
    Includes: 'Bohr', 'Ang', 'Angstrom'

K_b: *float*
    Boltzmann's Constant in Joules

h: *float*
    Plank's Constant in J*s

hbar: *float*
    Reduced Plank's Constant in J*s

amu: *float*
    Atomic mass unit in Kg

c: *float*
    Speed of light in m/s

Na: *float*
    Avogadro's Number in mol^-1

pi = PI = Pi: *float*
    The constant pi (3.141592...).

PERIODIC_TABLE: *list, dict*
    A list of dictionaries holding information pertaining to each element.
    Note, PERIODIC_TABLE[0] is a list of entries for each element's dictionary.
    This dictionary has not been fully verified, but gives appropriate values
    for rough calculations as needed (although weight should be correct for
    all elements). Contains the following:

        atomic_num: *int*
            The element's atomic number (equivalent to the index
            in `PERIODIC_TABLE`).
        weight: *float*
            The weight of the element in Atomic Mass Units.
        name: *str*
            The name of the element (such as Hydrogen).
        sym: *str*
            The symbolic name of the element (such as 'H' for Hydrogen).
        mp: *float*
            The melting point of the element in Celsius.
        bp: *float*
            The boiling point of the element in Celsius.
        density: *float*
            The density of the element in g/cm^3.
        group: *int*
            Which column of the periodic table the element resides in.
        econfig: *list, str*
            The electronic configuration of the element
            (such as ['[Ne]', '3s2', '3p5'] for Chlorine).
        ionization: *float*
            The ionization energy in eV for the element.
        vdw_r: *float*
            The van der waals radius of the element.

COLOUR = COLOR: *dict*
    A list of escape sequences for terminal output colouring on linux.
    Contains the following:

        "BLUE", "GREEN", "YELLOW", "RED", "BOLD", "UNDERLINE", and "ENDC".

------------

'''

import math
# Everything is converted to Joules
# i.e. 1 Ha = 4.359744E-18 Joules
ENERGY = {'Ha': 4.359744E-18,
          'eV': 1.602E-19,
          'J': 1.,
          'kcal': 4184. / 6.022E23,
          'kcal/mol': 4184. / 6.022E23,
          'kJ/mol': 1000. / 6.022E23,
          'kT_300': 4.14195E-21,
          'Ry': 2.1798741e-18}

# Everything is converted to atm
# i.e. 1 GPa =  9869.23 atm
PRESSURE = {'atm': 1., 'bar': 0.986923, 'Pa': 9.86923E-6, 'GPa': 9869.23}

# Everything is converted to Angstroms
# 1.e. 1 Bohr = 0.529177 Angstrom
DISTANCE = {'Bohr': 0.529177, 'Ang': 1.0, 'Angstrom': 1.0}

K_b = 1.38064852E-23  # Boltzmann's Constant in J
h = 6.62607004E-34  # Plank's Constant (from NIST) in J*s
pi = math.pi
PI = pi
Pi = pi
hbar = h / (2.0 * pi)
amu = 1.661E-27  # Kg
c = 299792458.0  # Speed of light (from NIST) in m/s
Na = 6.0221409 * 10 ** 23 # Avogadro's Number in mol^-1

# Constants for linux colors
COLOUR = {"BLUE": '\033[94m',
          "GREEN": '\033[92m',
          "YELLOW": '\033[93m',
          "RED": '\033[91m',
          "BOLD": '\033[1m',
          "UNDERLINE": '\033[4m',
          "ENDC": '\033[0m'}
COLOR = COLOUR

# OPTIMIZER CONSTANTS
STEP_SIZE_TOO_SMALL = -2
FAIL_CONVERGENCE = -1
MAXITER_CONVERGENCE = 0
G_MAX_CONVERGENCE = 1
G_RMS_CONVERGENCE = 2

# Periodic Table data from http://www.science.co.il/PTelements.asp
# Van der Waals Radii from
#   http://link.springer.com/article/10.1023%2FA%3A1011625728803
PERIODIC_TABLE = [['atomic_num', 'weight (AMU)', 'name', 'sym', 'mp (C)', 'bp (C)', 'density (g/cm3)', 'group', 'econfig', 'ionization (eV)', 'vdw_r (ang)'],
                  {'atomic_num': 1, 'weight': 1.0079, 'sym': 'H', 'ionization': 13.5984, 'bp': -253.0, 'group': 1, 'name': 'Hydrogen', 'density': 0.09, 'vdw_r': None, 'mp': -259.0, 'econfig': ['1s1']},
                  {'atomic_num': 2, 'weight': 4.0026, 'sym': 'He', 'ionization': 24.5874, 'bp': -269.0, 'group': 18, 'name': 'Helium', 'density': 0.18, 'vdw_r': None, 'mp': -272.0, 'econfig': ['1s2']},
                  {'atomic_num': 3, 'weight': 6.941, 'sym': 'Li', 'ionization': 5.3917, 'bp': 1347.0, 'group': 1, 'name': 'Lithium', 'density': 0.53, 'vdw_r': 2.63, 'mp': 180.0, 'econfig': ['[He]', '2s1']},
                  {'atomic_num': 4, 'weight': 9.0122, 'sym': 'Be', 'ionization': 9.3227, 'bp': 2970.0, 'group': 2, 'name': 'Beryllium', 'density': 1.85, 'vdw_r': 2.23, 'mp': 1278.0, 'econfig': ['[He]', '2s2']},
                  {'atomic_num': 5, 'weight': 10.811, 'sym': 'B', 'ionization': 8.298, 'bp': 2550.0, 'group': 13, 'name': 'Boron', 'density': 2.34, 'vdw_r': 2.05, 'mp': 2300.0, 'econfig': ['[He]', '2s2', '2p1']},
                  {'atomic_num': 6, 'weight': 12.0107, 'sym': 'C', 'ionization': 11.2603, 'bp': 4827.0, 'group': 14, 'name': 'Carbon', 'density': 2.26, 'vdw_r': 1.96, 'mp': 3500.0, 'econfig': ['[He]', '2s2', '2p2']},
                  {'atomic_num': 7, 'weight': 14.0067, 'sym': 'N', 'ionization': 14.5341, 'bp': -196.0, 'group': 15, 'name': 'Nitrogen', 'density': 1.25, 'vdw_r': 1.79, 'mp': -210.0, 'econfig': ['[He]', '2s2', '2p3']},
                  {'atomic_num': 8, 'weight': 15.9994, 'sym': 'O', 'ionization': 13.6181, 'bp': -183.0, 'group': 16, 'name': 'Oxygen', 'density': 1.43, 'vdw_r': 1.71, 'mp': -218.0, 'econfig': ['[He]', '2s2', '2p4']},
                  {'atomic_num': 9, 'weight': 18.9984, 'sym': 'F', 'ionization': 17.4228, 'bp': -188.0, 'group': 17, 'name': 'Fluorine', 'density': 1.7, 'vdw_r': 1.65, 'mp': -220.0, 'econfig': ['[He]', '2s2', '2p5']},
                  {'atomic_num': 10, 'weight': 20.1797, 'sym': 'Ne', 'ionization': 21.5645, 'bp': -246.0, 'group': 18, 'name': 'Neon', 'density': 0.9, 'vdw_r': None, 'mp': -249.0, 'econfig': ['[He]', '2s2', '2p6']},
                  {'atomic_num': 11, 'weight': 22.9897, 'sym': 'Na', 'ionization': 5.1391, 'bp': 883.0, 'group': 1, 'name': 'Sodium', 'density': 0.97, 'vdw_r': 2.77, 'mp': 98.0, 'econfig': ['[Ne]', '3s1']},
                  {'atomic_num': 12, 'weight': 24.305, 'sym': 'Mg', 'ionization': 7.6462, 'bp': 1090.0, 'group': 2, 'name': 'Magnesium', 'density': 1.74, 'vdw_r': 2.42, 'mp': 639.0, 'econfig': ['[Ne]', '3s2']},
                  {'atomic_num': 13, 'weight': 26.9815, 'sym': 'Al', 'ionization': 5.9858, 'bp': 2467.0, 'group': 13, 'name': 'Aluminum', 'density': 2.7, 'vdw_r': 2.4, 'mp': 660.0, 'econfig': ['[Ne]', '3s2', '3p1']},
                  {'atomic_num': 14, 'weight': 28.0855, 'sym': 'Si', 'ionization': 8.1517, 'bp': 2355.0, 'group': 14, 'name': 'Silicon', 'density': 2.33, 'vdw_r': 2.26, 'mp': 1410.0, 'econfig': ['[Ne]', '3s2', '3p2']},
                  {'atomic_num': 15, 'weight': 30.9738, 'sym': 'P', 'ionization': 10.4867, 'bp': 280.0, 'group': 15, 'name': 'Phosphorus', 'density': 1.82, 'vdw_r': 2.14, 'mp': 44.0, 'econfig': ['[Ne]', '3s2', '3p3']},
                  {'atomic_num': 16, 'weight': 32.065, 'sym': 'S', 'ionization': 10.36, 'bp': 445.0, 'group': 16, 'name': 'Sulfur', 'density': 2.07, 'vdw_r': 2.06, 'mp': 113.0, 'econfig': ['[Ne]', '3s2', '3p4']},
                  {'atomic_num': 17, 'weight': 35.453, 'sym': 'Cl', 'ionization': 12.9676, 'bp': -35.0, 'group': 17, 'name': 'Chlorine', 'density': 3.21, 'vdw_r': 2.05, 'mp': -101.0, 'econfig': ['[Ne]', '3s2', '3p5']},
                  {'atomic_num': 18, 'weight': 39.948, 'sym': 'Ar', 'ionization': 15.7596, 'bp': -186.0, 'group': 18, 'name': 'Argon', 'density': 1.78, 'vdw_r': None, 'mp': -189.0, 'econfig': ['[Ne]', '3s2', '3p6']},
                  {'atomic_num': 19, 'weight': 39.0983, 'sym': 'K', 'ionization': 4.3407, 'bp': 774.0, 'group': 1, 'name': 'Potassium', 'density': 0.86, 'vdw_r': 3.02, 'mp': 64.0, 'econfig': ['[Ar]', '4s1']},
                  {'atomic_num': 20, 'weight': 40.078, 'sym': 'Ca', 'ionization': 6.1132, 'bp': 1484.0, 'group': 2, 'name': 'Calcium', 'density': 1.55, 'vdw_r': 2.78, 'mp': 839.0, 'econfig': ['[Ar]', '4s2']},
                  {'atomic_num': 21, 'weight': 44.9559, 'sym': 'Sc', 'ionization': 6.5615, 'bp': 2832.0, 'group': 3, 'name': 'Scandium', 'density': 2.99, 'vdw_r': 2.62, 'mp': 1539.0, 'econfig': ['[Ar]', '3d1', '4s2']},
                  {'atomic_num': 22, 'weight': 47.867, 'sym': 'Ti', 'ionization': 6.8281, 'bp': 3287.0, 'group': 4, 'name': 'Titanium', 'density': 4.54, 'vdw_r': 2.44, 'mp': 1660.0, 'econfig': ['[Ar]', '3d2', '4s2']},
                  {'atomic_num': 23, 'weight': 50.9415, 'sym': 'V', 'ionization': 6.7462, 'bp': 3380.0, 'group': 5, 'name': 'Vanadium', 'density': 6.11, 'vdw_r': 2.27, 'mp': 1890.0, 'econfig': ['[Ar]', '3d3', '4s2']},
                  {'atomic_num': 24, 'weight': 51.9961, 'sym': 'Cr', 'ionization': 6.7665, 'bp': 2672.0, 'group': 6, 'name': 'Chromium', 'density': 7.19, 'vdw_r': 2.23, 'mp': 1857.0, 'econfig': ['[Ar]', '3d5', '4s1']},
                  {'atomic_num': 25, 'weight': 54.938, 'sym': 'Mn', 'ionization': 7.434, 'bp': 1962.0, 'group': 7, 'name': 'Manganese', 'density': 7.43, 'vdw_r': 2.25, 'mp': 1245.0, 'econfig': ['[Ar]', '3d5', '4s2']},
                  {'atomic_num': 26, 'weight': 55.845, 'sym': 'Fe', 'ionization': 7.9024, 'bp': 2750.0, 'group': 8, 'name': 'Iron', 'density': 7.87, 'vdw_r': 2.27, 'mp': 1535.0, 'econfig': ['[Ar]', '3d6', '4s2']},
                  {'atomic_num': 27, 'weight': 58.9332, 'sym': 'Co', 'ionization': 7.881, 'bp': 2870.0, 'group': 9, 'name': 'Cobalt', 'density': 8.9, 'vdw_r': 2.25, 'mp': 1495.0, 'econfig': ['[Ar]', '3d7', '4s2']},
                  {'atomic_num': 28, 'weight': 58.6934, 'sym': 'Ni', 'ionization': 7.6398, 'bp': 2732.0, 'group': 10, 'name': 'Nickel', 'density': 8.9, 'vdw_r': 2.23, 'mp': 1453.0, 'econfig': ['[Ar]', '3d8', '4s2']},
                  {'atomic_num': 29, 'weight': 63.546, 'sym': 'Cu', 'ionization': 7.7264, 'bp': 2567.0, 'group': 11, 'name': 'Copper', 'density': 8.96, 'vdw_r': 2.27, 'mp': 1083.0, 'econfig': ['[Ar]', '3d10', '4s1']},
                  {'atomic_num': 30, 'weight': 65.39, 'sym': 'Zn', 'ionization': 9.3942, 'bp': 907.0, 'group': 12, 'name': 'Zinc', 'density': 7.13, 'vdw_r': 2.24, 'mp': 420.0, 'econfig': ['[Ar]', '3d10', '4s2']},
                  {'atomic_num': 31, 'weight': 69.723, 'sym': 'Ga', 'ionization': 5.9993, 'bp': 2403.0, 'group': 13, 'name': 'Gallium', 'density': 5.91, 'vdw_r': 2.41, 'mp': 30.0, 'econfig': ['[Ar]', '3d10', '4s2', '4p1']},
                  {'atomic_num': 32, 'weight': 72.64, 'sym': 'Ge', 'ionization': 7.8994, 'bp': 2830.0, 'group': 14, 'name': 'Germanium', 'density': 5.32, 'vdw_r': 2.32, 'mp': 937.0, 'econfig': ['[Ar]', '3d10', '4s2', '4p2']},
                  {'atomic_num': 33, 'weight': 74.9216, 'sym': 'As', 'ionization': 9.7886, 'bp': 613.0, 'group': 15, 'name': 'Arsenic', 'density': 5.72, 'vdw_r': 2.25, 'mp': 81.0, 'econfig': ['[Ar]', '3d10', '4s2', '4p3']},
                  {'atomic_num': 34, 'weight': 78.96, 'sym': 'Se', 'ionization': 9.7524, 'bp': 685.0, 'group': 16, 'name': 'Selenium', 'density': 4.79, 'vdw_r': 2.18, 'mp': 217.0, 'econfig': ['[Ar]', '3d10', '4s2', '4p4']},
                  {'atomic_num': 35, 'weight': 79.904, 'sym': 'Br', 'ionization': 11.8138, 'bp': 59.0, 'group': 17, 'name': 'Bromine', 'density': 3.12, 'vdw_r': 2.1, 'mp': -7.0, 'econfig': ['[Ar]', '3d10', '4s2', '4p5']},
                  {'atomic_num': 36, 'weight': 83.8, 'sym': 'Kr', 'ionization': 13.9996, 'bp': -153.0, 'group': 18, 'name': 'Krypton', 'density': 3.75, 'vdw_r': None, 'mp': -157.0, 'econfig': ['[Ar]', '3d10', '4s2', '4p6']},
                  {'atomic_num': 37, 'weight': 85.4678, 'sym': 'Rb', 'ionization': 4.1771, 'bp': 688.0, 'group': 1, 'name': 'Rubidium', 'density': 1.63, 'vdw_r': 3.15, 'mp': 39.0, 'econfig': ['[Kr]', '5s1']},
                  {'atomic_num': 38, 'weight': 87.62, 'sym': 'Sr', 'ionization': 5.6949, 'bp': 1384.0, 'group': 2, 'name': 'Strontium', 'density': 2.54, 'vdw_r': 2.94, 'mp': 769.0, 'econfig': ['[Kr]', '5s2']},
                  {'atomic_num': 39, 'weight': 88.9059, 'sym': 'Y', 'ionization': 6.2173, 'bp': 3337.0, 'group': 3, 'name': 'Yttrium', 'density': 4.47, 'vdw_r': 2.71, 'mp': 1523.0, 'econfig': ['[Kr]', '4d1', '5s2']},
                  {'atomic_num': 40, 'weight': 91.224, 'sym': 'Zr', 'ionization': 6.6339, 'bp': 4377.0, 'group': 4, 'name': 'Zirconium', 'density': 6.51, 'vdw_r': 2.57, 'mp': 1852.0, 'econfig': ['[Kr]', '4d2', '5s2']},
                  {'atomic_num': 41, 'weight': 92.9064, 'sym': 'Nb', 'ionization': 6.7589, 'bp': 4927.0, 'group': 5, 'name': 'Niobium', 'density': 8.57, 'vdw_r': 2.46, 'mp': 2468.0, 'econfig': ['[Kr]', '4d4', '5s1']},
                  {'atomic_num': 42, 'weight': 95.94, 'sym': 'Mo', 'ionization': 7.0924, 'bp': 4612.0, 'group': 6, 'name': 'Molybdenum', 'density': 10.22, 'vdw_r': 2.39, 'mp': 2617.0, 'econfig': ['[Kr]', '4d5', '5s1']},
                  {'atomic_num': 43, 'weight': 98.0, 'sym': 'Tc', 'ionization': 7.28, 'bp': 4877.0, 'group': 7, 'name': 'Technetium', 'density': 11.5, 'vdw_r': 2.37, 'mp': 2200.0, 'econfig': ['[Kr]', '4d5', '5s2']},
                  {'atomic_num': 44, 'weight': 101.07, 'sym': 'Ru', 'ionization': 7.3605, 'bp': 3900.0, 'group': 8, 'name': 'Ruthenium', 'density': 12.37, 'vdw_r': 2.37, 'mp': 2250.0, 'econfig': ['[Kr]', '4d7', '5s1']},
                  {'atomic_num': 45, 'weight': 102.9055, 'sym': 'Rh', 'ionization': 7.4589, 'bp': 3727.0, 'group': 9, 'name': 'Rhodium', 'density': 12.41, 'vdw_r': 2.32, 'mp': 1966.0, 'econfig': ['[Kr]', '4d8', '5s1']},
                  {'atomic_num': 46, 'weight': 106.42, 'sym': 'Pd', 'ionization': 8.3369, 'bp': 2927.0, 'group': 10, 'name': 'Palladium', 'density': 12.02, 'vdw_r': 2.35, 'mp': 1552.0, 'econfig': ['[Kr]', '4d10']},
                  {'atomic_num': 47, 'weight': 107.8682, 'sym': 'Ag', 'ionization': 7.5762, 'bp': 2212.0, 'group': 11, 'name': 'Silver', 'density': 10.5, 'vdw_r': 2.37, 'mp': 962.0, 'econfig': ['[Kr]', '4d10', '5s1']},
                  {'atomic_num': 48, 'weight': 112.411, 'sym': 'Cd', 'ionization': 8.9938, 'bp': 765.0, 'group': 12, 'name': 'Cadmium', 'density': 8.65, 'vdw_r': 2.37, 'mp': 321.0, 'econfig': ['[Kr]', '4d10', '5s2']},
                  {'atomic_num': 49, 'weight': 114.818, 'sym': 'In', 'ionization': 5.7864, 'bp': 2000.0, 'group': 13, 'name': 'Indium', 'density': 7.31, 'vdw_r': 2.53, 'mp': 157.0, 'econfig': ['[Kr]', '4d10', '5s2', '5p1']},
                  {'atomic_num': 50, 'weight': 118.71, 'sym': 'Sn', 'ionization': 7.3439, 'bp': 2270.0, 'group': 14, 'name': 'Tin', 'density': 7.31, 'vdw_r': 2.46, 'mp': 232.0, 'econfig': ['[Kr]', '4d10', '5s2', '5p2']},
                  {'atomic_num': 51, 'weight': 121.76, 'sym': 'Sb', 'ionization': 8.6084, 'bp': 1750.0, 'group': 15, 'name': 'Antimony', 'density': 6.68, 'vdw_r': 2.41, 'mp': 630.0, 'econfig': ['[Kr]', '4d10', '5s2', '5p3']},
                  {'atomic_num': 52, 'weight': 127.6, 'sym': 'Te', 'ionization': 9.0096, 'bp': 990.0, 'group': 16, 'name': 'Tellurium', 'density': 6.24, 'vdw_r': 2.36, 'mp': 449.0, 'econfig': ['[Kr]', '4d10', '5s2', '5p4']},
                  {'atomic_num': 53, 'weight': 126.9045, 'sym': 'I', 'ionization': 10.4513, 'bp': 184.0, 'group': 17, 'name': 'Iodine', 'density': 4.93, 'vdw_r': 2.22, 'mp': 114.0, 'econfig': ['[Kr]', '4d10', '5s2', '5p5']},
                  {'atomic_num': 54, 'weight': 131.293, 'sym': 'Xe', 'ionization': 12.1298, 'bp': -108.0, 'group': 18, 'name': 'Xenon', 'density': 5.9, 'vdw_r': None, 'mp': -112.0, 'econfig': ['[Kr]', '4d10', '5s2', '5p6']},
                  {'atomic_num': 55, 'weight': 132.9055, 'sym': 'Cs', 'ionization': 3.8939, 'bp': 678.0, 'group': 1, 'name': 'Cesium', 'density': 1.87, 'vdw_r': 3.3, 'mp': 29.0, 'econfig': ['[Xe]', '6s1']},
                  {'atomic_num': 56, 'weight': 137.327, 'sym': 'Ba', 'ionization': 5.2117, 'bp': 1140.0, 'group': 2, 'name': 'Barium', 'density': 3.59, 'vdw_r': 3.05, 'mp': 725.0, 'econfig': ['[Xe]', '6s2']},
                  {'atomic_num': 57, 'weight': 138.9055, 'sym': 'La', 'ionization': 5.5769, 'bp': 3469.0, 'group': 3, 'name': 'Lanthanum', 'density': 6.15, 'vdw_r': 2.81, 'mp': 920.0, 'econfig': ['[Xe]', '5d1', '6s2']},
                  {'atomic_num': 58, 'weight': 140.116, 'sym': 'Ce', 'ionization': 5.5387, 'bp': 3257.0, 'group': 101, 'name': 'Cerium', 'density': 6.77, 'vdw_r': None, 'mp': 795.0, 'econfig': ['[Xe]', '4f1', '5d1', '6s2']},
                  {'atomic_num': 59, 'weight': 140.9077, 'sym': 'Pr', 'ionization': 5.473, 'bp': 3127.0, 'group': 101, 'name': 'Praseodymium', 'density': 6.77, 'vdw_r': None, 'mp': 935.0, 'econfig': ['[Xe]', '4f3', '6s2']},
                  {'atomic_num': 60, 'weight': 144.24, 'sym': 'Nd', 'ionization': 5.525, 'bp': 3127.0, 'group': 101, 'name': 'Neodymium', 'density': 7.01, 'vdw_r': None, 'mp': 1010.0, 'econfig': ['[Xe]', '4f4', '6s2']},
                  {'atomic_num': 61, 'weight': 145.0, 'sym': 'Pm', 'ionization': 5.582, 'bp': 3000.0, 'group': 101, 'name': 'Promethium', 'density': 7.3, 'vdw_r': None, 'mp': 1100.0, 'econfig': ['[Xe]', '4f5', '6s2']},
                  {'atomic_num': 62, 'weight': 150.36, 'sym': 'Sm', 'ionization': 5.6437, 'bp': 1900.0, 'group': 101, 'name': 'Samarium', 'density': 7.52, 'vdw_r': None, 'mp': 1072.0, 'econfig': ['[Xe]', '4f6', '6s2']},
                  {'atomic_num': 63, 'weight': 151.964, 'sym': 'Eu', 'ionization': 5.6704, 'bp': 1597.0, 'group': 101, 'name': 'Europium', 'density': 5.24, 'vdw_r': None, 'mp': 822.0, 'econfig': ['[Xe]', '4f7', '6s2']},
                  {'atomic_num': 64, 'weight': 157.25, 'sym': 'Gd', 'ionization': 6.1501, 'bp': 3233.0, 'group': 101, 'name': 'Gadolinium', 'density': 7.9, 'vdw_r': None, 'mp': 1311.0, 'econfig': ['[Xe]', '4f7', '5d1', '6s2']},
                  {'atomic_num': 65, 'weight': 158.9253, 'sym': 'Tb', 'ionization': 5.8638, 'bp': 3041.0, 'group': 101, 'name': 'Terbium', 'density': 8.23, 'vdw_r': None, 'mp': 1360.0, 'econfig': ['[Xe]', '4f9', '6s2']},
                  {'atomic_num': 66, 'weight': 162.5, 'sym': 'Dy', 'ionization': 5.9389, 'bp': 2562.0, 'group': 101, 'name': 'Dysprosium', 'density': 8.55, 'vdw_r': None, 'mp': 1412.0, 'econfig': ['[Xe]', '4f10', '6s2']},
                  {'atomic_num': 67, 'weight': 164.9303, 'sym': 'Ho', 'ionization': 6.0215, 'bp': 2720.0, 'group': 101, 'name': 'Holmium', 'density': 8.8, 'vdw_r': None, 'mp': 1470.0, 'econfig': ['[Xe]', '4f11', '6s2']},
                  {'atomic_num': 68, 'weight': 167.259, 'sym': 'Er', 'ionization': 6.1077, 'bp': 2510.0, 'group': 101, 'name': 'Erbium', 'density': 9.07, 'vdw_r': None, 'mp': 1522.0, 'econfig': ['[Xe]', '4f12', '6s2']},
                  {'atomic_num': 69, 'weight': 168.9342, 'sym': 'Tm', 'ionization': 6.1843, 'bp': 1727.0, 'group': 101, 'name': 'Thulium', 'density': 9.32, 'vdw_r': None, 'mp': 1545.0, 'econfig': ['[Xe]', '4f13', '6s2']},
                  {'atomic_num': 70, 'weight': 173.04, 'sym': 'Yb', 'ionization': 6.2542, 'bp': 1466.0, 'group': 101, 'name': 'Ytterbium', 'density': 6.9, 'vdw_r': None, 'mp': 824.0, 'econfig': ['[Xe]', '4f14', '6s2']},
                  {'atomic_num': 71, 'weight': 174.967, 'sym': 'Lu', 'ionization': 5.4259, 'bp': 3315.0, 'group': 101, 'name': 'Lutetium', 'density': 9.84, 'vdw_r': None, 'mp': 1656.0, 'econfig': ['[Xe]', '4f14', '5d1', '6s2']},
                  {'atomic_num': 72, 'weight': 178.49, 'sym': 'Hf', 'ionization': 6.8251, 'bp': 5400.0, 'group': 4, 'name': 'Hafnium', 'density': 13.31, 'vdw_r': 2.52, 'mp': 2150.0, 'econfig': ['[Xe]', '4f14', '5d2', '6s2']},
                  {'atomic_num': 73, 'weight': 180.9479, 'sym': 'Ta', 'ionization': 7.5496, 'bp': 5425.0, 'group': 5, 'name': 'Tantalum', 'density': 16.65, 'vdw_r': 2.42, 'mp': 2996.0, 'econfig': ['[Xe]', '4f14', '5d3', '6s2']},
                  {'atomic_num': 74, 'weight': 183.84, 'sym': 'W', 'ionization': 7.864, 'bp': 5660.0, 'group': 6, 'name': 'Tungsten', 'density': 19.35, 'vdw_r': 2.36, 'mp': 3410.0, 'econfig': ['[Xe]', '4f14', '5d4', '6s2']},
                  {'atomic_num': 75, 'weight': 186.207, 'sym': 'Re', 'ionization': 7.8335, 'bp': 5627.0, 'group': 7, 'name': 'Rhenium', 'density': 21.04, 'vdw_r': 2.35, 'mp': 3180.0, 'econfig': ['[Xe]', '4f14', '5d5', '6s2']},
                  {'atomic_num': 76, 'weight': 190.23, 'sym': 'Os', 'ionization': 8.4382, 'bp': 5027.0, 'group': 8, 'name': 'Osmium', 'density': 22.6, 'vdw_r': 2.33, 'mp': 3045.0, 'econfig': ['[Xe]', '4f14', '5d6', '6s2']},
                  {'atomic_num': 77, 'weight': 192.217, 'sym': 'Ir', 'ionization': 8.967, 'bp': 4527.0, 'group': 9, 'name': 'Iridium', 'density': 22.4, 'vdw_r': 2.34, 'mp': 2410.0, 'econfig': ['[Xe]', '4f14', '5d7', '6s2']},
                  {'atomic_num': 78, 'weight': 195.078, 'sym': 'Pt', 'ionization': 8.9587, 'bp': 3827.0, 'group': 10, 'name': 'Platinum', 'density': 21.45, 'vdw_r': 2.37, 'mp': 1772.0, 'econfig': ['[Xe]', '4f14', '5d9', '6s1']},
                  {'atomic_num': 79, 'weight': 196.9665, 'sym': 'Au', 'ionization': 9.2255, 'bp': 2807.0, 'group': 11, 'name': 'Gold', 'density': 19.32, 'vdw_r': 2.41, 'mp': 1064.0, 'econfig': ['[Xe]', '4f14', '5d10', '6s1']},
                  {'atomic_num': 80, 'weight': 200.59, 'sym': 'Hg', 'ionization': 10.4375, 'bp': 357.0, 'group': 12, 'name': 'Mercury', 'density': 13.55, 'vdw_r': 2.25, 'mp': -39.0, 'econfig': ['[Xe]', '4f14', '5d10', '6s2']},
                  {'atomic_num': 81, 'weight': 204.3833, 'sym': 'Tl', 'ionization': 6.1082, 'bp': 1457.0, 'group': 13, 'name': 'Thallium', 'density': 11.85, 'vdw_r': 2.53, 'mp': 303.0, 'econfig': ['[Xe]', '4f14', '5d10', '6s2', '6p1']},
                  {'atomic_num': 82, 'weight': 207.2, 'sym': 'Pb', 'ionization': 7.4167, 'bp': 1740.0, 'group': 14, 'name': 'Lead', 'density': 11.35, 'vdw_r': 2.53, 'mp': 327.0, 'econfig': ['[Xe]', '4f14', '5d10', '6s2', '6p2']},
                  {'atomic_num': 83, 'weight': 208.9804, 'sym': 'Bi', 'ionization': 7.2856, 'bp': 1560.0, 'group': 15, 'name': 'Bismuth', 'density': 9.75, 'vdw_r': 3.52, 'mp': 271.0, 'econfig': ['[Xe]', '4f14', '5d10', '6s2', '6p3']},
                  {'atomic_num': 84, 'weight': 209.0, 'sym': 'Po', 'ionization': 8.417, 'bp': 962.0, 'group': 16, 'name': 'Polonium', 'density': 9.3, 'vdw_r': None, 'mp': 254.0, 'econfig': ['[Xe]', '4f14', '5d10', '6s2', '6p4']},
                  {'atomic_num': 85, 'weight': 210.0, 'sym': 'At', 'ionization': 9.3, 'bp': 337.0, 'group': 17, 'name': 'Astatine', 'density': None, 'vdw_r': None, 'mp': 302.0, 'econfig': ['[Xe]', '4f14', '5d10', '6s2', '6p5']},
                  {'atomic_num': 86, 'weight': 222.0, 'sym': 'Rn', 'ionization': 10.7485, 'bp': -62.0, 'group': 18, 'name': 'Radon', 'density': 9.73, 'vdw_r': None, 'mp': -71.0, 'econfig': ['[Xe]', '4f14', '5d10', '6s2', '6p6']},
                  {'atomic_num': 87, 'weight': 223.0, 'sym': 'Fr', 'ionization': 4.0727, 'bp': 677.0, 'group': 1, 'name': 'Francium', 'density': None, 'vdw_r': None, 'mp': 27.0, 'econfig': ['[Rn]', '7s1']},
                  {'atomic_num': 88, 'weight': 226.0, 'sym': 'Ra', 'ionization': 5.2784, 'bp': 1737.0, 'group': 2, 'name': 'Radium', 'density': 5.5, 'vdw_r': None, 'mp': 700.0, 'econfig': ['[Rn]', '7s2']},
                  {'atomic_num': 89, 'weight': 227.0, 'sym': 'Ac', 'ionization': 5.17, 'bp': 3200.0, 'group': 3, 'name': 'Actinium', 'density': 10.07, 'vdw_r': None, 'mp': 1050.0, 'econfig': ['[Rn]', '6d1', '7s2']},
                  {'atomic_num': 90, 'weight': 232.0381, 'sym': 'Th', 'ionization': 6.3067, 'bp': 4790.0, 'group': 102, 'name': 'Thorium', 'density': 11.72, 'vdw_r': 2.75, 'mp': 1750.0, 'econfig': ['[Rn]', '6d2', '7s2']},
                  {'atomic_num': 91, 'weight': 231.0359, 'sym': 'Pa', 'ionization': 5.89, 'bp': None, 'group': 102, 'name': 'Protactinium', 'density': 15.4, 'vdw_r': None, 'mp': 1568.0, 'econfig': ['[Rn]', '5f2', '6d1', '7s2']},
                  {'atomic_num': 92, 'weight': 238.0289, 'sym': 'U', 'ionization': 6.1941, 'bp': 3818.0, 'group': 102, 'name': 'Uranium', 'density': 18.95, 'vdw_r': 2.65, 'mp': 1132.0, 'econfig': ['[Rn]', '5f3', '6d1', '7s2']},
                  {'atomic_num': 93, 'weight': None, 'sym': None, 'ionization': None, 'bp': None, 'group': None, 'name': None, 'density': None, 'vdw_r': None, 'mp': None, 'econfig': None},
                  {'atomic_num': 94, 'weight': None, 'sym': None, 'ionization': None, 'bp': None, 'group': None, 'name': None, 'density': None, 'vdw_r': None, 'mp': None, 'econfig': None},
                  {'atomic_num': 95, 'weight': 243.0, 'sym': 'Am', 'ionization': None, 'bp': 2284.0, 'group': 102, 'name': 'Americium', 'density': 12, 'vdw_r': None, 'mp': 1449.0, 'econfig': ['[Rn]', '5f7', '7s2']},
]
