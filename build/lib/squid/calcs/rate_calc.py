'''
The Rate Calculation module takes themochemistry information from an Orca 
single point calculation and calculates the exponential pre-factor of the 
Arrhenius equation at a specified temperature. The algorithm uses Transition
State Theory (TST) as the basis for mathematically defining the reaction
kinetics.

- :func:`translation`
- :func:`vibration`
- :func:`rotation`
- :func:`activation_energy`
- :func:`get_rate`

------------

'''

# System imports
import sys
import math
# Squid imports
from squid import constants


def translation(molecule, T):
    '''
    Calculates the translational partition function for a reactant.

    **Parameters**

        molecule: *str*
            A string for the name of the simulation containing the
            optimized reactant.

    **Returns**

        qtrans: *float*
            The partition function for translation of a reactant.
    '''

    # Constants
    k_b = constants.K_b  # m^2kg/s^2K
    h = constants.h  # m^2kg/s
    Na = constants.Na  # mol^-1

    lines = []
    mass = 0
    f = open('orca/%s/%s.out' % (molecule, molecule), 'r')
    for line in f:
        lines.append(line.strip().split())
    f.close()

    for j in range(len(lines)):
        if len(lines[j]) > 1:
            if lines[j][0] == 'Total' and lines[j][1] == 'Mass':
                mass = float(lines[j][3]) / 1000.0
                break

    qtrans = ((2.0 * math.pi * (mass / Na) * k_b * T) / (h ** 2)) ** 1.5

    return qtrans


def vibration(molecule, T):
    '''
    Calculates the vibrational partition function for a reactant.

    **Parameters**

        molecule: *str*
            A string for the name of the simulation containing
            the optimized reactant.

    **Returns**

        qvib: *float*
            The partition function for vibration of a reactant.
    '''

    raise Exception("Error - Forgot to define c in this code. What is it?")

    # Constants
    k_b = constants.K_b  # m^2kg/s^2K
    h = constants.h  # m^2kg/s

    lines = []
    vib = []
    qvib = 1.0

    f = open('orca/%s/%s.out' % (molecule, molecule), 'r')
    for line in f:
        lines.append(line.strip().split())
    f.close()

    for j in range(len(lines)):
        if len(lines[j]) > 1:
            if lines[j][0] == 'VIBRATIONAL' and lines[j][1] == 'FREQUENCIES':
                vib_start = j + 3
            elif lines[j][0] == 'NORMAL' and lines[j][1] == 'MODES':
                vib_end = j - 3
                break

    vib_list = lines[vib_start:vib_end]

    for j in range(len(vib_list)):
        if float(vib_list[j][1]) != 0.00 and len(vib_list[j]) == 3:
            vib.append(float(vib_list[j][1]))

    vib = [c * x for x in vib]

    for j in range(len(vib)):
            qvib *= abs((1.0 / (1.0 - math.exp(-((h * vib[j]) / (k_b * T))))))

    return qvib


def rotation(molecule):
    '''
    Calculates the rotational partition function for a reactant.

    **Parameters**

        molecule: *str*
            A string for the name of the simulation containing
            the optimized reactant.

    **Returns**

        qrot: *float*
            The partition function for rotation of a reactant.
    '''

    lines = []
    qrot = 0
    rot_temp = False

    f = open('orca/%s/%s.out' % (molecule, molecule), 'r')
    for line in f:
        lines.append(line.strip().split())
    f.close()

    for j in range(len(lines)):
        if len(lines[j]) > 1:
            if lines[j][0] == 'THERMOCHEMISTRY' and lines[j][2] == '%sK' % T:
                rot_temp = True

            if lines[j][0] == 'qrot' and len(lines[j]) == 3 and rot_temp:
                qrot = float(lines[j][2]) / sn
                break

    if not rot_temp:
        print('Error - Data for input temperature not available.')
        sys.exit()

    return qrot


def activation_energy(molecule):
    '''
    Extracts the final energy of a reactant from an Orca output file.

    **Parameters**

        molecule: *str*
            A string for the name of the simulation containing
            the optimized reactant.

    **Returns**

        E_tmp: *float*
            The energy of the optimized reactant.
    '''

    lines = []
    E_tmp = 0

    f = open('orca/%s/%s.out' % (molecule, molecule), 'r')
    for line in f:
        lines.append(line.strip().split())
    f.close()

    for j in range(len(lines)):
        if len(lines[j]) > 4:
            if lines[j][0] == 'FINAL' and lines[j][1] == 'SINGLE' and lines[j][2] == 'POINT' and lines[j][3] == 'ENERGY':
                E_tmp = float(lines[j][4])
                break

    return E_tmp


def get_rate():
    '''
    Calculate the pre-exponential factor of the Arrenhius equation
    for a reaction, using Transition State Theory (TST)
    '''

    elems = [x['sym'] for x in constants.PERIODIC_TABLE[1:]]
    mw = [x['weight'] for x in constants.PERIODIC_TABLE[1:]]

    # Defaults
    (T, R_num, P_num, sn) = (296, 2, 1, 1)

    # Constants
    k_b = constants.K_b # m^2kg/s^2K
    h = constants.h # m^2kg/s
    c = 100.0 * constants.c # cm/s
    Na = constants.Na # mol^-1

    # Partition Functions
    Q_trans = [0.0 for _ in range(R_num+P_num)]
    Q_vib = [1.0 for _ in range(R_num+P_num)]
    Q_rot = [0.0 for _ in range(R_num+P_num)]

    # Simulation Names
    R, P = [], []

    # Thermodynamic Variables
    Ea = 0.0
    k = 0.0

    # User Inputs
    plotting_flags = raw_input('What is the temperature (K) for the rate calculation?: ') # K
    if plotting_flags != "":
        T = plotting_flags
    else:
        print('\nDefaulting to 296 K')
    T = '%.2f' % float(T)
    T_num = float(T)
    print(T)
    print(type(T))

    plotting_flags = raw_input('How many reactants are there?: ')
    if plotting_flags != "":
        R_num = int(plotting_flags)
    else:
        print('\nDefaulting to 2 reactants')

    plotting_flags = raw_input('How many products are there?: ')
    if plotting_flags != "":
        P_num = int(plotting_flags)
    else:
        print('\nDefaulting to 1 product')

    plotting_flags = raw_input('What is the desired rotational symmetry number? (default is 1): ')
    if plotting_flags != "":
        sn = int(plotting_flags)
    else:
        print('\nDefaulting to sn = 1')

    # Molecule Names

    i = 0
    while i < R_num:
        R.append(raw_input('What is the name of the simulation for reactant %d?: ' % (i+1)))
        i += 1

    i = 0
    while i < P_num:
        P.append(raw_input('What is the name of the simulation for product %d?: ' % (i+1)))
        i += 1

    R_comp = raw_input('What is the name of the simulation for the reactant complex?: ')

    # Calculate the rate constant
    k = ((Na * k_b * T_num) / h)

    i = 0
    while i < R_num:
        Q_trans[i] = translation(R[i], T_num)
        Q_vib[i] = vibration(R[i], T_num)
        Q_rot[i] = rotation(R[i])
        k = k / (Q_trans[i] * Q_vib[i] * Q_rot[i])
        i += 1
    while (i - R_num) < P_num:
        Q_trans[i] = translation(P[i-R_num], T_num)
        Q_vib[i] = vibration(P[i-R_num], T_num)
        Q_rot[i] = rotation(P[i-R_num])
        k = k * (Q_trans[i] * Q_vib[i] * Q_rot[i])
        i += 1

    Ea -= activation_energy(R_comp)
    i = 0
    while i < P_num:
        Ea += activation_energy(P[i])
        i += 1

    Ea = Ea * Na * constants.ENERGY['Ha']

    k = k * math.exp(-Ea / (Na * k_b * T_num))
    print('The pre-exponential factor of the reaction at T=%s is %f' % (T, k))
