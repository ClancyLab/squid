# import os
# from squid.utils.units import elem_weight


# # Imports the full atom style lammps data file.
# # Default is to import everything but you might get better performance if you
# # turn off the data you do not need
# def read_lammps_data(name, read_atoms=True, read_bonds=True,
#                      read_angles=True, read_dihedrals=True):
#     '''
#     Helper function for read_lammpstrj to read in larger files.

#     **Parameters**

#         name: *str*
#             Name of the .lammpstrj file to be read in. NOTE, this file MUST
#             have the extension .lammpstrj.
#         read_atoms: *bool, optional*
#             Whether to read in the atomic information (True), or not (False).
#         read_bonds: *bool, optional*
#             Whether to read in the bond information (True), or not (False).
#         read_angles: *bool, optional*
#             Whether to read in the angle information (True), or not (False).
#         read_dihedrals: *bool, optional*
#             Whether to read in the dihedral information (True), or not
#             (False).

#     **Returns**

#         atoms: *list,* :class:`squid.structures.atom.Atom`
#             A list of atoms read in from the data file.
#         bonds: *list,* :class:`structures.Bond`
#             A list of bonds read in from the data file.
#         angles: *list,* :class:`structures.Angle`
#             A list of angles read in from the data file.
#         dihedrals: *list,* :class:`structures.Dihedral`
#             A list of dihedrals read in from the data file.
#     '''
#     if not name.endswith('.data') and '.' not in name:
#         name += '.data'
#     # If file does not exist, return empty lammpstrj object
#     if not os.path.isfile(name):
#         raise Exception('Expected lammps data file does not exist at %s/%s'
#                         % (os.getcwd(), name))

#     # Initialize variables
#     atom_types, bond_types, angle_types, dihedral_types = [], [], [], []
#     atoms, bonds, angles, dihedrals = [], [], [], []

#     section_names = ['Masses', 'Pair Coeffs', 'Bond Coeffs', 'Angle Coeffs',
#                      'Dihedral Coeffs', 'Atoms', 'Bonds',
#                      'Angles', 'Dihedrals']
#     section_flags = {key: False for key in section_names}

#     # Iterate file line by line. This reduces the memory required since it
#     # only loads the current line
#     with open(name) as f:
#         for line in f:
#             info = line.split()

#             # Check for new section
#             if line[:-1] in section_names:
#                 for key in section_flags:
#                     section_flags[key] = False

#                 section_flags[line[:-1]] = True
#                 continue

#             if section_flags['Masses'] and len(info) > 0:
#                 atom_type = structures.Struct(index=int(info[0]),
#                                               mass=float(info[1]),
#                                               style='lj/cut')
#                 atom_types.append(atom_type)

#             if section_flags['Pair Coeffs'] and len(info) > 0:
#                 atom_type_index = int(info[0]) - 1
#                 atom_types[atom_type_index].vdw_e = float(info[1])
#                 atom_types[atom_type_index].vdw_r = float(info[2])

#             if section_flags['Bond Coeffs'] and len(info) > 0:
#                 bond_types.append(structures.Struct(e=float(info[1]),
#                                                     r=float(info[2]),
#                                                     style='harmonic'))

#             if section_flags['Angle Coeffs'] and len(info) > 0:
#                 angle_types.append(structures.Struct(e=float(info[1]),
#                                                      angle=float(info[2]),
#                                                      style='harmonic'))

#             if section_flags['Dihedral Coeffs'] and len(info) > 0:
#                 dihedral_types.append(
#                     structures.Struct(e=tuple([float(s) for s in info[1:5]]),
#                                       style='opls')
#                 )

#             if section_flags['Atoms'] and read_atoms:
#                 # Check if atom has expected number of characteristics
#                 if len(info) >= 7:
#                     # Add charge to appropriate atom_type
#                     atom_type_index = int(info[2]) - 1
#                     atom_types[atom_type_index].charge = float(info[3])

#                     # Create new atom and assign atom_type
#                     new_atom = structures.Atom(element=info[2],
#                                                x=float(info[4]),
#                                                y=float(info[5]),
#                                                z=float(info[6]),
#                                                index=int(info[0]),
#                                                bonded=[],
#                                                molecule_index=int(info[1]))
#                     new_atom.type = atom_types[atom_type_index]
#                     atoms.append(new_atom)

#                 elif len(info) > 0:
#                     print('Atom skipped due to missing information')

#             if section_flags['Bonds'] and read_bonds:
#                 # Check if bond has expected number of characteristics
#                 if len(info) == 4:
#                     a, b = int(info[2]), int(info[3])
#                     bonds.append(structures.Bond(atoms[a - 1],
#                                                  atoms[b - 1],
#                                                  type=bond_types[
#                                                  int(info[1]) - 1]
#                                                  )
#                                  )
#                     atoms[a - 1].bonded.append(atoms[b - 1])
#                     atoms[b - 1].bonded.append(atoms[a - 1])
#                 elif len(info) > 0:
#                     print('Bond skipped due to missing information')

#             if section_flags['Angles'] and read_angles:
#                 # Check if angle has expected number of characteristics
#                 if len(info) == 5:
#                     a, b, c = int(info[2]), int(info[3]), int(info[4])
#                     angles.append(structures.Angle(atoms[a - 1],
#                                   atoms[b - 1],
#                                   atoms[c - 1],
#                                   type=angle_types[int(info[1]) - 1]))
#                 elif len(info) > 0:
#                     print('Angle skipped due to missing information')

#             if section_flags['Dihedrals'] and read_dihedrals:
#                 # Check if angle has expected number of characteristics
#                 if len(info) == 6:
#                     a, b = int(info[2]), int(info[3])
#                     c, d = int(info[4]), int(info[5])
#                     dihedrals.append(
#                         structures.Dihedral(atoms[a - 1],
#                                             atoms[b - 1],
#                                             atoms[c - 1],
#                                             atoms[d - 1],
#                                             type=dihedral_types[
#                                             int(info[1]) - 1]
#                                             )
#                     )
#                 elif len(info) > 0:
#                     print('Dihedral skipped due to missing information')

#     # Create atoms and bonds
#     return atoms, bonds, angles, dihedrals


def write_lammps_data(system, **kwargs):
    '''
    Writes a lammps data file from the given system.

    **Parameters**

        system: :class:`squid.structures.system.System`
            Atomic system to be written to a lammps data file.
        pair_coeffs_included: *bool, optional*
            Whether to write pair coefficients into the data file (True),
            or not (False).

    **Returns**

        None
    '''
    pair_coeffs_included = True
    if "pair_coeffs_included" in kwargs:
        pair_coeffs_included = kwargs["pair_coeffs_included"]

    # Ensure a system name exists
    assert system.name is not None,\
        "Error - System name cannot be None!"

    # start writing file
    f = open(system.name + '.data', 'w')
    f.write('LAMMPS Description\n\n%d atoms\n%d bonds\n%d angles\n\
%d dihedrals\n0  impropers\n\n'
            % (len(system.atoms),
               len(system.bonds),
               len(system.angles),
               len(system.dihedrals)))

    f.write('%d atom types\n%d bond types\n%d angle types\n%d dihedral types\n\
0  improper types\n'
            % (len(system.atom_labels),
               len(system.parameters.bond_params),
               len(system.parameters.angle_params),
               len(system.parameters.dihedral_params)))

    f.write('%3.5f %3.5f xlo xhi\n' % (system.xlo, system.xhi))
    f.write('%3.5f %3.5f ylo yhi\n' % (system.ylo, system.yhi))
    f.write('%3.5f %3.5f zlo zhi\n' % (system.zlo, system.zhi))

    # If the system is triclinic box
    if (abs(system.box_angles[0] - 90) > 0.001 or
            abs(system.box_angles[1] - 90) > 0.001 or
            abs(system.box_angles[2] - 90) > 0.001):
        f.write('%3.5f %3.5f %3.5f xy xz yz\n'
                % (system.xy, system.xz, system.yz))

    f.write('''
Masses

''' + ('\n'.join(["%d\t%f" % (i + 1, m)
                  for i, m in enumerate(system.get_atom_masses())])) + '\n')

    if pair_coeffs_included:
        f.write('\nPair Coeffs\n\n')

        # Assume lj/cut potential since no hybrids are included
        f.write(system.dump_pair_coeffs_data())

    if system.bonds:
        f.write("\n\nBond Coeffs\n\n")
        f.write(system.parameters.dump_bonds(in_input_file=False))
    if system.angles:
        f.write("\n\nAngle Coeffs\n\n")
        f.write(system.parameters.dump_angles(in_input_file=False))
    if system.dihedrals:
        f.write("\n\nDihedral Coeffs\n\n")
        f.write(system.parameters.dump_dihedrals(in_input_file=False))

    f.write("\n\nAtoms\n\n")
    f.write(system.dump_atoms_data())
    if system.bonds:
        f.write('\n\nBonds\n\n')
        f.write(system.dump_bonds_data())
    if system.angles:
        f.write('\n\nAngles\n\n')
        f.write(system.dump_angles_data())
    if system.dihedrals:
        f.write('\n\nDihedrals\n\n')
        f.write(system.dump_dihedrals_data())
    f.write('\n\n')
    f.close()
