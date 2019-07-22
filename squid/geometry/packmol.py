import os
from squid.utils import units
from squid.files.xyz_io import *
from squid.files.misc import which
from squid.utils.cast import is_array


def get_packmol_obj():
    '''
    This function will find the packmol executable and handle errors
    accordingly.

    **Returns**

        packmol_path: *str*
            The path to packmol.
    '''
    packmol_path = which("packmol")

    assert packmol_path is not None and packmol_path != "",\
        "Error - Unable to access Packmol.  Please ensure it is in your PATH \
environment variable!"

    return packmol_path


def packmol(system_obj, molecules, molecule_ratio=(1,),
            density=1.0, seed=1, persist=True, number=None,
            additional="", custom=None, extra_block_at_beginning='',
            extra_block_at_end='', tolerance=2.0):
    '''
    Given a list of molecules, pack this system appropriately.  Note,
    we now will pack around what is already within the system!  This is
    done by first generating a packmol block for the system at hand,
    followed by a block for the solvent.

    A custom script is also allowed; however, if this path is chosen, then
    ensure all file paths for packmol exist.  We change directories within
    this function to a sys_packmol folder, where all files are expected to
    reside.

    **Parameters**

        system_obj: :class:`structures.system.System`
            The system object to pack the molecules into.
        molecules: *list,* :class:`structures.molecule.Molecule`
            Molecules to be added to this system.
        molecule_ratio: *tuple, float, optional*
            The ration that each molecule in *molecules* will be added to
            the system.
        density: *float, optional*
            The density of the system in g/mL
        seed: *float, optional*
            Seed for random generator.
        persist: *bool, optional*
            Whether to maintain the generated sys_packmol directory or
            not.
        number: *int or list, int, optional*
            Overide density and specify the exact number of molecules to
            pack. When using a list of molecules, you must specify each
            in order within a list.
        additional: *str, optional*
            Whether to add additional constraints to the standard packmol
            setup.
        custom: *str, optional*
            A custom packmol script to run for the given input molecules.
            Note, you should ensure all necessary files are within the
            sys_packmol folder if using this option.
        extra_block_at_beginning: *str, optional*
            An additional block to put prior to the standard block.
        extra_block_at_end: *str, optional*
            An additional block to put after the standard block.
        tolerance: *float, optional*
            The tolerance around which we allow atomic overlap/proximity.

    **Returns**

        None

    **References**

        * Packmol - http://www.ime.unicamp.br/~martinez/packmol/home.shtml

    '''
    if not os.path.exists('sys_packmol'):
        os.mkdir('sys_packmol')
    os.chdir('sys_packmol')

    assert is_array(molecules),\
        "Error - Molecules must be an array!"

    packmol_path = get_packmol_obj()

    f = open(system_obj.name + '.packmol', 'w')

    if custom is not None:
        f.write(custom)
        f.close()
    else:
        f.write('''
tolerance ''' + str(tolerance) + '''
filetype xyz
output ''' + system_obj.name + '''.packed.xyz
seed ''' + str(seed) + '''
''')

        # If the system already has atoms, then set them
        if system_obj.atoms is not None and len(system_obj.atoms) > 0:
            write_xyz(system_obj.atoms, "%s_fixed.xyz" % system_obj.name)
            f.write('''
structure %s_fixed.xyz
number 1
fixed 0. 0. 0. 0. 0. 0.
centerofmass
end structure
''' % system_obj.name)

        # convert density to amu/angstrom^3. 1 g/mL = 0.6022 amu/angstrom^3
        density *= 0.6022
        average_molecular_weight = sum(
            [(units.elem_weight(a.element)) *
             molecule_ratio[i]
             for i in range(len(molecules))
             for a in molecules[i].atoms]) / sum(molecule_ratio)
        count = (density *
                 system_obj.box_size[0] *
                 system_obj.box_size[1] *
                 system_obj.box_size[2] /
                 average_molecular_weight)
        molecule_counts = [int(round(count * x / sum(molecule_ratio)))
                           for x in molecule_ratio]
        if number is not None:
            molecule_counts = number
            if type(number) is not list:
                molecule_counts = [molecule_counts]

        f.write(extra_block_at_beginning)
        lower = tuple([-x / 2.0 for x in system_obj.box_size])
        upper = tuple([x / 2.0 for x in system_obj.box_size])
        for i, m in enumerate(molecules):
            xyz_file = open('%s_%d.xyz' % (system_obj.name, i), 'w')
            xyz_file.write(str(len(m.atoms)) + '\nAtoms\n')
            for a in m.atoms:
                xyz_file.write('%s%d %f %f %f\n'
                               % (a.element, i, a.x, a.y, a.z))
            xyz_file.close()

            f.write('''
structure %s_%d.xyz
number %d
inside box %f %f %f %f %f %f''' % ((system_obj.name, i, molecule_counts[i]) + lower + upper) + '''
''' + additional + '''
end structure
''')
        f.write(extra_block_at_end)
        f.close()

    # Run packmol
    os.system(packmol_path +
              ' < ' +
              system_obj.name +
              '.packmol > packmol.log')
    atoms = read_xyz(system_obj.name + '.packed.xyz')
    os.chdir('..')

    # Now have a list of atoms with element = H0 for molecule 0,
    # H1 for molecule 1, etc
    i = 0
    offset = len(system_obj.atoms)
    while i < len(atoms):
        ints_in_element = [
            j for j, h in enumerate(atoms[i].element) if h.isdigit()]
        if len(ints_in_element) == 0:
            # This is the fixed molecule that already exists in the system
            i += offset
            continue
        # More robust, now we handle > 10 molecules!
        molecule_number = int(atoms[i].element[min(ints_in_element):])
        molecule = molecules[molecule_number]
        system_obj.add(molecule)
        # Update positions of latest molecule
        for a in system_obj.atoms[-len(molecule.atoms):]:
            a.x, a.y, a.z = atoms[i].x, atoms[i].y, atoms[i].z
            i += 1

    if not persist:
        os.system("rm -rf sys_packmol")


def run_unit_tests():
    # This unit test exists in the structures.system.System unit test.
    pass


if __name__ == "__main__":
    run_unit_tests()
