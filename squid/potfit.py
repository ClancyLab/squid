# System imports
import os
import shutil

# Squid imports
import g09
import orca
import units


def _read_orca(name):
    """
    A method for reading in the output of Orca single point calculations to
    get the atomic positions with forces.  Further, energy is also returned.

    **Parameters**

        name: *str*
            The name of the Orca simulation in questions.

    **Returns**

        new_energy: *float*
            The energy of the system in Hartree (Ha).
        new_atoms: *list,* :class:`structures.Atom`
            A list of atoms with the forces attached in units of Hartree per
            Angstrom (Ha/Ang).
    """
    read_data = orca.engrad_read(name,
                                 force='Ha/Ang',
                                 pos='Ang')
    new_atoms, new_energy = read_data

    return new_energy, new_atoms


def _read_g09(name):
    """
    A method for reading in the output of g09 single point calculations to
    get the atomic positions with forces.  Further, energy is also returned.

    **Parameters**

        name: *str*
            The name of the g09 simulation in questions.

    **Returns**

        new_energy: *float*
            The energy of the system in Hartree (Ha).
        new_atoms: *list,* :class:`structures.Atom`
            A list of atoms with the forces attached in units of Hartree per
            Angstrom (Ha/Ang).
    """
    result = g09.parse_atoms(name,
                             check_convergence=False,
                             parse_all=False)
    if not result:
        raise Exception('parse_atoms failed')
    new_energy, new_atoms = result

    for a in new_atoms:
        a.fx = units.convert('Ha/Bohr', 'Ha/Ang', a.fx)
        a.fy = units.convert('Ha/Bohr', 'Ha/Ang', a.fy)
        a.fz = units.convert('Ha/Bohr', 'Ha/Ang', a.fz)

    return new_energy, new_atoms


def _wrap(atoms, buf):
    """
    A function to get the bounding box of a system of atoms if it is not
    provided.  Buffer let's us add a small vaccum if we want along the
    x, y, and z.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            A list of atoms to get the box of.
        buf: *list, float*
            A buffer for the x, y, and z directions. This is an array
            holding only 3 values.

    **Returns**
        box: *list, list, float*
            A list holding the bounding x, y, and z dimensions.
    """
    xs = [a.x for a in atoms]
    ys = [a.y for a in atoms]
    zs = [a.z for a in atoms]

    xlen = max(xs) - min(xs) + buf[0]
    ylen = max(ys) - min(ys) + buf[1]
    zlen = max(zs) - min(zs) + buf[2]

    return [[xlen, 0.0, 0.0],
            [0.0, ylen, 0.0],
            [0.0, 0.0, zlen]]


def _unique_types(atoms):
    """
    A function to return a unique type list for the given atoms.

    **Parameters**

        atoms: *list,* :class:`structures.Atom`
            A list of atoms.

    **Returns**

        atom_types: *list, int*
            A list of types corresponding to each atom.
    """
    unique_types = {}
    for a in atoms:
        if a.element not in unique_types:
            unique_types[a.element] = len(unique_types)

    return [unique_types[a.element] for a in atoms]


def _buf_replace(buf, sid, val):
    """
    A function to find all instances of sid in buf and replace it.

    **Parameters**

        buf: *str*
            A string with instances of sid in it.
        sid: *str*
            Something to be removed from buf.
        val: *str*
            Something to take the place of sid in buf.

    **Returns**

        replaced_buf: *str*
            A string with all instances of sid replaced by val.
    """
    while str(sid) in buf:
        buf = buf.replace(str(sid), str(val))

    return buf


def generate_config(name, ts_names,
                    dft="orca", box=None, box_buffer=[10.0, 10.0, 10.0],
                    atom_types=None, n_types=False):
    """
    A function to generate a config file from a list of DFT simulations.
    The configuration file is used in potfit as the training set to
    parameterize against.

    NOTE! CURRENTLY WE ARE SETTING THE COHESIVE ENERGY TO TOTAL ENERGY! THIS
    IS NOT ACCURATE, BUT AS WE PLAN TO LEAVE eng_weight AS 0, IT DOESN'T
    MATTER!

    **Parameters**

        name: *str*
            The name to give the configuration file. Note, this is assumed to
            be exact, so no .config or .config will be added to the file name.
        ts_names: *list, str*
            A list of dft simulation names to be used for the training set.
            Note, it is HIGHLY recommended that this list be long, as we need
            sufficient information to allow for accurate training.
        dft: *str, optional*
            A string specifying the dft package used. By default this is set
            to orca, but it also allows for g09.
        box: *list, list, list, float, optional*
            A list, holding a list of lists to specify box vectors of each
            training set.  Note, if None is used then the box will be wrapped
            by the max and min atomic coordinates plus the box_buffer.
        box_buffer: *list, float, optional*
            A buffer to apply to each boundary (x, y, and z) for when box is
            set to None.
        atom_types: *list, list, int, optional*
            If you want to parameterize for varying atomic types, that is you
            want to parameterize one carbon atom as type 0 and another as type
            1, you should add in the atom_types list.  This would be a list
            for each training set, specifying the type for each atom.
        n_types: *bool, optional*
            Whether to return the number of types or not.

    **Returns**

        number_of_types: *int*
            If n_types is True, the number of types is returned
    """
    # Ensure variables passed appropriately
    dft = dft.lower()
    if (box is not None and
            len(box) != len(ts_names)):
        raise Exception("Error - If specifying box dimensions, you must \
do so for every dft output.")

    # Get read function
    if dft == "orca":
        read_ts = _read_orca
    elif dft == "g09":
        read_ts = _read_g09
    else:
        raise Exception("Error - Only allowable dft methods are orca and g09")

    # Read in training sets
    training_set = []
    energies = []
    for ts in ts_names:
        energy, atoms = read_ts(ts)
        training_set.append(atoms)
        energies.append(energy)
        if not hasattr(training_set[-1][0], 'fx'):
            raise Exception("Cannot find forces in training set %s" % ts)

    # Write out config file
    fptr_config = open(name, 'w')
    number_of_types = 0
    for i in range(len(training_set)):
        # Get variables to write out
        ts_atoms = training_set[i]
        ts_energy = energies[i]
        if box is not None:
            ts_box = box[i]
        else:
            ts_box = _wrap(ts_atoms, box_buffer)
        if atom_types is not None:
            ts_types = atom_types[i]
        else:
            ts_types = _unique_types(ts_atoms)

        if (max(ts_types) + 1) > number_of_types:
            number_of_types = max(ts_types) + 1

        # Write out header
        # Using the 1 here tells potfit to use force.
        fptr_config.write("#N %d 1\n" % len(ts_atoms))
        fptr_config.write("#X %.14e %.14e %.14e\n" % tuple(ts_box[0]))
        fptr_config.write("#Y %.14e %.14e %.14e\n" % tuple(ts_box[1]))
        fptr_config.write("#Z %.14e %.14e %.14e\n" % tuple(ts_box[2]))
        fptr_config.write("#E %.14e\n" % ts_energy)
        fptr_config.write("#F\n")

        for typ, ato in zip(ts_types, ts_atoms):
            fptr_config.write(
                "%d %.14e %.14e %.14e %.14e %.14e %.14e\n" %
                (typ, ato.x, ato.y, ato.z, ato.fx, ato.fy, ato.fz)
            )

    fptr_config.close()

    if n_types:
        return number_of_types


def run(name, ts_names, method=None,
        dft="orca", box=None, box_buffer=[10.0, 10.0, 10.0],
        atom_types=None, anneal_temp=300):
    """
    An api for potfit, allowing us to run parameterizations from python.

    NOTE! CURRENTLY WE REQUIRE THAT YOU HAVE MANUALLY MADE A POTENTIAL FILE
    IN THE CURRENT FOLDER! IF NOT THEN AN ERROR WILL BE THROWN! AS WE VALIDATE
    POTFIT FURTHER, IF IT WORKS WELL WE WILL ADD IN AUTOMATIC PARAMETER FILE
    GENERATION FROM LATIN HYPER CUBE CODE.

    **Parameters**

        name: *str*
            The name to give the generated files.
        ts_names: *list, str*
            A list of dft simulation names to be used for the training set.
            Note, it is HIGHLY recommended that this list be long, as we need
            sufficient information to allow for accurate training.
        method: *str*
            What method to parameterize for.  This needs to follow the potfit
            naming convention of potfit_method1_method2_etc. For instance,
            to do mpi (parallel parameterization) for a lennard jones system,
            you should do potfit_mpi_apot_pair.
        dft: *str, optional*
            A string specifying the dft package used. By default this is set
            to orca, but it also allows for g09.
        box: *list, list, list, float, optional*
            A list, holding a list of lists to specify box vectors of each
            training set.  Note, if None is used then the box will be wrapped
            by the max and min atomic coordinates plus the box_buffer.
        box_buffer: *list, float, optional*
            A buffer to apply to each boundary (x, y, and z) for when box is
            set to None.
        atom_types: *list, list, int, optional*
            If you want to parameterize for varying atomic types, that is you
            want to parameterize one carbon atom as type 0 and another as type
            1, you should add in the atom_types list.  This would be a list
            for each training set, specifying the type for each atom.

    **Returns**
    """

    if not os.path.exists("potfit"):
        os.mkdir("potfit")
    if not os.path.exists("potfit/%s" % name):
        os.mkdir("potfit/%s" % name)

    if not os.path.exists("%s_start.pot" % name):
        raise Exception("Could not find parameter file %s_start.pot" % name)
    else:
        shutil.copy("%s_start.pot" % name, "potfit/%s/%s_start.pot" % (name, name))

    # Generate the config file and move it into the potfit folder to run
    number_of_types = generate_config(
        name + ".config", ts_names,
        dft=dft, box=box, box_buffer=box_buffer,
        atom_types=atom_types, n_types=True)
    shutil.move(name + ".config", "potfit/%s/%s.config" % (name, name))

    # Generate the parameter file
    param_file = '''ntypes\t$NUMBER_OF_TYPES$
config\t$NAME$.config
startpot\t$NAME$_start.pot
endpot\t$NAME$_end.pot
tempfile\t$NAME$.tmp
plotfile\t$NAME$.plot
flagfile\tSTOP

write_pair\t1
write_lammps\t1

output_prefix\t$NAME$

opt\t1
anneal_temp\t$ANNEAL_TEMP$
eng_weight\t0
seed\t42'''

    param_file = _buf_replace(param_file, "$NAME$", name)
    param_file = _buf_replace(param_file, "$NUMBER_OF_TYPES$", number_of_types)
    param_file = _buf_replace(param_file, "$ANNEAL_TEMP$", anneal_temp)

    fptr_param = open("potfit/%s/%s.param" % (name, name), 'w')
    fptr_param.write(param_file)
    fptr_param.close()
