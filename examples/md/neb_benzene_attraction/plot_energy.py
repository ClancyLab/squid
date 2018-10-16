import os
from squid import units
from matplotlib import pyplot as plt


def read_energy(name, index, iteration):
    '''
    This function simply reads the energy of a given path.
    '''
    path = "lammps/%s-%d-%d/energy.profile" % (name, index, iteration)

    energy = open(path, 'r').read().strip()
    energy = open(path, 'r').read().strip().split("\n")[-1]
    return float(energy.strip().split()[-1].strip())


def get_energy_path(name, index, unit_to_use="kcal/mol"):
    '''
    This function will compile together the energies of a given index from
    MD NEB.

    **Parameters**

        name: *str*
            The NEB simulation name.
        index: *int*
            The index you wish to combine.

    **Returns**

        frames: *list, float*
            The energies with this iteration.
    '''
    # Step 1 - Find how many frames exist
    fptrs = [f for f in os.listdir("lammps") if f.startswith("%s-%d-" % (name, index))]
    if index == 0:
        fptrs = fptrs[1:-1]

    E0 = read_energy(name, 0, 0)
    EN = read_energy(name, 0, len(fptrs) + 1)

    energies = [E0]
    for iteration, _ in enumerate(fptrs):
        energies.append(read_energy(name, index, iteration + 1))
    energies.append(EN)

    return [units.convert_energy("kcal/mol", unit_to_use, e - E0) for e in energies]    


unit_to_use = "kcal/mol"
for i in range(9, 15):
    plt.plot(get_energy_path("solv_box", i, unit_to_use=unit_to_use), label=str(i))

plt.xlabel("Reaction Coordinate")
plt.ylabel("Energy (%s)" % unit_to_use)
plt.legend()
plt.show()

