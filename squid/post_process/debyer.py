'''
Python hooks for the debyer code.
Link: https://debyer.readthedocs.io/en/latest/

- :func:`get_pdf`

------------

'''

import os
import random

from squid.files.misc import which
from squid.files.xyz_io import write_xyz


def get_pdf(frames, start=0.0, stop=5.0, step=0.1, cutoff=10.0,
            rho=1.0, quanta=0.001, output=None, persist=False):
    '''
    Obtain the pair distribution function of a list of atoms using
    the Debyer code.

    **Parameters**

        frames: *str or list,* :class:`squid.structures.atom.Atom`
            An xyz file name (with or without the .xyz extension) or
            an input frame to calculate the pdf for.
        start: *float, optional*
            The starting radial distance in Angstroms for the
            calculated pattern.
        stop: *float, optional*
            The ending radial distance in Angstroms for the calculated pattern.
        step: *float, optional*
            Step in Angstroms for the calculated pattern.
        cutoff: *float, optional*
            Cutoff distance in Angstroms for Interatomic Distance
            (ID) calculations.
        rho: *float, optional*
            Numeric density of the system.
        quanta: *float, optional*
            Interatomic Distance (ID) discritization quanta.
        output: *str, optional*
            Output file name with NO extension given
        persist: *bool, optional*
            Whether to persist made .g and .xyz files (True), or remove
            them (False)

    **Returns**

        pdf: *list, tuple, float*
            A list of tuples holding the pdf data (distance in
            Angstroms and Intensity).

    **References**

        * https://debyer.readthedocs.io/en/latest/
    '''

    debyer_path = which("debyer")
    assert debyer_path is not None,\
        "Error - Cannot find debyer in the PATH env variable."

    # If passed frames and not an xyz file name, write to xyz
    append = str(int(random.random() * 1E12))
    if type(frames) is not str:
        write_xyz(frames, "tmp_for_pdf_%s" % append)
        file_name = "tmp_for_pdf_%s" % append
    else:
        file_name = frames

    # Else, we want to ensure file_name is correct
    if file_name.endswith(".xyz"):
        file_name = file_name.split(".xyz")[0]
    if output is None:
        output = file_name

    if stop > cutoff:
        raise Exception("Stopping position should be larger less than \
or equal to the cutoff.")

    # Make command for debyer
    cmd = debyer_path
    cmd += " --cutoff=%.2f" % cutoff
    cmd += " --quanta=%.2f" % quanta
    cmd += " -g"
    cmd += " -f%.2f" % start
    cmd += " -t%.2f" % stop
    cmd += " -s%.2f" % step
    cmd += " --ro=%.2f" % rho
    cmd += " -o %s.g" % output
    cmd += " %s.xyz" % file_name

    # Run debyer and read in the pdf
    os.system(cmd)
    fptr_pdf = open("%s.g" % output, 'r').read().split("\n")
    i = 0
    while fptr_pdf[i].strip().startswith("#"):
        i += 1
    j = len(fptr_pdf) - 1
    while fptr_pdf[j].strip() == "":
        j -= 1
    fptr_pdf = fptr_pdf[i:j + 1]

    pdf = [(float(a.split()[0]), float(a.split()[1])) for a in fptr_pdf]

    if not persist:
        os.system("rm %s.g" % output)
        os.system("rm %s.xyz" % file_name)

    return pdf
