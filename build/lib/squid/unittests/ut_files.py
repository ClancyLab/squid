from squid import files, structures


def test_read_cml():

    # atoms, bonds, angles, dihedrals = files.read_cml("misc/acetone.cml", new_method=False)

    P, molecs = files.read_cml("misc/acetone.cml", new_method=True)
    P.set_opls_mask()
    # print P.dump_style("opls")

    # print P.dihedral_params
    # print P.dihedral_params.index(('4', '3', '3', '13'))

    # print molecs[0].bonds[0].type

    test_sys = structures.System()
    test_sys.add(molecs[0])

    files.write_lammps_data(test_sys, new_method=True, params=P)

    return False


if __name__ == "__main__":
    assert test_read_cml(), "Error - read cml failed!"
