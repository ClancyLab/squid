from squid.ff_params import Parameters


def test_param_obj():
    '''
    When using the Parameters object, you first ensure you read it in with the
    necessary restrictions.  Afterwards, you can sub-restrict as needed.
    '''
    RESTRICT_LIST = ["13", "24", "46", "48", "1000", "1001"]

    P = Parameters(
        fptr=[
            ("SMRFF", "misc/junk.smrff"),
            ("OPLS", "misc/oplsaa.prm")
        ],
        restrict=RESTRICT_LIST
    )

    ###############################################################################
    print("Test 1 - We print out all parameters in a readable type format.\n\n")

    # This is what we read in
    P.lj_mask = True
    P.coul_mask = True
    P.morse_mask = True
    # P.tersoff_mask = True
    # P.bond_mask = True
    # P.angle_mask = True
    # P.dihedral_mask = True
    print str(P) + "\n\n------------------------\n\n"

    ###############################################################################
    print("Test 2 - We restrict to only 1000 and 1001 and re-print.\n\n")

    # Now, let's restrict for a system that's only 1000 and 1001
    P.restrict = ["1000", "1001"]
    print str(P) + "\n\n------------------------\n\n"

    ###############################################################################
    print("Test 3 - We unrestrict and print in a lammps input type format.\n\n")

    # Now, let's unrestrict and show all the parameters
    P.restrict = RESTRICT_LIST
    print P.dump_style(style="lj/cut/coul/cut")
    print P.dump_style(style="morse")
    print P.dump_style(style="tersoff")
    print P.dump_style(style="opls")
    print "\n\n------------------------\n\n"

    ###############################################################################
    print("Test 4 - We re-restrict and print in a lammps input type format.\n\n")

    # Finally, let's re-restrict and only show some
    P.restrict = ["1000", "1001"]
    print P.dump_style(style="lj/cut/coul/cut")
    print P.dump_style(style="morse")
    print P.dump_style(style="tersoff")
    print P.dump_style(style="opls")
    print "\n\n------------------------\n\n"

    return True


if __name__ == '__main__':
    assert test_param_obj(), "Error - ff_params failed!"
