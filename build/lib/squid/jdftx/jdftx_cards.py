class QeControl:
    '''
    The X card.  This is where X

    **Parameters**

        calculation: *str*
            What type of calculation is to be run.  Options are:

            - scf:
            - nscf:
            - bands:
            - relax:
            - md:
            - vc-relax:
            - vc-md:
    '''

    def __init__(self, calculation="scf", title="", verbosity="low",
                 restart_mode="from_scratch", wf_collect=False, nstep=None,
                 iprint=None, tstress=None):
        raise Exception("THE CARD IS NOT DONE!")
