from squid import g09
from squid import orca
from squid import files


# Run water simulation
def opt_water():
    frames = files.read_xyz('water.xyz')
    return g09.job('water',
                   'HSEH1PBE/cc-pVTZ OPT=() SCRF(Solvent=Toluene)',
                   atoms=frames,
                   queue=None,
                   force=True)


def opt_water_orca():
    frames = files.read_xyz('water.xyz')
    return orca.job('water',
                    '! pw6b95 def2-TZVP D3BJ OPT NumFreq',
                    atoms=frames,
                    queue=None)


# job = opt_water()
# job.wait()
# g09.cubegen_analysis("water", orbital=3)

opt_water_orca()
