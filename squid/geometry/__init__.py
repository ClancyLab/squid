from squid.geometry.misc import run_unit_tests as run_misc
from squid.geometry.packmol import run_unit_tests as run_packmol
from squid.geometry.spatial import run_unit_tests as run_spatial
from squid.geometry.transform import run_unit_tests as run_transform


def run_all_unit_tests():
    run_misc()
    run_packmol()
    run_spatial()
    run_transform()
