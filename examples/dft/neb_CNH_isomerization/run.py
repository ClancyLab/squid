# Squid imports
from squid import neb
from squid import files
from squid import units


def bfgs(frames=files.read_xyz("CNH_HCN.xyz"), maxiter=1000):
    new_opt_params = {'step_size': 0.1,
                      'step_size_adjustment': 0.5,
                      'max_step': 0.2,
                      'linesearch': 'backtrack',
                      'accelerate': True,
                      'reset_step_size': 5,
                      'start_hess': 1.0,
                      'maxiter': maxiter}
    return neb.NEB("neb_test",
                   frames,
                   "! HF-3c",
                   opt="BFGS",
                   new_opt_params=new_opt_params)


def lbfgs(frames=files.read_xyz("CNH_HCN.xyz"), maxiter=1000):
    new_opt_params = {'step_size': 1.0,
                      'step_size_adjustment': 0.5,
                      'max_step': 0.04,
                      'linesearch': None,
                      'accelerate': True,
                      'reset_step_size': 20,
                      'maxiter': maxiter}
    return neb.NEB("neb_test",
                   frames,
                   "! HF-3c",
                   opt="LBFGS",
                   new_opt_params=new_opt_params,
                   ci_neb=True,
                   no_energy=False)


def sd(frames=files.read_xyz("CNH_HCN.xyz"), maxiter=1000):
    new_opt_params = {'step_size': 0.1,
                      'step_size_adjustment': 0.5,
                      'max_step': 0.2,
                      'linesearch': 'backtrack',
                      'accelerate': True,
                      'maxiter': maxiter}
    return neb.NEB("neb_test",
                   frames,
                   "! HF-3c",
                   opt="SD",
                   new_opt_params=new_opt_params)


def quick_min(frames=files.read_xyz("CNH_HCN.xyz"), maxiter=1000):
    new_opt_params = {'dt': 0.1,
                      'verbose': False,
                      'euler': False,
                      'maxiter': maxiter}
    return neb.NEB("neb_test",
                   frames,
                   "! HF-3c",
                   opt="QM",
                   new_opt_params=new_opt_params)


def fire(frames=files.read_xyz("CNH_HCN.xyz"), maxiter=1000):
    new_opt_params = {'dt': 0.1,
                      'euler': False,
                      'maxiter': maxiter}
    return neb.NEB("neb_test",
                   frames,
                   "! HF-3c",
                   opt="FIRE",
                   new_opt_params=new_opt_params)


def conjugate_gradient(frames=files.read_xyz("CNH_HCN.xyz"), maxiter=1000):
    new_opt_params = {'step_size': 0.1,
                      'step_size_adjustment': 0.5,
                      'max_step': 0.2,
                      'linesearch': 'backtrack',
                      'accelerate': True,
                      'maxiter': maxiter,
                      'method': 'FR',
                      'linesearch': 'backtrack'}
    return neb.NEB("squid_cg",
                   frames,
                   "! HF-3c",
                   opt="CG",
                   new_opt_params=new_opt_params)


def scipy(method="slsqp", frames=files.read_xyz("CNH_HCN.xyz"), maxiter=1000):
    new_opt_params = {'disp': False,
                      'gtol': units.convert_energy("Ha", "eV", 1e-3),
                      'eps': 1.4901161193847656e-08,
                      'return_all': False,
                      'maxiter': maxiter,
                      'norm': 2.0,
                      'ftol': 0.0}
    return neb.NEB("scipy",
                   frames,
                   "! HF-3c",
                   opt="scipy_" + method,
                   new_opt_params=new_opt_params)


test = lbfgs()
test.optimize()
