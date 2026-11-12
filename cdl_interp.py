"""
cdl_interp.py  - run a phenix program with CDL's interpolated targets

phenix's conformation-dependent library normally takes backbone bond and
angle ideals from the phi/psi bin a residue is in, so the targets jump at bin
edges.  mmtbx.conformation_dependent_library.update_restraints() can instead
interpolate between bins (interpolate=True, the "slow mode"), but phenix 2.2
exposes no parameter for it.  This wrapper turns it on for every call, at
restraint setup and in each macro-cycle's CDL update, then runs the program.

usage: phenix.python cdl_interp.py <module> [program arguments]
  e.g. phenix.python cdl_interp.py phenix.command_line.phenix_refine model.pdb data.mtz ...
       phenix.python cdl_interp.py mmtbx.command_line.geometry_minimization model.pdb ...
"""
import sys
import runpy
import mmtbx.conformation_dependent_library as cdl

_update_restraints = cdl.update_restraints

def update_restraints(*args, **kwds):
    kwds["interpolate"] = True
    return _update_restraints(*args, **kwds)

# callers import update_restraints from the module when they call it, so they get this one
cdl.update_restraints = update_restraints

if len(sys.argv) < 2:
    sys.exit(__doc__)
module = sys.argv[1]
sys.argv = [module] + sys.argv[2:]
print("cdl_interp: CDL interpolation on for %s" % module)
runpy.run_module(module, run_name="__main__", alter_sys=True)
