"""
Example converting from cti to chemkin.
"""

import cantera as ct
import cantera2chemkin as ct2ck

in_mech   = 'example_mech.cti'
out_mech  = 'example_mech.inp'
out_therm = 'example_therm.dat'
out_trans = 'example_trans.dat'

sol = ct.Solution(in_mech)
ct2ck.convertMech(sol, mechOut=out_mech, thermOut=out_therm, tranOut=out_trans)

