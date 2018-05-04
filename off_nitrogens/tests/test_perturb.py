from off_nitrogens.calc_improper import *
from off_nitrogens.perturb_angle import *
import numpy as np

def test_move_valence():
    """Test a (changed valence, same improper) move on simple coordinates."""
    # atom0 is central, atom3 is being moved
    atom0 = np.asarray([0,0,0])
    atom1 = np.asarray([0,1,-1])
    atom2 = np.asarray([-0.866,-0.5,-1])
    atom3 = np.asarray([0.866,-0.5,-1])
    atom0, atom1, atom2, atom3_rot = perturb_valence(atom0, atom1, atom2, atom3, 60., True)

    # check improper but make sure the central is first and moved is last
    imp_before = calc_improper_angle(atom0, atom1, atom2, atom3)
    imp_after = calc_improper_angle(atom0, atom1, atom2, atom3_rot)
    if abs(imp_after-imp_before) > 0.01:
        raise Exception("After moving valence angle, improper angle is NOT maintained.")
    # check valence angles after move (before move, all are about 75.52 degrees)
    ang_012 = calc_valence_angle(atom0, atom1, atom2)
    ang_013 = calc_valence_angle(atom0, atom1, atom3_rot)
    ang_023 = calc_valence_angle(atom0, atom2, atom3_rot)
    if abs(ang_012-75.522) > 0.01:
        raise Exception("After moving valence angle, the unaffected valence has changed.")
    if abs(ang_013-89.999) > 0.01:
        raise Exception("After moving valence angle by 60 deg., the 0-1-3 final valence is wrong ({}).".format(ang_013))
    if abs(ang_023-41.408) > 0.01:
        raise Exception("After moving valence angle by 60 deg., the 0-2-3 final valence is wrong ({}).".format(ang_023))

