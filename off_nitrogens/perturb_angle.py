#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
perturb_angle.py

Perturb geometry around an improper dihedral angle. The options are:
 1. Change valence angle without changing improper
 2. Change improper angle without changing valence angles

For use in generating geometries to parameterize valence and improper angles

By: Victoria Lim and Jessica Maat

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
import math
import sys

from calc_improper import *

#=============================================================================================
# PRIVATE SUBROUTINES
#=============================================================================================

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.

    Parameters
    ----------
    axis : tuple, list, or numpy array
        vector normal to plane in which rotation is performed
    theta : float
        angle of rotation amount in degrees

    Returns
    -------
    numpy array
        three-dimension rotation matrix

    Reference
    ---------
    https://tinyurl.com/yam5hu78
    """
    axis = np.asarray(axis)
    theta = np.radians(theta)

    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def perturb_valence(atom0, atom1, atom2, atom3, theta, verbose=False):
    """
    Rotate atom3 in the plane of atom1, atom2, and atom3 by theta degrees.
    This involves a counterclockwise rotation relative to the normal vector.
    I.e., if normal vector is (0,0,1) the rotation will be counterclockwise
    looking from +z onto the xy plane. But if the normal vector is (0,0,-1)
    the rotation will be clockwise from the same viewpoint.

    Parameters
    ----------
    atom0 : numpy array
        CENTRAL atom coordinates
    atom1 : numpy array
        outer atom coordinates
    atom2 : numpy array
        outer atom coordinates
    atom3 : numpy array
        outer atom coordinates TO BE MOVED
    theta : float
        how many degrees by which to rotate

    Returns
    -------
    atom* : numpy arrays
        the coordinates in the same order as given in parameters
    rot_mat : numpy array
        three-dimension rotation matrix, returned so that the same
        matrix can be applied to any atoms attached to atom3.

    """

    # outer atoms are atom1 atom2 atom3. get normal vector to that plane.
    v1 = atom2-atom1
    v2 = atom2 - atom3
    w2 = np.cross(v1, v2)
    # calculate rotation matrix
    rot_mat = rotation_matrix(w2, theta)
    atom3_rot = np.dot(rot_mat, atom3)
    # print details of geometry
    if verbose:
        print("\n>>> Perturbing valence while maintaining improper angle...")
        # atom being moved is atom3.
        print("\natom0 original coords:\t",atom0)
        print("atom1 original coords:\t",atom1)
        print("atom2 original coords:\t",atom2)
        print("atom3 original coords:\t",atom3)
        print("atom3 {:.2f} deg rotated: {}".format(60.,atom3_rot))

        # check improper but make sure the central is first and moved is last
        print("\nImproper angle, before: ", calc_improper_angle(atom0, atom1, atom2, atom3))
        print("Improper angle, before: ", calc_improper_angle(atom0, atom2, atom3, atom1))
        print("Improper angle, before: ", calc_improper_angle(atom0, atom3, atom1, atom2))
        print("Improper angle, after: ", calc_improper_angle(atom0, atom1, atom2, atom3_rot))
        print("Improper angle, after: ", calc_improper_angle(atom0, atom2, atom3_rot, atom1))
        print("Improper angle, after: ", calc_improper_angle(atom0, atom3_rot, atom1, atom2))

        ## check valences
        print("\nvalence angle 1, before: ", calc_valence_angle(atom0, atom1, atom2))
        print("valence angle 2, before: ", calc_valence_angle(atom0, atom1, atom3))
        print("valence angle 3, before: ", calc_valence_angle(atom0, atom2, atom3))
        print("valence angle 1, after: ", calc_valence_angle(atom0, atom1, atom2))
        print("valence angle 2, after: ", calc_valence_angle(atom0, atom1, atom3_rot))
        print("valence angle 3, after: ", calc_valence_angle(atom0, atom2, atom3_rot))
        print()

    return atom0, atom1, atom2, atom3_rot, rot_mat

def oemol_perturb_valence(mol, central_atom, outer_atom, theta):
    """
    From an OpenEye OEMol, specify the improper angle and the specific atom of
    that improper that should be perturbed. The improper angles are obtained
    from the find_improper_angles function in the calc_improper script.

    Parameters
    ---------
    mol : OpenEye OEMol
        molecule from which to generate perturbed geometry
    central_atom : string
        atom name in the mol which is central to the improper of interest
        Ex., "N1"
    outer_atom : string
        atom name in the mol which is to be rotated
        Ex., "N1"
    theta : float
        how many degrees by which to rotate

    [TODO]

    """
    # todo [1]
    # calc_improper.py: change find_improper_angles to return the name of the central atom as 5th element in crdlist
    #  - under aidx, add something like: aname = atom.GetName()
    #  - then update this line: crdlist.append((crd0, crd1, crd2, crd3, aname))
    #  - update returns section in docstring
    #  - make sure tests are still passing

    # call find_improper_angle function to get coordinates for some specified improper angle
    # todo [2]

    # call perturb_valence on those coordinates
    # todo [3]



    return # placeholder

# todo [4]
# write a new function to set the new coordinates back to the OEMol
# def update_oemol_coordinates(mol, moved_atom, move_matrix)
#   1. use the name moved_atom to get the OEAtom (atom_moved)
#   2. loop over all atoms connected to moved_atom
#   3. apply the move_matrix on those coordinates: something like
#
#       for connected_atom in ___:
#           old_coords = connected_atom.GetCoords()
#           connected_atom.SetCoords(np.dot(move_matrix, old_coords))
#
#      check the syntax though, I'm just writing pseudo-code
#

def perturb_improper(atom0, atom1, atom2, atom3, theta, verbose=False):
    """
    [TODO] [6]

    """
    # todo [5]
    return # placeholder

# todo [7] write a test for your function

