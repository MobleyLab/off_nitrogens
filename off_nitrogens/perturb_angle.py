#!/usr/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
perturb_angle.py

Perturb geometry around an improper dihedral angle. The options are:
 1. Change valence angle without changing improper
 2. Change improper angle without changing valence angles

For use in generating geometries to parameterize valence and improper angles for an oemol

By: Victoria Lim and Jessica Maat

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import numpy as np
import math
import sys
from openeye import oeomega
from openeye import oechem
from calc_improper import *
#from off_nitrogens.calc_improper import *


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
    v1 = np.asarray(atom2)-np.asarray(atom1)
    v2 = np.asarray(atom2) - np.asarray(atom3)
    w2 = np.cross(v1, v2)
    # calculate rotation matrix
    rot_mat = rotation_matrix(w2, theta)
    length = np.linalg.norm(atom0-atom3)
    atom3_rot = np.dot(rot_mat, atom3)
    new_length = np.linalg.norm(atom0-atom3_rot)

    return atom0, atom1, atom2, atom3_rot






def perturb_improper(atom0, atom1, atom2, atom3, theta, verbose=False):
    """
    Rotate atom3 out of the plane of atom1, atom2  by theta degrees while keeping valence angle the same.

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
        how many degrees by which to move atom upwards

    Returns
    -------
    atom* : numpy arrays
        the coordinates in the same order as given in parameters
    rot_mat : numpy array
        three-dimension rotation matrix, returned so that the same
        matrix can be applied to any atoms attached to atom3.

    """

    # get the vector between atom2 and atom1 and then rotate atom3 about that plane by angle theta.
    v1 = atom2-atom1
    # calculate rotation matrix
    rot_mat = rotation_matrix(v1, theta)
    atom3_rot = np.dot(rot_mat, atom3)

    return atom0, atom1, atom2, atom3_rot





def oemol_perturb(mol, central_atom, outer_atom, angle_type,  theta):
    """
    From an OpenEye OEMol, specify the improper angle and the specific atom of
    that improper that should be perturbed. The improper angles are obtained
    from the find_improper_angles function in the calc_improper script.

    Parameters
    ---------
    mol : OpenEye OEMol
        molecule from which to generate perturbed geometry
    central_atom : indices
        atom indice in the mol which is central to the improper of interest
        Ex., "3"
    outer_atom : indices
        atom indice in the mol which is to be rotated
        Ex., "7"
    True: Boolean
        True = Improper perturbation
        False = Valence perturbation
        Ex., True
    theta : float
        how many degrees by which to rotate

    [TODO]

    """
    #Set perturbation type
    if angle_type == True:
        angle_type = perturb_improper
    else:
        angle_type = perturb_valence


    #determine which improper angle on the atom is of interest and store the coordinates of the improper in  various variables
    cmol = oechem.OEMol(mol)
    center_atom = cmol.GetAtom(oechem.OEHasAtomIdx(central_atom))
    move_atom = cmol.GetAtom(oechem.OEHasAtomIdx(outer_atom))
    center_coord = np.array(cmol.GetCoords(center_atom))
    move_coord = np.array(cmol.GetCoords(move_atom))

    other_coords = list()
    for neighbor in center_atom.GetAtoms():
        if neighbor.GetIdx() !=  outer_atom:
            new_coord = np.array(cmol.GetCoords(neighbor))
            other_coords.append(new_coord)

    atom0, atom1, atom2, atom3_rot = angle_type(center_coord, other_coords[0], other_coords[1], move_coord, theta)
    cmol.SetCoords(move_atom, oechem.OEFloatArray(atom3_rot))

    # Create an output file for visualization: .pdb file
    angle = str(theta)
    perturbation_type = angle_type.__name__

    #adding the index of the trivalent center to the SD tags
    lisltofIdx = "%d, %d, %d, %d" % (indexList[0], indexList[1], indexList[2], indexList[3])
    oechem.OEAddSDData(oemol, "Index list around trivalent nitrogen center: ", listofIdx)


    ofile = oechem.oemolostream(angle + '_' + perturbation_type + 'molecule.mol2')
    oechem.OEWriteMolecule(ofile, cmol)

    ofile.close()

    return cmol

