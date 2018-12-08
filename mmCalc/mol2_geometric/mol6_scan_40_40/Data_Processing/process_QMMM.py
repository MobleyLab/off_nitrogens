

###WIP:
### This script will take an .xyz (scan-final.xyz) from the finished QM torsion scan and generate/process data to create a plot of QM and MM energies.


#Imports





#Functions
def SDF2oemol(sdffile):
    """
    Description:
    Takes in the sdf file and converts it to an oemol

    Input:
    sdffile: .sdf file

    Returns:
    Oemol: Oemol object of the molecule stored in the input .sdf file
    """
    #WIP

    return oemol


def QM2oemol(xyzfile, oemol):
    """
    Description:
    Takes in an xyz file and creates oemol objects with the coordinates in the original .xyz file
    with QM energies. The .xyz files are formatted as an output geomeTRIC xyz file, scan-final.xyz.    The original energies are stored in Hartree but are converted to kcal/mol in the QM energy tag.
    Input:
    xyzfile:Take in scan-final.xyz file which contains the geometry and energy outputs from a
    geomeTRIC torsion scan
    oemol:The oemol of the molecule involved in the torsion scan of the .xyz file

    Returns:
    oemolList: A list of OEmols that contain the QM data tagged as "QMEng"
    """
    #WIP


    return oemolList

def GetMM(oemol, FF, FFRmNit):
    """
    Description:
    Takes in an oemol and calculates the MM energies with two forcefields, one which has removed
    nitrogen improper parameters. The data from these calculations is stored in the oemol object
    as MM and MMRmNit. The units of energy are KCal/mol.

    Input:
    oemol: A single oemol object
    FF: .offxml file of smirnoff99Frosst.offxml
    FFRmNit: .offxml file of smirnoff99Frosst.offxml with the removed nitrogen improper parameter

    Return:
    oemolMM: An oemol with the MM energies stored in tags MM and MMRmNit in kcal/mol
    """

    #WIP

    return oemolMM

def findAngles(oemol, constraint, Scan=True):
    """
    Description:
    Calculates the improper (and potentially valence) angles of an oemol and stores the value in
    tag "improper" and "valence" for an oemol.

    Input:
    oemol: An oemol object
    constraint: An constraint.txt file used in an input for an geomeTRIC scan
    Scan: A boolean, True = 1d torsion scan, False = 2d torsion scan. If 2d torsion scan, the
    valence paramters will also be calculated.

    Return:
    oemolImp: An oemol with the calculated improper or valence angles with corresponding tags
    "improper" or "valence"
    """
    #WIP
    return oemolImp

def makeOEB(oemolList, title):
    """
    Description:
    Takes in an oemol list and creates an output OEB file.

    Input:
    oemolList: A list of oemols
    title: The title of the OEB file.

    Return:
    """
    #WIP
    return



#Main
