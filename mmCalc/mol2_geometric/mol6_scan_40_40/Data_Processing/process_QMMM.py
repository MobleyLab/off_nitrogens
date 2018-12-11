###WIP:
### This script will take an .xyz (scan-final.xyz) from the finished QM torsion scan and generate/process data to create a plot of QM and MM energies.


#Imports
from openforcefield.typing.engines.smirnoff import *
from openforcefield.utils import get_data_filename, extractPositionsFromOEMol, generateTopologyFromOEMol
from openeye.oechem import *
#import oenotebook as oenb
from openeye.oeomega import * # conformer generation
from openeye.oequacpac import * #for partial charge assignment
from calc_improper import *




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
    ifs = oechem.oemolistream()
    ifs.open(sdffile)
    oemol = ifs.GetOEGraphMols()
    print(oemol)
    coords = oemol.GetCoords()
    oemol.SetCoords(coords)
    return oemol

def smiles2oemol(smiles):
    """
    Description:
    Takes in an SMILES string and returns an oemol.

    Input:
    smiles: SMILES string for molecule

    Returns:
    mol1: OEMol object
    """
    mol1 = OEMol()
    OESmilesToMol(mol1, smiles)

    #assign charges, necessary to keep?
    chargeEngine = OEAM1BCCCharges()
    OEAssignCharges(mol1, chargeEngine)

    return mol1


def QM2Oemol(xyzfile, oemol):
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

    #open input xyz file
    ifs = oechem.oemolistream()
    ifs.open(xyzfile)

    #generate empty list of oemols
    oemolList=list()

    #iterate through oemols and set coordinates to original oemol
    for mol in ifs.GetOEGraphMols():
        #get coordinates
        print(mol.GetTitle())
        coords = mol.GetCoords()
        print(coords)
        #set coordinates of original oemol with correct bond connectivity
        newmol = copy.deepcopy(oemol)
        newmol.SetCoords(coords)
        #get the energies from the .xyz file title and save in oemol data
        title = mol.GetTitle()
        energy = float(str.split(title)[-1])
        energy_kcal= energy * 627.509
        #set the QM energy in a tag called qm
        newmol.SetData("QM", energy_kcal)

        #append the list with the new mol that has stored QM energies
        oemolList.append(newmol)

        #test energy Data tags
        for k in oemolList:
            print(k.GetData("QM"))
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
    oemol: An oemol with the MM energies stored in tags MM and MMRmNit in kcal/mol
    """

    #prep both force fields, create system and topology, get MM energies and store in data tag
    ff = ForceField(FF)
    topology = generateTopologyFromOEMol(oemol)
    system = ff.createSystem(topology, [oemol])
    positions = extractPositionsFromOEMol(oemol)
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    state = simulation.context.getState(getEnergy = True, getPositions=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    oemol.SetData("MM", energy)

    #repeat with removed nitrogen
    ffn = ForceField(FFRmNit)
    topology = generateTopologyFromOEMol(oemol)
    system = ffn.createSystem(topology, [oemol])
    positions = extractPositionsFromOEMol(oemol)
    simulation = app.Simulation(topology, system, integrator)
    simulation.context.setPositions(positions)
    state = simulation.context.getState(getEnergy = True, getPositions=True)
    energy_rm = state.getPotentialEnergy() / unit.kilocalories_per_mole
    oemol.SetData("MMRmNit", energy_rm)

    return oemol

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
    oemol: An oemol with the calculated improper or valence angles with corresponding tags
    "improper" or "valence"
    """

    #determine the constraints
    #open the constraint file
    constraintfile = open(constraint, "r")



    imp = find_improper_angles(k)[0][0]
    imp_ang = calc_improper_angle(imp[0], imp[2], imp[1], imp[3], True)

    #WIP
    return oemol

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
def processData(sdffile, xyzfile, constraintFile, moltitle, FF, FFRmNit):
    """
    Description:
    Takes in innitial .sdf file, xyzfile from geomeTRIC output, constraint file and a title to
    create an .oeb file with oemols with QM, and MM data.
    """
# fix, out of order function inputs    oemolList = QM2Oemol(SDF2oemol(sdffile), xyzfile)

    for mol in oemolList:
        GetMM(mol, FF, FFRmNit)
        findAngles(oemol, constraint, Scan=True)

    makeOEB(oemolList, title)

#plotting





#test
#QM2Oemol('scan-final.xyz', SDF2oemol('molClass_pyrnit_molecule_6.mol2'))
QM2Oemol('scan-final.xyz', smiles2oemol('CS(=O)(=O)Nc1ncncc1'))
