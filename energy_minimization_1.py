from openforcefield.typing.engines.smirnoff import *
from openforcefield.utils import get_data_filename, extractPositionsFromOEMol, generateTopologyFromOEMol
from openeye.oechem import *
import oenotebook as oenb
from openeye.oeomega import * # conformer generation
from openeye.oequacpac import * #for partial charge assignment



def energyminimization(oemol, nsteps, outprefix='molecule'):
    """ Energy minimization calculation on oemol
    Arguments:
    _________
        oemol: (OpenEye OEMol object)
            OpenEye OEMol of the molecule to simulate. Must have all hydrogens and have 3D conformation.
        nsteps: integer
            Number of 2 femtosecond timesteps to take
        outprefix (optional, default 'molecule'): string
            Prefix for output files (trajectory/Topology/etc.).
    Returns:
    ________
        status: (bool)
            True/False depending on whether task succeeded
        system: (OpenMM System)
            system as simulated
        topology: (OpenMM Topology)
            Topology for system
        positions: (simtk.unit position array)
            final position array

    """

    # Prep forcefield, create system and Topology
    #ff = ForceField('remove_nonbonded.ffxml')
    ff = ForceField('remove_almosteverything.ffxml')

    #topology = generateTopologyFromOEMol(oemol)
    from oeommtools.utils import oemol_to_openmmTop
    topology, positions = oemol_to_openmmTop(oemol)
    system = ff.createSystem(topology, [oemol])

    positions = extractPositionsFromOEMol(oemol)

    # Energy minimize
    # Even though we're just going to minimize, we still have to set up an integrator, since a Simulation needs one
    integrator = openmm.VerletIntegrator(2.0*unit.femtoseconds)
    # Prep the Simulation using the parameterized system, the integrator, and the topology
    simulation = app.Simulation(topology, system, integrator)
    # Copy in the positions
    simulation.context.setPositions(positions)

    # Get initial state and energy; print
    state = simulation.context.getState(getEnergy = True, getPositions=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    print(outprefix + " Energy before minimization (kcal/mol): %.2g" % energy)

    # Write out a PDB
    from oeommtools.utils import openmmTop_to_oemol
    outmol = openmmTop_to_oemol( topology, state.getPositions())
    ofile = oemolostream(outprefix+'_initial_newfield.pdb')
    OEWriteMolecule(ofile, outmol)
    ofile.close()



    # Minimize, get final state and energy and print
    simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    print("Energy after minimization (kcal/mol): %.2g" % energy)
    newpositions = state.getPositions()



    # Write out a PDB
    from oeommtools.utils import openmmTop_to_oemol
    outmol = openmmTop_to_oemol( topology, state.getPositions())
    ofile = oemolostream(outprefix+'newfield.pdb')
    OEWriteMolecule(ofile, outmol)
    ofile.close()

    # Return
    return True, system, topology, state.getPositions()





def input_energy_minimization(smiles, name):
    """Reads in SMILE strings of molecules, creates OEmols and runs energy minimization
    """
    mol = OEMol()
    OESmilesToMol(mol, smiles)
    omega = OEOmega()
    omega.SetMaxConfs(100) #Generate up to 100 conformers since we'll use for docking
    omega.SetIncludeInput(False)
    omega.SetStrictStereo(False) #Pick random stereoisomer if stereochemistry not provided

    #Initialize charge generation
    chargeEngine = OEAM1BCCCharges()

    status = omega(mol)

    # Assign atom names; currently there is a limitation where molecules must have atom names
    OETriposAtomTypes(mol)
    OETriposAtomNames(mol)


    if not status:
        print("error generating conformers.")
    #OEAssignCharges(mol, chargeEngine)

    energyminimization(OEMol(mol), 100000, outprefix=name)

#list of my molecules
molecules = {1:'CNC', 2:'CNC(=O)C', 3:'CNC(=O)OC', 4:'CNC(=O)NC', 5:'CNS(=O)(=O)C', 6:'CS(=O)(=O)Nc1ncncc1', 7:'CS(=O)(=O)Nc1ccccc1', 8:'CNc1ccc([O-])cc1', 9:'CNc1ccc(N)cc1', 10:'CNc1ccccc1', 11:'CNc1ncncc1', 12:'CNc1ccncc1'}


#for molecules, use input energy minimization and run energy minization
for k in molecules:
	input_energy_minimization(molecules[k] , name = str(k)+ '_' +  molecules[k] )




#input_energy_minimization('CNc1ccccc1', name = 'aniline')
