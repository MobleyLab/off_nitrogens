import os
import sys

# Import stuff
from openeye.oechem import *
#import oenotebook as oenb
from openeye.oeomega import * # conformer generation
from openeye.oequacpac import * #for partial charge assignment
from openforcefield.typing.engines.smirnoff import *
from openforcefield.utils import get_data_filename, extractPositionsFromOEMol, generateTopologyFromOEMol
from mdtraj.reporters import DCDReporter

def runMD(mol_smiles, ffxml, outdcd):

    # Create empty OEMol
    mol = OEMol()
    # Convert SMILES
    OESmilesToMol(mol, mol_smiles)
    # Draw
    #oenb.draw_mol(mol)

    #initialize omega for conformer generation
    omega = OEOmega()
    omega.SetMaxConfs(100) #Generate up to 100 conformers since we'll use for docking
    omega.SetIncludeInput(False)
    omega.SetStrictStereo(True) #Refuse to generate conformers if stereochemistry not provided

    #Initialize charge generation
    chargeEngine = OEAM1BCCCharges()

    # Set to use a simple neutral pH model
    OESetNeutralpHModel(mol)

    # Generate conformers with Omega; keep only best conformer
    status = omega(mol)
    if not status:
        print("Error generating conformers for %s." % (guest_smiles))

    # Assign AM1-BCC charges
    OEAssignCharges(mol, chargeEngine)

    # Write out PDB of molecule
    ofile = oemolostream('mymolecule.pdb')
    OEWriteMolecule(ofile, mol)
    ofile.close()


    ff = ForceField(ffxml)
    topology = generateTopologyFromOEMol(mol)
    system = ff.createSystem(topology, [mol])

    positions = extractPositionsFromOEMol(mol)

    # Even though we're just going to minimize, we still have to set up an integrator, since a Simulation needs one
    integrator = openmm.VerletIntegrator(2.0*unit.femtoseconds)
    # Prep the Simulation using the parameterized system, the integrator, and the topology
    simulation = app.Simulation(topology, system, integrator)
    # Copy in the positions
    simulation.context.setPositions( positions)

    # Get initial state and energy; print
    state = simulation.context.getState(getEnergy = True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    print("Energy before minimization (kcal/mol): %.2g" % energy)

    # Minimize, get final state and energy and print
    simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    print("Energy after minimization (kcal/mol): %.2g" % energy)
    newpositions = state.getPositions()

    # Set up DCD reporter for storing trajectory; prep for Langevin dynamics
    integrator = openmm.LangevinIntegrator(300*unit.kelvin, 1./unit.picosecond, 2.*unit.femtoseconds)

    # Prep Simulation
    simulation = app.Simulation(topology, system, integrator)
    # Copy in minimized positions
    simulation.context.setPositions(newpositions)

    # Initialize velocities to correct temperature
    simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
    # Set up to write trajectory file to DCD file in data directory every 100 frames
    dcd_reporter = DCDReporter(os.path.join('.', outdcd), 100) #Store every 100 frames
    # Initialize reporters, including a CSV file to store certain stats every 100 frames
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(app.StateDataReporter(os.path.join('.', 'data.csv'), 100, step=True, potentialEnergy=True, temperature=True, density=True))

    # Run the simulation and print start info; store timing
    print("Starting simulation")
    start = time.clock()
    simulation.step(1000) #1000 steps of dynamics
    end = time.clock()

    # Print elapsed time info, finalize trajectory file
    print("Elapsed time %.2f seconds" % (end-start))
    dcd_reporter.close()
    print("Done!")




if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--smiles", required=True,
                        help="SMILES string of the molecule, in quotation marks.")
    parser.add_argument("-f", "--ffxml", required=True,
                        help="Name of the force field to do the MD simulation.")
    parser.add_argument("-o", "--outdcd", required=True,
                        help="Name of the force field to do the MD simulation.")
    args = parser.parse_args()
    opt = vars(args)

    # make sure input file is specified AND exists
    if not opt['ffxml'] or not os.path.exists(opt['ffxml']):
        sys.exit("\nERROR: Specify a valid force field file.\n")

    #insert output file

    runMD(opt['smiles'], opt['ffxml'], opt['outdcd'])
