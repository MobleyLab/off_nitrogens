from openforcefield.typing.engines.smirnoff import *
from openforcefield.utils import get_data_filename, extractPositionsFromOEMol, generateTopologyFromOEMol
from openeye.oechem import *
import oenotebook as oenb
from openeye.oeomega import * # conformer generation
from openeye.oequacpac import * #for partial charge assignment
from oeommtools.utils import openmmTop_to_oemol
import argparse
import math
import os, sys




#Create output files
#output file with energies and improper smirks labels
energyf = open('energy.txt','w').close()
smirkf =  open('improper_labels.txt','w').close()
anglef = open('bondsums.txt', 'w').close()
energyf = open('energy.txt','a+')
smirkf =  open('improper_labels.txt','a+')
anglef = open('bondsums.txt', 'a+')


def energyminimization(smiles, nsteps, forcefield, geometry):
    """ Energy minimization calculation on oemol
    Arguments:
    _________
        smiles: (Smiles syntax molecule, can be a list of molecules)
            This will be converted to oemols.
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
    f = open(forcefield + '_' + geometry +' improper_labels.txt','a+')
    for k in smiles:
        outprefix = forcefield +'_' + geometry + '_' + str(k) + '_' +  smiles[k]
        oemol = OEMol()
        OESmilesToMol(oemol, smiles[k])
        omega = OEOmega()
        omega.SetMaxConfs(100) #Generate up to 100 conformers since we'll use for docking
        omega.SetIncludeInput(False)
        omega.SetStrictStereo(False) #Pick random stereoisomer if stereochemistry not provided
         
        
        #Initialize charge generation
        #chargeEngine = OEAM1BCCCharges()
        status = omega(oemol)
        if not status:
            print("error generating conformers.")
        #OEAssignCharges(oemol, chargeEngine)
        
        # Prep forcefield, create system and Topology
        ff = ForceField('forcefield/'+forcefield)    
        topology = generateTopologyFromOEMol(oemol)
        system = ff.createSystem(topology, [oemol])
        positions = extractPositionsFromOEMol(oemol)

        #find smirks labels and write out to file
        labels= ff.labelMolecules( [oemol], verbose = False )
        smirkf.write(outprefix)
        for mol_entry in range(len(labels)):
            for force in labels[mol_entry].keys():
                #print("\n%s:" % force)
                #f.write("\n%s:" % force)
                for (atom_indices, pid, smirks) in labels[mol_entry][force]:
                    atomstr=''
                    for idx in atom_indices:
                        atomstr += '%6s' % idx
                    if pid in ['i1', 'i2', 'i3', 'i4']:
                        print("%s : %s \t smirks %s" % (atomstr, pid, smirks))
                        smirkf.write("\n %s : %s \t smirks %s" % (atomstr, pid, smirks))
    
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
        print(outprefix + " Energy before minimization (kcal/mol) with %s : %.2g" %(forcefield, energy))

        # Write out a PDB
        from oeommtools.utils import openmmTop_to_oemol
        outmol = openmmTop_to_oemol(topology, state.getPositions())
        ofile = oemolostream(outprefix+'_initial.pdb')
        OEWriteMolecule(ofile, outmol)
        ofile.close()

        # Minimize, get final state and energy and print
        simulation.minimizeEnergy()
        state = simulation.context.getState(getEnergy=True, getPositions=True)
        energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
        print(outprefix + "Energy after minimization (kcal/mol): %.2g" % energy)
        newpositions = state.getPositions()

        # Write out a PDB
        from oeommtools.utils import openmmTop_to_oemol
        outmol = openmmTop_to_oemol( topology, state.getPositions())
        ofile = oemolostream(outprefix+'newfield.pdb')
        OEWriteMolecule(ofile, outmol)
        ofile.close()

        #find the angle sum of nitrogen centers and write out to file
        angList = []
        labelList = []
        for atom in outmol.GetAtoms(oechem.OEIsInvertibleNitrogen()):
            aidx = atom.GetIdx()
            nbors = list(atom.GetAtoms())
            ang1 = math.degrees(oechem.OEGetAngle(oemol, nbors[0],atom,nbors[1]))
            ang2 = math.degrees(oechem.OEGetAngle(oemol, nbors[1],atom,nbors[2]))
            ang3 = math.degrees(oechem.OEGetAngle(oemol, nbors[2],atom,nbors[0]))
            ang_sum = math.fsum([ang1,ang2,ang3])
            anglef.write("\n\n%s: sum of angles for N, index %d: %f" % (outprefix, aidx, ang_sum))
            angList.append(ang_sum)
            print(ang_sum)
            labelList.append("{}_{}_{}".format(outprefix,k,aidx))

            





#list of my molecules
molecules = {1:'CNC', 2:'CNC(=O)C', 3:'CNC(=O)OC', 4:'CNC(=O)NC', 5:'CNS(=O)(=O)C', 6:'CS(=O)(=O)Nc1ncncc1', 7:'CS(=O)(=O)Nc1ccccc1', 8:'CNc1ccc([O-])cc1', 9:'CNc1ccc(N)cc1', 10:'CNc1ccccc1', 11:'CNc1ncncc1', 12:'CNc1ccncc1'}

part = {4:'CNC(=O)NC', 6:'CS(=O)(=O)Nc1ncncc1', 9:'CNc1ccc(N)cc1', 10:'CNc1ccccc1', 12:'CNc1ccncc1'}
pla = {2:'CNC(=O)C', 3:'CNC(=O)OC', 11:'CNc1ncncc1'}
pyr = {1:'CNC', 5:'CNS(=O)(=O)C', 7:'CS(=O)(=O)Nc1ccccc1', 8:'CNc1ccc([O-])cc1'}


#original field 
energyminimization(part, 100000, 'jessica.ffxml', 'part')
energyminimization(pla, 100000, 'jessica.ffxml', 'pla')
energyminimization(pyr, 100000, 'jessica.ffxml', 'pyr')



