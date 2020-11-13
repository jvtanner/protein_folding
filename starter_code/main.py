"""
This is the master file, which you should use to set up and run the simulations.
You may define functions or classes as necessary

For an input sequence, do the following (see project page for details):
	1. load sequence from fasta file and initialize protein into extended configuration
	2. Run a series of simulations (n=5 or 10):
		- Perform MCMC sampling with 9-mer fragments from kT=100 to kT=1 (assembly stage)
		- Perform MCMC sampling with 3-mer fragments from kT=1 to kT=0.1 (refinement stage)
		- Take best (lowest-energy) structure after refinement and perform energy minimization (see utils.relax)
		- Log energy and RMSD to native structure after minimization
	3. Visualize lowest-RMSD structure in PyMol

"""

from pyrosetta import *
from pyrosetta.rosetta import *
init(extra_options='-mute all -constant_seed')

from Protein import Protein
from FragmentSampler import MCMCSampler
from FragmentSet import FragmentSet

import matplotlib.pyplot as plt
import utils
import argparse
import numpy as np
import time
import os

def graphics(total_energy, anneal_rate):
    print('Anneal rate: {}'.format(anneal_rate))

    plt.figure()
    plt.title('Anneal Rate: {}'.format(anneal_rate))
    # plt.xlim([0, 100])
    # plt.ylim([0, 5])
    plt.hist(total_energy)
    plt.show()

def main():

    parser = argparse.ArgumentParser(description='Optimize protein 3D structure')
    parser.add_argument('--fasta', type=str, help='sequence data')
    parser.add_argument('--logdir', type=str, help='the log outputs')
    parser.add_argument('--nsims', type=int, help='number of simulations')
    parser.add_argument('--nfrags', type=int, help='number of fragments to sample from at each iteration')
    parser.add_argument('--anneal_rate', type=float, help='temperature annealing parameter')

    args = parser.parse_args()

    # Read in the sequence in the fasta file
    with open(args.fasta, 'r') as f:
        seq = f.readlines()[1].strip()

    protein_name = args.fasta.split('.')[0]

    # pbd file
    pdb_native = protein_name + '.pdb'

    # create the references to the other files from the protein name
    frag_3mer = protein_name + '_3mers.frag'
    frag_9mer = protein_name + '_9mers.frag'
    rmsd_3mer = protein_name + '_3mers.rmsd'
    rmsd_9mer = protein_name + '_9mers.rmsd'

    # Create fragsets, 3 and 9mers.
    fragset_3mer = FragmentSet(frag_3mer, rmsd_3mer)
    fragset_9mer = FragmentSet(frag_9mer, rmsd_9mer)

    # Create instance of Protein and extend
    seq_protein = Protein(sequence=seq)
    seq_protein.init_extended()

    args.nfrags = 3
    args.anneal_rate = .999
    args.nsims = 1

    total_energy = []
    with open(os.path.join(args.logdir, 'simulation_summary_1fw4.txt'), 'w') as f:
        for i in range(args.nsims):
            # Run simulation for 9-mer.
            MCMC_9 = MCMCSampler(seq_protein, 9, 100, 1, fragset_9mer, args.nfrags, args.anneal_rate)
            MCMC_9.simulate()
            # For graphing
            energy_9 = MCMC_9.lst_energy
            # Run simulation for 3-mer, pass in the best scoring protein from MCMC_9 simulation.
            MCMC_3 = MCMCSampler(Protein(pose=MCMC_9.best_pose), 3, 1, 0.1, fragset_3mer, args.nfrags, args.anneal_rate)
            MCMC_3.simulate()
            # For graphing
            energy_3 = MCMC_3.lst_energy
            total_energy = energy_9 + energy_3
            # Relax the protein
            MCMC_3.best_pose.dump_pdb('1ubq_best_3mer.pdb')
            final_protein, final_rmsd, final_score = utils.relax('1ubq_best_3mer.pdb', pdb_native)
            # Keep track of the simulation, best scores, and best rmsd values
            f.write(str(i) + '\t' + str(final_score) + '\t' + str(final_rmsd) + '\n')

    # graphics(total_energy, args.anneal_rate)

if __name__ == '__main__':
    main()
