"""
This file contains the main fragment sampling class, which performs a Monte Carlo simulated annealing procedure to fold a protein.
"""

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring import *
init(extra_options='-mute all  -constant_seed')

import numpy as np
import utils
from Protein import Protein
from FragmentSet import FragmentSet

import os
import random

def outfile(filename, output):
    action = 'a+'
    if not os.path.isfile(filename):
        action = 'w'
    with open(filename, action) as f:
        f.write(output)

def debug(input, current_iteration, stop_iteration):
    if False:
        if current_iteration < stop_iteration:
            print(input)


class MCMCSampler(object):
    def __init__(self, seq_protein, k, T_start, T_end, fragset, N, anneal_rate):
        """
        TO DO: initialize necessary variables
        The score function is given to you (Rosetta centroid score function)
        """
        self.scorefxn = create_score_function('score3')

        self.prev_protein = seq_protein
        self.best_pose = Pose()
        self.best_pose.assign(self.prev_protein.pose)

        self.fragset = fragset
        self.k = k
        self.N = N
        self.anneal_rate = anneal_rate

        self.best_score = self.compute_energy(seq_protein)
        self.old_energy = self.compute_energy(seq_protein)

        self.T = T_start
        self.T_end = T_end

        self.prob_accept = 0
        self.iter = 0

        # For graphing
        self.lst_energy = []

    def compute_energy(self, protein):
        """
        TO DO
        Compute energy of protein.
        Hint: look at utils.py
        --------
        Params:
            - protein (Protein object): protein to score
        Return:
            - energy of conformation (float)
        """
        return self.scorefxn(protein.pose)

    def perturb_fragment(self, sample_index, random_fragment):
        """
        TO DO
        Sample from possible fragments for a position, and replace torsion angles of that fragment in the protein.
        ---------
        Params:
            - TO DO
        Returns:
            - TO DO
        """
        # Create a new copy of the original protein to experiment on
        new_pose = Pose()
        new_pose.assign(self.prev_protein.pose)
        self.next_protein = Protein(pose=new_pose)

        # Print the Previous Protein
        debug('Previous Protein:', self.iter, 5)
        for pos in range(1, self.prev_protein.length+1):
            debug('\tPosition: {}\tAngle: {}'.format(pos, self.prev_protein.get_torsion(pos)), self.iter, 5)

        # For each residue in fragment, replace the corresponding torsion angles in copy of protein
        pos = sample_index
        debug('random fragment: {}\nPosition: {}'.format(random_fragment, pos), self.iter, 5)
        for phi, psi in random_fragment:
            self.next_protein.set_torsion(pos, phi, psi)
            pos += 1
        # Print Protein after angle replacement
        debug('New Protein:', self.iter, 5)
        for pos in range(1, self.prev_protein.length+1):
            debug('\tPosition: {}\tAngle: {}'.format(pos, self.next_protein.get_torsion(pos)), self.iter, 5)

    def metropolis_accept(self): # you may want to add more arguments
        """
        TO DO
        Calculate probability of accepting or rejecting move based on Metropolis criterion.
        --------
        Params:
            - TO DO
        Returns:
            - TO DO
        """
        # Compute change in energy
        delta_e = self.new_energy - self.old_energy
        # Calculate a random number between zero and one
        rand_num = np.random.rand()

        # Passes if energy change is negative
        if delta_e <= 0:
            self.prob_accept = 1
            return True
        else:
            # If energy is positive, calculate probability of accepting
            self.prob_accept = np.exp(-delta_e / self.T)
            if self.prob_accept > rand_num:
                return True
            else:
                return False

    def anneal_temp(self):
        """
        TO DO
        Anneal temperature using exponential annealing schedule. Consider kT to be a single variable (i.e. ignore Boltzmann constant)
        --------
        Params:
            - TO DO
        Returns:
            - TO DO
        """
        self.T = self.anneal_rate * self.T

    def step(self):
        """
        TO DO
        Take a single MCMC step. Each step should do the following:
        1. sample position in chain
            - Note: think about positions you can sample a k-mer fragment from. 
              For example, you cannot sample from position 1 because there is no phi angle
        2. sample fragment at that position and replace torsions in a *copied version* of the protein
        3. measure energy after replacing fragment
        4. accept or reject based on Metropolis criterion
            - if accept: incorporate proposed insertion and anneal temperature
            - if reject: sample new fragment (go to step 3)
        """
        # Sample an eligible position within the original protein
        sample_index = random.randint(1, (int(self.prev_protein.length) - self.k))
        # Candidate fragments with lowest rmsd values
        debug('\nn value: {}, sample index: {}'.format(self.N, sample_index), self.iter, 5)
        candidate_fragments = self.fragset.get_lowRMS_fragments(sample_index, self.N)
        debug('\ncandidate_fragments: {}\n'.format(candidate_fragments), self.iter, 5)

        # Run through all possible options of frag candidates before moving on
        fragment_indices = set()
        while len(fragment_indices) < len(candidate_fragments):
            # From list, choose a random fragment
            fragment_index = random.randint(0, (len(candidate_fragments)-1))
            random_fragment = candidate_fragments[fragment_index]

            # Run through all possible options of frag candidates before moving on
            fragment_indices.add(fragment_index)

            # Sample from possible fragments for a position, and replace torsion angles of that fragment in the protein.
            self.perturb_fragment(sample_index, random_fragment)

            # Test the energy of the changed protein copy
            self.new_energy = self.compute_energy(self.next_protein)
            # Metropolis test to see if we should stick with the changed version
            if self.metropolis_accept():
                debug('Passed Metropolis!!\nProbability: {}'.format(self.prob_accept), self.iter, 5)
                # Anneal temp
                self.anneal_temp()
                # Accept the protein changes
                new_pose = Pose()
                new_pose.assign(self.next_protein.pose)
                self.prev_protein = Protein(pose=new_pose)
                # Update energy
                self.old_energy = self.new_energy
                # Update best pose and score if energy is better than previous
                if self.new_energy < self.best_score:
                    self.best_score = self.new_energy
                    self.best_pose = Pose()
                    self.best_pose.assign(self.next_protein.pose)
                return
            debug('Failed Metropolis...\nProbability: {}'.format(self.prob_accept), self.iter, 5)

    def simulate(self):
        """
        TO DO
        Run full MCMC simulation from start_temp to end_temp. 
        Be sure to save the best (lowest-energy) structure, so you can access it after.
        It is also a good idea to track certain variables during the simulation (temp, energy, and more).
        -------- 
        Params:
            - TO DO
        Returns:
            - TO DO
        """
        outfile('kmer_stats.txt', 'iter: \ttemp: \t\t\t\tenergy:\n')

        # Take as many steps as necessary until we reach the target temp
        while self.T >= self.T_end:
            self.step()
            outfile('kmer_stats.txt', '{} \t\t{} \t{}\n'.format(self.iter, self.T, self.old_energy))
            self.lst_energy.append(self.new_energy)
            self.iter += 1


