import utils
import numpy as np
import os
import json


def outfile(filename, output):
    action = 'a+'
    if not os.path.isfile(filename):
        action = 'w'
    with open(filename, action) as f:
        f.write(output)

def debug(input):
	if False:
		print(input)

class FragmentSet(object):
	def __init__(self, fragfile, rmsdfile):
		"""
		This class contains the fragment library for the input protein. It must do the following:
		- Read in fragment file and parse fragments at each position. Fragment files are of the form <protein>_<frag_len>mers.frag
		- Read in RMSD file containing pre-calculated RMSD to native structure for each fragment at each position.
		- Based on fragments and their corresponding RMSDs, rank fragments at each position by RMSD
		"""

		# Parse the RMSD into a dictionary of dictionaries
		rmsd_vals = {}

		# Keep track of the length of our k-mer
		self.k = int(fragfile[-10])

		with open(rmsdfile, 'r') as r:
			r_lines = r.readlines()

		for r_line in r_lines:
			r_line = [x.strip() for x in r_line.split()]
			# Set the key to be the position number, value an empty dict
			if not int(r_line[0]) in rmsd_vals:
				rmsd_vals[int(r_line[0])] = {}
			# Each inner dict's key will be frag number, value will be RMSD
			rmsd_vals[int(r_line[0])][int(r_line[1])] = float(r_line[2])

		# Create a position dictionary of fragment dictionaries of list of
		# phi/psy tuples and RMSD value.
		with open(fragfile, 'r') as f:
			f_lines = f.readlines()

		position = 0
		fragment = 0
		master_d = {}
		angles = []

		for f_line in f_lines:
			# Update position variable and add empty dictionary to master dictionary for each position
			if 'position:' in f_line:
				f_stripped = [x.strip() for x in f_line.split()]
				position = int(f_stripped[1])
				if position not in master_d:
					master_d[position] = {}
				continue
			# When empty line, add the coordinate list of tuples and RMSD.
			if f_line == '\n':
				# The end of each file is >1 empty line. Catches this.
				# Checks if the fragment is in the position dictionary.
				if fragment > 0:
					master_d[position][fragment] = [angles, rmsd_vals[position][fragment]]
				angles = []
				continue
			# Update the phi, psi values and fragment number
			else:
				# debug(('else statement', f_line))
				phi = float(f_line[18:27].strip())
				psi = float(f_line[28:36].strip())
				angles.append((phi, psi))
				fragment = int(f_line[91:].strip())


		# Convert inner dictionaries into a list which can be sorted
		convert = lambda y: [[v[1], k, v[0]] for k, v in y.items()]
		# Sort the list by the rmsd value
		for pos in master_d:
			# Put all values of master_d in a list so that it can be sorted
			master_d[pos] = convert(master_d[pos])
			# Sort by rmsd value
			master_d[pos].sort(key=lambda x: x[0])

		# {Position: [rmsd, fragment, [(phi, psi)]]} sorted by rmsd values.
		self.master_d = master_d

		outfile(fragfile + '_fragset.txt', json.dumps(master_d))

	def get_lowRMS_fragments(self, pos, N):
		"""
		Returns the top-ranked fragments by RMSD at a defined position in the chain
		--------
		Params
			- pos (int): fragment position in chain (1-indexed)
			- N (int): number of fragments to return
		Returns
			- lowRMS_fragments (list): top N fragments at pos by RMSD. This should be a list of lists of (phi, psi) tuples. 
			  For example, a 3-mer fragment could be represented as the following: [(-60.892, 142.456), (-72.281, 128.933), (-132.337, -175.477)]
		"""
		return [x[2] for x in self.master_d[pos][:N]]

