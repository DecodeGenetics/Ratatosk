#!/usr/bin/env python

import sys
import os
import random

def replaceIUPAC(seq=None):

	dna_alpha=['A', 'C', 'G', 'T']
	
	return ''.join([c if c in dna_alpha else random.choice(dna_alpha) for c in seq.upper()])

if __name__ == "__main__":

	lines_per_record = 4 # Assume FASTQ input

	try:
		fn = sys.argv[1]
	except IndexError as ie:
		raise SystemError("Error: Specify file name\n")

	if not os.path.exists(fn):
		raise SystemError("Error: File does not exist\n")

	with open(fn, 'r') as fh:

		record = []

		for line in fh:

			record.append(line.rstrip())

			if (len(record) == lines_per_record):

				record[1] = replaceIUPAC(record[1])

				for i in range(0, lines_per_record): print(record[i])

				record = []
