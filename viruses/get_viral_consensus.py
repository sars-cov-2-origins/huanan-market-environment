import sys
import glob
import pysam
import numpy as np
import operator
from collections import defaultdict
from Bio import SeqIO

## Fasta file
refs = {}
for r in SeqIO.parse('./viral_dereplicated.fasta', 'fasta'):
	refs[r.id] = str(r.seq)

coverages = []

sequence = {}
i = 0 
for b in refs[sys.argv[1]]:
	sequence[i] = {"A":0,"C":0,"T":0,"G":0}
	i += 1

for bam_file_name in glob.glob('../mappings/*.bam'):
	samfile = pysam.AlignmentFile(bam_file_name, "rb")
	
	
	covered_bases = 0
	sequence_sample = {}
	for pileupcolumn in samfile.pileup(sys.argv[1]): ## Contig name
		bases = {"A":0,"C":0,"T":0,"G":0}
		for pileupread in pileupcolumn.pileups:
			if not pileupread.is_del and not pileupread.is_refskip:
				base = pileupread.alignment.query_sequence[pileupread.query_position]
				if base in bases:
					bases[base] += 1
		if sum(bases.values()) > 0:
			covered_bases += 1

		sequence_sample[pileupcolumn.pos] = bases


	if covered_bases >= int(sys.argv[2]): # add it
		print(bam_file_name)
		i = 0
		for b in refs[sys.argv[1]]:
			if i in sequence_sample:
				for base in bases:
					sequence[i][base] += sequence_sample[i][base]
			i += 1

final = ''
i = 0 
for b in refs[sys.argv[1]]:
	best_base = max(sequence[i].items(), key=operator.itemgetter(1))[0]
	if sequence[i][best_base] >= 1:
		final += best_base 
	else:
		final += 'N'
	i += 1
print(final)
# print(sequence)
