import sys
import glob
import pysam
import numpy as np
import operator
from collections import defaultdict
from Bio import SeqIO

### SRRA to LabCode + Stall
codes = {}
f = open('sra.tsv')
for line in f.readlines():
	codes[line.split("\t")[-1].strip()] = line.split("\t")[0] + "_" + line.split("\t")[1]
f.close()

## Read in organism list
names = {}
f = open('organisms.txt')
for line in f.readlines():
	names[line.split("\t")[1].strip()] = line.split("\t")[0].replace(" ","_")
f.close()

## Fasta file
refs = {}
for r in SeqIO.parse('./reference/mitochondria_dereplicated93_final.fasta', 'fasta'):
	if r.id in names:
		refs[r.id] = str(r.seq)



for bam_file_name in glob.glob('./sorted/*.bam'):
	samfile = pysam.AlignmentFile(bam_file_name, "rb")


	
	for ref_name in names:
		coverages = []
		sequence = defaultdict(lambda: "N")
		for pileupcolumn in samfile.pileup(ref_name): ## Contig name
			bases = {"A":0,"C":0,"T":0,"G":0}
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					base = pileupread.alignment.query_sequence[pileupread.query_position]
					if base in bases:
						bases[base] += 1

			coverages.append(sum(bases.values()))

			best_base = max(bases.items(), key=operator.itemgetter(1))[0]
			tie_found = False
			for base in bases:
				if base != best_base and bases[base] == bases[best_base]:
					tie_found = True # ties count as N

			if bases[best_base] >= 3 and not tie_found: # At least one base - change to 3 for high confidence calls
				sequence[pileupcolumn.reference_pos] = best_base

		i = 0
		final = ''
		for b in refs[ref_name]:
			final += sequence[i]
			i += 1

		completeness = (len(final) - final.count("N")) / len(final)

		if completeness >= 0.50:
			srr_name = bam_file_name.split("/")[-1].split(".sort.bam")[0]
			print(">" + names[ref_name] + "_" + srr_name + "_" + codes[srr_name] + "_cov" + str(round(np.mean(coverages),1)))
			print(final)