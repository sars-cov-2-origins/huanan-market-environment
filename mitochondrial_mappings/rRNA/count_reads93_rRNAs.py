import glob
import sys
import pysam
import operator
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
seqs = {}
lengths = {}

for record in SeqIO.parse(sys.argv[1], 'fasta'):
	seqs[record.id] = record.description.replace(" ","_")
	lengths[record.id] = len(record.seq)

rrna_positions = defaultdict(list)
f = open('rrna_positions.tsv')
for line in f.readlines():
	rrna_positions[line.split()[0]].append((int(line.split()[1]), int(line.split()[2].strip())))

results = defaultdict(list)
for fn in glob.glob(sys.argv[2].rstrip("/") + '/*.bam'):

	samfile = pysam.AlignmentFile(fn, "rb")
	fn_name = fn.split(".bam")[0].split("_")[-2].split("/")[-1]
	print(fn_name)
	for seq in rrna_positions:
		read_count = 0
		covered_positions = set()
		reads = set()

		for read in samfile.fetch(seq, 200, lengths[seq]-200):
			# check if it in rRNA
			in_rrna = False
			read_start = min(read.get_reference_positions())
			read_stop = max(read.get_reference_positions())
			for rrna in rrna_positions[seq]:
				range_rrna = set(range(rrna[0],rrna[1]+1))
				range_read = set(range(read_start,read_stop+1))
				if len(range_rrna.intersection(range_read)) > 0:
					in_rrna = True
#			print(rrna_positions[seq])
#			print((read_start,read_stop))
#			print(in_rrna)
			if not in_rrna:
				read_len = len(read.get_reference_positions())
				read_pid = 1 - (int(read.get_tag("NM")) / read_len)
				if read_pid >= 0.95:
					if len(read.get_reference_positions()) >= 40 and read.mapping_quality >= 20:
						covered_positions.update(read.get_reference_positions())
						if read.query_name not in reads: # Only count pairs once
							read_count += 1
							reads.add(read.query_name)

		results['Sample_Name'].append(fn_name)
		results['reference'].append(seq)
		results['read_count'].append(read_count)
		results['covered_bases'].append(len(covered_positions))

results = pd.DataFrame(results)

counts = results.pivot(columns='reference',index='Sample_Name',values='read_count')
print(counts.shape)
counts = counts[counts.columns[counts.sum() > 0]]
print(counts.shape)
counts = counts.reset_index()
counts.to_csv("mt93_counts_rrna.tsv", sep="\t", index=None)

breadth = results.pivot(columns = 'reference', index='Sample_Name', values='covered_bases')
print(breadth.shape)
breadth = breadth[breadth.columns[breadth.sum() > 0]]
print(breadth.shape)
breadth = breadth.reset_index()
breadth.to_csv("mt93_coveredbases_rrna.tsv", sep="\t", index=None)
