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

results = defaultdict(list)
for fn in glob.glob(sys.argv[2].rstrip("/") + '/*.bam'):

	samfile = pysam.AlignmentFile(fn, "rb")
	fn_name = fn.split(".bam")[0].split("_")[-2].split("/")[-1]
	print(fn_name)
	for seq in seqs:
		read_count = 0
		covered_positions = set()
		reads = set()

		for read in samfile.fetch(seq, 200, lengths[seq]-200):
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
counts.to_csv("mt98_counts.tsv", sep="\t", index=None)

breadth = results.pivot(columns = 'reference', index='Sample_Name', values='covered_bases')
print(breadth.shape)
breadth = breadth[breadth.columns[breadth.sum() > 0]]
print(breadth.shape)
breadth = breadth.reset_index()
breadth.to_csv("mt98_coveredbases.tsv", sep="\t", index=None)
