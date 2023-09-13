import glob
import sys
import pysam
import operator
import pandas as pd
from collections import defaultdict
from Bio import SeqIO

## Important constants

CONTIG_EDGE = 200 
READ_ALN_FRACTION = 0.95
MIN_MAPQ = 20
MIN_READ_ANI = 0.97

plant_insect_phage_exclude = set()

keep = 'Pepper mild mottle virus' #keep pmmov as fecal indicator

f = open('plant_insect_phage_exclude.txt')
for line in f.readlines():
	plant_insect_phage_exclude.add(line.strip())
f.close()

amplicon_samples = ['SRR23971533', 'SRR23971484', 'SRR23971591', 'SRR23971473', 'SRR23971416'] # don't include amplicon samples

seqs = {}
for record in SeqIO.parse(sys.argv[1], 'fasta'):
	if len(record.seq) >= 1000:
		dont_include = False
		for b in plant_insect_phage_exclude:
			if b in record.description and keep not in record.description:
				dont_include = True
		if not dont_include:
			seqs[record.id] = [len(record.seq), "_".join(record.description.split()[1:])]

seq_masked = defaultdict(set)

f = open('viral_dereplicated_masked.tsv')
for line in f.readlines():
	seq = line.split()[0].split(">")[1]
	start = int(line.split("\t")[-2])
	stop = int(line.split("\t")[-1])
	if seq in seqs:
		seq_masked[seq].update([x for x in range(start,stop+1)])
f.close()

results = []

for fn in glob.glob(sys.argv[2].rstrip("/") + '/*.bam'):
	print(fn)
	samfile = pysam.AlignmentFile(fn, "rb")
	fn_name = fn.split(".bam")[0].split("_")[-2].split("/")[-1]
	if fn_name not in amplicon_samples:
		for seq in seqs:
			read_count = 0
			covered_positions = set()
			reads = set()
			for read in samfile.fetch(seq, CONTIG_EDGE, seqs[seq][0]-CONTIG_EDGE):

				read_len = len(set(read.get_reference_positions()).difference(seq_masked[seq])) # ignore masked sites
				if read_len > 0:
					read_pid = 1 - (int(read.get_tag("NM")) / read_len)
					read_aln_len = float(read_len) / float(read.infer_query_length())
					if read_pid >= MIN_READ_ANI:
						if read_aln_len >= READ_ALN_FRACTION and read.mapping_quality >= MIN_MAPQ:
							covered_positions.update(read.get_reference_positions())
							if read.query_name not in reads: # Only count pairs once
								read_count += 1
								reads.add(read.query_name)
			if read_count > 0:
				result = {"file":fn.split(".bam")[0].split("_")[-2].split("/")[-1],"contig":seq,"name":seqs[seq][1],"viral_length":seqs[seq][0],"read_count":read_count,"covered_bases":len(covered_positions)}
				results.append(result)

results = pd.DataFrame(results)

viruses_with_coverage = results.query("covered_bases >= 500").contig.unique()

results = results[(results.contig.isin(viruses_with_coverage)) | (results.name.str.contains("Influenza"))] # Make an exception here for Influenza because it is fragmented. Without this we just see the PB2 fragment

print(results)

results.to_csv("filtered_viral_counts_97_95_20_200.tsv", sep="\t", index=False)
