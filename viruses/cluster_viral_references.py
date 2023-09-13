
final = set()
lens = {}
from Bio import SeqIO

descs = {}
seqs = {}
for record in SeqIO.parse('./refseq_and_marketpapers.fna', 'fasta'):
	lens[record.id] = len(record.seq)
	descs[record.id] = record.description
	seqs[record.id] = record.seq

f = open('all_vs_all.mash.tsv')

from collections import defaultdict
v2v = defaultdict(list)
for line in f.readlines():
	query = line.split("\t")[0]
	ref = line.split("\t")[1]
	dist = float(line.split("\t")[2])
	afrac = line.split("\t")[-1].strip()
	afrac = float(afrac.split("/")[0]) / float(afrac.split("/")[1])
	if afrac >= 0.25 and dist <= 0.05:
		v2v[ref].append(query)
		v2v[query].append(ref)

seqs_ordered = sorted(lens, key=lens.get, reverse=True)

for seq in seqs_ordered:
	found = False
	for cluster_other in v2v[seq]:
		if cluster_other in final:
			found = True
	if not found:
		final.add(seq)
f.close()
for f in final:
	print(">" + descs[f])
	print(seqs[f])
