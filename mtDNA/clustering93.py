from collections import defaultdict
from Bio import SeqIO
import pandas as pd

d = pd.read_csv('mt98_coveredbases.tsv', sep="\t")
d = d.set_index("Sample_Name").sum().sort_values()
abunds = defaultdict(int)

for genome in d.index:
    abunds[genome] = int(d[genome])

lengths = {}
seqs = {}
for record in SeqIO.parse('reference/mt_derep98.fasta', 'fasta'):
    seqs[record.id] = [record.description, record.seq]
    lengths[record.id] = len(record.seq)

f = open('genbank_and_refseq_mt.mash')

dists = defaultdict(list)
for line in f.readlines():
    dist = float(line.split("\t")[2])
    query = line.split("\t")[0]
    ref = line.split("\t")[1]
    if query in seqs and ref in seqs:
        if dist <= 0.07: # 93% clustering
            dists[query].append(ref)
            dists[ref].append(query)
f.close()

## Chaining clusterings
for s1 in dists:
    for s2 in dists[s1]:
        for s3 in dists[s2]:
            if s3 not in dists[s1]:
                dists[s1].append(s3)

## get all distances
derep = set()

already_clustered = set()

for cluster in seqs:
    found = False
    for cluster2 in dists[cluster]:
        if cluster2 in already_clustered:
            found = True
            break

    # Take largest
    if not found:
        choices = dists[cluster]
        choices.append(cluster)

        highest_abundance = -1
        most_abundant_centroid = ''
        for c in choices:
            already_clustered.add(c)
            if abunds[c] > highest_abundance:
                most_abundant_centroid = c
                highest_abundance = abunds[c]
        derep.add(most_abundant_centroid)

for sequence in derep:
    print(">" + seqs[sequence][0])
    print(seqs[sequence][1])