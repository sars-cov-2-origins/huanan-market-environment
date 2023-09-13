from collections import defaultdict
from Bio import SeqIO

lengths = {}
seqs = {}
for record in SeqIO.parse('genbank_and_refseq_mt.fasta', 'fasta'):
	seqs[record.id] = [record.description, record.seq]
	lengths[record.id] = len(record.seq)


f = open('genbank_and_refseq_mt.mash')

dists = defaultdict(list)
for line in f.readlines():
    dist = float(line.split("\t")[2])
    query = line.split("\t")[0]
    ref = line.split("\t")[1]
    if dist <= 0.02: # 98% clustering
        dists[query].append(ref)
        dists[ref].append(query)
f.close()

## get all distances
derep = set()

high_priority = {'NC_012920.1': 'Homo sapiens mitochondrion, complete genome',
'NC_002008.4': 'Canis lupus familiaris mitochondrion, complete genome', 
'NC_006853.1': 'Bos taurus mitochondrion, complete genome'
}

already_clustered = set()

for key in high_priority.keys():
    derep.add(key)
    already_clustered.add(key)


for cluster in dists:
    found = False
    for cluster2 in dists[cluster]:
        if cluster2 in already_clustered:
            found = True
            break

    # Take largest
    if not found:
    	choices = dists[cluster]
    	choices.append(cluster)

    	refseq = list()
    	for cluster2 in choices:
    		if "NC_" in cluster2:
    			refseq.append(cluster2)

    	if len(refseq) > 0:
    		choices = refseq

    	longest_length = 0
    	longest = ''
    	for c in choices:
    		already_clustered.add(c)
    		if lengths[c] > longest_length:
    			longest = c
    			longest_length = lengths[c]
    	derep.add(longest)

for sequence in derep:
   	print(">" + seqs[sequence][0])
   	print(seqs[sequence][1])