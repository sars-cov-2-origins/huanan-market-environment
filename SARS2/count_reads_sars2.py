import glob
import sys
import pysam
import operator
from collections import defaultdict
from Bio import SeqIO

CONTIG_EDGE = 200
READ_ALN_FRACTION = 0.95
MIN_MAPQ = 20
MIN_READ_ANI = 0.97


for record in SeqIO.parse(sys.argv[1], 'fasta'):
	ref = record.seq
	ref_id = record.id
	break

print("Sample\tRead_count\tCovered_bases")
for fn in glob.glob(sys.argv[2].rstrip("/") + '/*.bam'):
	samfile = pysam.AlignmentFile(fn, "rb")

	read_count = 0
	reads = set()
	covered_positions = set()
	for read in samfile.fetch(record.id, CONTIG_EDGE, len(ref)-CONTIG_EDGE):
		read_len = len(read.get_reference_positions())
		read_pid = 1 - (int(read.get_tag("NM")) / read_len)
		read_aln_len = float(read_len) / float(read.infer_query_length())
		if read_pid >= MIN_READ_ANI:
			if read_aln_len >= READ_ALN_FRACTION and read.mapping_quality >= MIN_MAPQ:
				covered_positions.update(read.get_reference_positions())
				if read.query_name not in reads: # Only count pairs once
					reads.add(read.query_name)
					read_count += 1
	print(fn.split(".bam")[0].split("_")[-2].split("/")[-1] + "\t" + str(read_count) + "\t" + str(len(covered_positions)))
