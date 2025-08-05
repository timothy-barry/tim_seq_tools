import re
import gzip
from Levenshtein import hamming

def fq(file):
	if re.search('.gz$', file):
		fastq = gzip.open(file, 'rt')
	else:
		fastq = open(file, 'r')
	with fastq as f:
		while True:
			l1 = f.readline().removesuffix('\n')
			if not l1:
				break
			l2 = f.readline().removesuffix('\n')
			l3 = f.readline().removesuffix('\n')
			l4 = f.readline().removesuffix('\n')
			yield [l1, l2, l3, l4]

# index_fp = "/Users/timbarry/research_offsite/external/crispr-quant/CD34_WT_Cas9_ELANE_e3SA_GSPplus_5uM_S86_S2_L001_I1_001.fastq"
# index_fp = "/Users/timbarry/research_offsite/external/crispr-quant/CD34_WT_Cas9_ELANE_e3SA_GSPplus_5uM_S86_S2_L001_I2_001.fastq"
# index_fp = "/Users/timbarry/research_code/guideseq/test/data/demultiplexed/control.i1.fastq"
# returns the number of reads that differ from the consensus read by a hamming distance of more than 2
def compute_n_different_from_consensus_read(index_fp, hamming_dist = 2):
	zeroith_line = next(fq(index_fp))
	barcode_len = len(zeroith_line[1])
	base_dict = [{"A" : 0, "T" : 0, "G" : 0, "C" : 0, "N" : 0} for i in range(barcode_len)]

	# iterate over lines, incrementing base in each position
	for line in fq(index_fp):
		curr_barcode = line[1]
		for i, base in enumerate(curr_barcode):
			base_dict[i][base] += 1

	# obtain the consensus barcode
	consensus_barcode = ""
	for curr_dict in base_dict:
		curr_base = max(curr_dict, key = curr_dict.get)
		consensus_barcode += curr_base

	# check for match of each read against the consensus barcode (hamming distance = 1)
	n_mismatches = 0
	for line in fq(index_fp):
		curr_barcode = line[1]
		if hamming(curr_barcode, consensus_barcode) >= hamming_dist:
			n_mismatches += 1

	return n_mismatches
