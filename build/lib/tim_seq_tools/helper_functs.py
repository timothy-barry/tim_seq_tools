import re
import gzip
import os
from Levenshtein import hamming
from pathlib import Path

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


# index_fp = "/Users/timbarry/research_offsite/external/crispr-quant/293T-SpRY-Cas9-1620-GSPneg-S79_S87_L001_/293T-SpRY-Cas9-1620-GSPneg-S79_S87_L001_I1_001.fastq.gz"
def simplify_read_ids(index_fp):
	# create new file
	input_path = Path(index_fp)
	new_name = input_path.stem + "_new" + input_path.suffix
	new_path = input_path.with_name(new_name)
	with open(new_path, 'w') as out_file:
		for record in fq(index_fp):  # record: list of 4 lines
			# Simplify the read ID (first line)
			header = record[0].split()[0]  # gets '@SEQ_ID'
			seq = record[1]
			plus = record[2]
			qual = record[3]
			# Write each to file, each with a newline
			out_file.write(f"{header}\n{seq}\n{plus}\n{qual}\n")
	os.remove(index_fp)
	os.rename(new_path, input_path)

def make_umi_table(index_fp, out_file_fp = "", umi_start_pos=0, umi_end_pos=8):
	# create new file
	with open(out_file_fp, 'w') as out_file:
		for record in fq(index_fp):
			curr_umi = record[1][umi_start_pos:umi_end_pos]
			curr_id = record[0]
			out_file.write(f"{curr_id}\t{curr_umi}\n")
