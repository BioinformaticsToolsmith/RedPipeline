#
# Author: Hani Zakaria Girgis, PhD
# The Bioinformatics Toolsmith Laboratory
# Purpose: This program runs a pipeline of tools for
# 	detecting repeats in a genome

import sys
import os
import subprocess
import multiprocessing
import pybedtools
from Bio import SeqIO
from Bio import Align
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from clustering import ClusterList
import numpy as np
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
import matplotlib.pyplot as plt
import glob

fa_dir  = sys.argv[1]
out_dir = sys.argv[2]
in_dir = out_dir + '/' + 'Fa'

# Parameters
# An LTR-RT must be coverage by Red's regions
ltr_coverage       = 0.7
tr_coverage        = 0.7
poly_a_min_score   = 6
poly_a_win         = 15
min_len            = 50 # Sequences less than this are excluded from this analysis
short_max_len      = 1000 # Maximum length of unclassified elements
tiny_max_len       = 100 # Macimum length of tiny unclassified elements
sine_max_len       = 1000 # Maximum length of a SINE
dna_min_score      = 10 # Minimum alignment score of inverted repeats of DNA transposons
dna_win            = 20 # The size of the left most window
mite_max_len       = 1000 # Maximum length of a MITE
meshclust_id       = 0.7 # Identity threshold defining sequences in the same cluster
meshclust_min_size = 5 # Minimum number of sequences in a cluster
fastcar_id         = 0.6 # Identity threshold defining similar sequences
helitron_win       = 50 # Search window for palindrom 
helitron_min_score = 16 # Size of palindrom

# Directories where programs outputs are stored
red_dir = out_dir + '/red' 
ltr_dir = out_dir + '/ltr'
tr_dir  = out_dir + '/tr' 

# Important output files
# This file contains clusters with their members (coordinates only)
all_out = out_dir + '/ALL.cluster'
# This file contains the sequences of all centers
# Should be deleted
all_fasta = out_dir + '/ALL.fa'
# This file contains id scores of all centers versus each other
# Should be deleted
all_vs_all = out_dir + '/ALL_vs_ALL.txt'

# Final library output
library = out_dir + '/library.fa'
# A file including information on matching LTR to their interiors
ltr_to_interior = out_dir + '/ltr_to_interior.txt'

# Dictionary: name -> file
# I want to know in which file a sequence can be found
# This is important when a file contains multiple sequences
header_file_dict = {}
ltr_id_dict = {}
ltr_interior_dict = {}

ltr_group          = 'LTR'
interior_group     = 'LTR_INT' 
dna_group          = 'DNA'
interspersed_group = 'INE'
helitron_group     = 'HAIRPIN'
unknown_group      = 'UNKNOWN'

# Credit: https://stackoverflow.com/questions/4417546/constantly-print-subprocess-output-while-process-is-running
def execute(cmd):
    popen = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        print(stdout_line, end='') 
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)

# Step 1
def run_Red():
	print('Making directory:', red_dir)
	if not os.path.isdir(red_dir):
		os.mkdir(red_dir, 0o700)
	print('Running Red')
	cmd = 'Red -gnm %s -frm 2 -rpt %s' %(in_dir, red_dir)
	execute(cmd)

# Step 2
def run_LtrDetector():
	print('Making directory:', ltr_dir)
	if not os.path.isdir(ltr_dir):
		os.mkdir(ltr_dir, 0o700)
	print('Running LtrDetector')
	cmd = 'LtrDetector -fasta %s -destDir %s -nThreads %d -bedFormat' %(in_dir,ltr_dir,multiprocessing.cpu_count())
	execute(cmd)

def format_input():
	"""
	Make sure the header of a fasta sequence does not end in space or includes spaces
	If it includes spaces, they will be replaced with underscores
	"""
	if not os.path.isdir(in_dir):
		os.mkdir(in_dir, 0o700)

	for file_name in os.listdir(fa_dir):
		if '.fa' in file_name:
			copy           = in_dir + '/' +  file_name
			copy_file      = open(copy, 'w')

			original       = fa_dir + '/' + file_name
			original_file  = open(original, 'r')

			for line in original_file:
				if line[0] == '>':
					new_header = '_'.join(line.split()) + '\n'
					copy_file.write(new_header)
				else:
					copy_file.write(line)

			copy_file.close()
			original_file.close()

# Step 3
def run_Look4TRs():
	print('Making directory:', tr_dir)	
	if not os.path.isdir(tr_dir):
		os.mkdir(tr_dir, 0o700)
	print('Running Look4TRs')
	cmd = 'Look4TRs --adr %s --out %s --default' %(in_dir, tr_dir)
	execute(cmd)

# Step 4
#	a. Filter FP by removing any LTR-RT that is not covered (<70%) by Red's repeats
def group_based_on_structure():
	ltr_counter = 0
	for file in os.listdir(red_dir):
		print('Processing file: ', file)
		a = pybedtools.BedTool(ltr_dir + '/' + file[0:file.find('.bed')] + 'Detector.bed')
	
		#ltr_file = file.find('.bed')
		b = pybedtools.BedTool(red_dir + '/' + file )

		c = a.intersect(b, wo=True, sorted=True)

		nick_name = os.path.splitext(file)[0];

		file_tiny  = open(out_dir + '/' + nick_name + "tiny.fa",  "w")
		file_short = open(out_dir + '/' + nick_name + "short.fa", "w")
		file_long  = open(out_dir + '/' + nick_name + "long.fa",  "w")

		# Calculate of Red's repeats within LTR-RT's
		covered_list = [];
		if len(c) > 0:
			p        = c[0]
			p_char   = p[0]
			p_start  = int(p[1])
			p_end    = int(p[2])
			coverage = 0.0

			for n in c[1:]:
				n_char  = n[0]
				n_start = int(n[1])
				n_end   = int(n[2])
				if(p_char == n_char and p_start == n_start and p_end == n_end):
					coverage += int(n[-1])
				else:
					if(coverage/(p_end - p_start) > ltr_coverage):
						covered_list.append(str(n))
					coverage = int(n[-1])
				p_char  = n_char 
				p_start = n_start 
				p_end   = n_end	
			# Handle last case
			if(coverage/(p_end - p_start) > ltr_coverage):
				covered_list.append(str(n))

		# Filtered LTR-RT
		d = pybedtools.BedTool(''.join(covered_list), from_string=True)

		# Non-LTR repeats
		e = b.intersect(d, sorted=True, v=True)
	
		# TR 
		f = pybedtools.BedTool(tr_dir + '/' + file)
		
		# No LTR and No Tr
		g = e.subtract(f, A=True, f=tr_coverage)
		print(len(b), len(e), len(g))
	
		#
		# Open FASTA file
		#
		seq_dict = SeqIO.index(in_dir + '/' + nick_name + '.fa', "fasta")

		# Print the LTRs
		if len(d) > 0:
			ltr_chr = d[0][0]
			ltr_seq = seq_dict[ltr_chr].seq
			for k in d[0:]:
				if ltr_chr != k[0]:
					ltr_chr = k[0]
					ltr_seq = seq_dict[ltr_chr]

				# Process left LTR
				elm_len = int(k[4]) - int(k[3])
				if elm_len < min_len:
					print('Left LTR is too small: ', elm_len)
				elif elm_len < tiny_max_len:
					file = file_tiny
				elif elm_len < short_max_len:
					file = file_short	
				else:
					file = file_long

				left_ltr_seq = str(ltr_seq[int(k[3]):int(k[4])]).lower()
				if not 'n' in left_ltr_seq:
					pos = ltr_chr + ":" + k[3] + "-" + k[4]
					lt_ltr_header = '>' +  pos + ':' + ltr_group  
					file.write(lt_ltr_header + '\n')
					file.write(left_ltr_seq + '\n')
					header_file_dict[pos] = nick_name
					ltr_id_dict[pos] = ltr_counter

				# Process right LTR
				elm_len = int(k[6]) - int(k[5])
				if elm_len < min_len:
					print('Right LTR is too small: ', elm_len)
				elif elm_len < tiny_max_len:
					file = file_tiny
				elif elm_len < short_max_len:
					file = file_short	
				else:
					file = file_long

				right_ltr_seq = str(ltr_seq[int(k[5]):int(k[6])]).lower()
				if not 'n' in right_ltr_seq:
					pos = ltr_chr + ":" + k[5] + "-" + k[6]
					rt_ltr_header = '>' +  pos + ':' + ltr_group
					file.write(rt_ltr_header + '\n')
					file.write( right_ltr_seq +'\n')
					header_file_dict[pos] = nick_name
					ltr_id_dict[pos] = ltr_counter

				# Process interior of LTR element
				int_start = int(k[4]) + 1
				int_end   = int(k[5]) - 1
				
				elm_len = int_end - int_start
				if elm_len < min_len:
					print('LTR interior is too small: ', elm_len)
				elif elm_len < tiny_max_len:
					file = file_tiny
				elif elm_len < short_max_len:
					file = file_short	
				else:
					file = file_long

				# Sometimes an LTR RT has no interior; it may be a false positive
				interior_seq = str(ltr_seq[int_start:int_end]).lower()
				# if not 'n' in interior_seq: 
				#if int_end > int_start:
				if (not 'n' in interior_seq) and elm_len >= min_len:
					pos = ltr_chr + ":" + str(int_start) + "-" + str(int_end) 
					int_header = '>' + pos + ':' + interior_group
					file.write(int_header + '\n')
					file.write(interior_seq +'\n')
					header_file_dict[pos] = nick_name
					ltr_id_dict[pos] = ltr_counter
				# else:
				# 	print('Skipping a predicted LTR element with a too-small interior.')	
				# else:
				# 	print('The sequence has unknown nucleotides –– skipped.')		

				ltr_counter += 1

		aligner = Align.PairwiseAligner()
		aligner.mode = 'local'
		aligner.open_gap_score = -1.0
		aligner.extend_gap_score = -1.0
		aligner.mismatch = -1.0

		if len(g) > 0:
			id  = g[0][0]
			seq = seq_dict[id].seq
			for n in g[0:]:
				if(id != n[0]):
					id  = n[0]
					seq = seq_dict[n[0]].seq

				start   = int(n[1])	
				end     = int(n[2])
				seq_len = end - start

				if seq_len < min_len:
					continue

				if seq_len < tiny_max_len:
					file = file_tiny
				elif seq_len < short_max_len:
					file = file_short
				else:
					file = file_long

				# chr:start-end:file
				# When a file contains more than one sequence, the file name is needed.	
				seq_header = '>' + n[0] + ':' + n[1] + '-' + n[2]
				
				header_file_dict[seq_header[1:]] = nick_name
				
				# Find putative DNA transposons
				rc = seq[end-dna_win:end].reverse_complement()
				score = aligner.score(str(seq[start:start+dna_win]), str(rc))

				if score >= dna_min_score:
					file.write(seq_header + ':' + dna_group + '\n')
					file.write(str(seq[start:end]).lower() + '\n')
					continue
				#
				# Find putative SINEs and LINEs
				#
				# Search for poly-A tail in the positive strand
				score = aligner.score(str(seq[end-poly_a_win:end]), poly_a_win * "a")
				if score < poly_a_min_score:
					# Search for poly-A tail in the negative strand
					score = aligner.score(str(seq[start:start+poly_a_win]), poly_a_win * "t")	

				if score >= poly_a_min_score:
					file.write(seq_header + ':' + interspersed_group + '\n')
					file.write(str(seq[start:end]).lower() + '\n')
					continue

				# Search for Palindrome	
				palindrome = find_longest_palindrome(str(seq[end-helitron_win:end]))
				if palindrome == ():
					score = 0
				else:	
					score = palindrome[1] - palindrome[0] + 1
			
				if score < helitron_min_score:
					palindrome = find_longest_palindrome(str(seq[start:start+helitron_win]))
					if palindrome == ():
						score = 0
					else:
						score = palindrome[1] - palindrome[0] + 1		

				if score >= helitron_min_score:
					file.write(seq_header + ':' + helitron_group + '\n')
					file.write(str(seq[start:end]).lower() + '\n')
					continue

				file.write(seq_header + ':' + unknown_group + '\n')
				file.write(str(seq[start:end]).lower()  + '\n')	

		file_tiny.close()
		file_short.close()
		file_long.close()		

def remove_empty_files(dir):
	# Remove 0-size files
	for nick_name in os.listdir(dir):
		full_name = dir + '/'+ nick_name
		if os.path.getsize(full_name) == 0:
			os.remove(full_name)
			print('Deleting empty file: ', full_name)

def run_meshclust():
	# MeShClust cannot process empty files
	remove_empty_files(out_dir)
	
	# Start with a new file
	# all_out_file       = open(all_out, 'w')
	all_cluster_list = []

	#for group in ['short']:
	for group in ['tiny', 'short', 'long']:
		print('Running MeShClust on ', group , ' ...')
		meshclust_out = out_dir + '/' + group + '.cluster'
		cmd = 'meshclust --id %f %s/*%s.fa --output %s' %(meshclust_id, out_dir, group, meshclust_out)
		print(cmd)
		execute(cmd)
		# all_out_file.write('\n' + group + '\n')
		cluster_maker = ClusterList(meshclust_out, meshclust_min_size, group)
		# cluster_list.print(all_out_file)

		print(group, ' number of clusters: ', len(cluster_maker.cluster_list))
		# Add the found centers to the list 
		for cluster in cluster_maker.cluster_list:
			c = cluster.center
			header = c.chr + ':' + str(c.start) + '-' + str(c.end)
			c.file_nick_name = header_file_dict[header]
			all_cluster_list.append(cluster)
	
	# all_out_file.close()
	all_cluster_list.sort(key = lambda cluster:cluster.center.file_nick_name)
	
	if len(all_cluster_list) > 0:
		# all_fasta_file = open(all_fasta, 'w')
		nick_name = all_cluster_list[0].center.file_nick_name
		seq_dict = SeqIO.index(in_dir + '/' + nick_name + '.fa', "fasta")
		chr = all_cluster_list[0].center.chr
		seq = seq_dict[chr].seq
		
		for cluster in all_cluster_list:
			# Load the FASTA file
			if cluster.center.file_nick_name != nick_name:
				nick_name = cluster.center.file_nick_name
				seq_dict = SeqIO.index(in_dir + '/' + nick_name + '.fa', "fasta")	
			# Load the sequence	
			if cluster.center.chr != chr:
				chr = cluster.center.chr
				seq = seq_dict[chr].seq

			cluster.set_center_seq(seq[cluster.center.start:cluster.center.end].lower())
	
	print('Number of clusters before merge:', len(all_cluster_list))			
	merge_clusters(all_cluster_list)
	print('Number of clusters after merge:', len(all_cluster_list))

	# Assign a transposon group/familty to a cluster
	for cluster in all_cluster_list:
		cluster.assign_label()

	# Print the final library	
	print_library(all_cluster_list)

	# Match LTR cluster to its interior cluster
	match_ltr_to_interior(all_cluster_list)

def calculate_identity(alignment):
	"""
	Calculates the identity between two aligned seqeucnes.
	The parameter is a tuple of the form (seq1, seq2, score, alignment_start, alignment_end).
	Returns a ratio —— not a percentage.
	"""
	seq1, seq2, align_start, align_end = alignment[0], alignment[1], alignment[-2], alignment[-1]
	match_count = 0.0
	for i in range(align_start, align_end):
		if seq1[i] == seq2[i]:
			match_count += 1.0
	return match_count / (align_end - align_start)

def merge_clusters(cluster_list):
	"""
	Merges clusters if their centers are similar. It accounts for similarity to
	the revese complement as well.
	"""
	can_merge_list = [True] * len(cluster_list)

	for i in range(0, len(cluster_list)-1):
		seq_i_len = len(cluster_list[i].seq)
		#print('.', end='')
		for j in range(i+1, len(cluster_list)):
			lower_limit = meshclust_id * seq_i_len
			upper_limit = (1 - meshclust_id) * seq_i_len + seq_i_len
			seq_j_len = len(cluster_list[j].seq)
			if can_merge_list[j] and seq_j_len >= lower_limit and seq_j_len <= upper_limit:
				print('Aligning ', str(seq_i_len), ' bp vs ', str(seq_j_len), ' bp')
				alignment_list = pairwise2.align.globalms(cluster_list[i].seq, cluster_list[j].seq, 1, -3, -5, -2)
				if calculate_identity(alignment_list[0]) >= meshclust_id:
					cluster_list[i].merge(cluster_list[j])
					can_merge_list[j] = False
				else:
					seq_j_rc = str(cluster_list[j].seq.reverse_complement())
					alignment_list = pairwise2.align.globalms(cluster_list[i].seq, seq_j_rc, 1, -3, -5, -2)
					if calculate_identity(alignment_list[0]) >= meshclust_id:
						cluster_list[i].merge(cluster_list[j])
						can_merge_list[j] = False		
	#print('\n')

	# Delete merged clusters					
	for i in reversed(range(0, len(can_merge_list))):
		if not can_merge_list[i]:
			del cluster_list[i]	


def print_library(cluster_list):
	"""
	Name each family with a unique number and the composition of its groups
	"""
	library_file = open(library, 'w')
	counter = 1
	for cluster in cluster_list:
		name = '>Familiy' + str(counter) + '_' + cluster.get_detailed_label()
		cluster.set_name(name)
		library_file.write(name + '\n')
		library_file.write(str(cluster.seq) + '\n')
		counter += 1
	library_file.close()

def calculate_cluster_intersection(c1, c2):
	"""
	Calculate the intersection between the LABELS between two cluster
	It does NOT claculate the intersection based on membership
	It is designed to be used in matching LTR to its interior
	"""
	s1 = set()
	for m in c1.member_list:
		pos = m.chr + ':' + str(m.start) + '-' + str(m.end) 
		if pos in ltr_id_dict.keys():
			s1.add(ltr_id_dict[pos]) 

	s2 = set()
	for m in c2.member_list:
		pos = m.chr + ':' + str(m.start) + '-' + str(m.end) 
		if pos in ltr_id_dict.keys():
			s2.add(ltr_id_dict[pos])

	i = s1 & s2		
	return len(i)

def match_ltr_to_interior(cluster_list):
	# Make a list of clusters predominantly LTR 
	# Make another list of clusters predominantly interiors
	ltr_list = []
	int_list = []	
	for c in cluster_list:
		label = c.get_label()
		if label == ltr_group:
			ltr_list.append(c)
		elif label == interior_group:
			int_list.append(c)

	# Match an LTR cluster to an Interior cluster	
	ltr_to_interior_file = open(ltr_to_interior, 'w')	
	for ltr_c in ltr_list:		
		max_i = 0
		max_c = None
		for int_c in int_list:
			i = calculate_cluster_intersection(ltr_c, int_c)
			if i > max_i:
				max_i = i
				max_c = int_c
		if max_c == None:
			print('No match found for', ltr_c.name)
		else:
			ltr_interior_dict[ltr_c] = max_c
			ltr_to_interior_file.write(ltr_c.name + ' --> ' + max_c.name + '\n')

	ltr_to_interior_file.close()		

def run_fastcar():
	# fastcar ALL.fa --all-vs-all --id 0.6 --datatype uint16_t --output fastcar_out
	output_root = out_dir + '/fastcarPart'
	cmd = 'fastcar %s --all-vs-all --id %f --output %s' %(all_fasta, fastcar_id, output_root)
	execute(cmd)
	
	# Collect FASTCAR's output files into one
	part_list = glob.glob(output_root + '*')
	file = open(all_vs_all, 'w')
	for part in [open(x, 'r') for x in part_list]:
		file.write(part.read())
	file.close()

def generate_final_tree():	
	file = open(all_vs_all, 'r')
	# List of labels
	label_dict = {}
	counter = 0
	for line in file:
		token_list = line.split()
		if not token_list[0] in label_dict:
			label_dict[token_list[0]] = counter
			counter += 1
		if not token_list[1] in label_dict:	
			label_dict[token_list[1]] = counter
			counter += 1

	label_count = len(label_dict.keys())

	# The matrix 
	matrix = np.zeros((label_count, label_count))
	# Go to the begining of the file
	file.seek(0)
	for line in file:
		token_list = line.split()
		distance = 100 - float(token_list[2])
		id_1 = label_dict[token_list[0]]
		id_2 = label_dict[token_list[1]]
		matrix[id_1][id_2] = distance
		matrix[id_2][id_1] = distance
	file.close()
	# Make a hierarchical cluster dendrogram 
	row_clusters = linkage(matrix, method='average')
	row_dendrogram = dendrogram(row_clusters, labels=list(label_dict.keys()), orientation='left')
	plt.tight_layout()
	plt.show()

def test_mite():
	aligner = Align.PairwiseAligner()
	aligner.mode = 'local'
	aligner.open_gap_score = -1.0
	aligner.extend_gap_score = -1.0
	aligner.mismatch = -1.0

	for ir_min_score in range(5,11): 
		for lt_win in range(10, 50, 5):
			dna_counter = 0
			tp_counter = 0
			fp_counter = 0
			
			rt_win = lt_win

			for record in SeqIO.parse("/Users/zakarota/Tools/repeatMasker/Libraries/RepeatMaskerLib.embl","embl"):
				mark = 'Type: DNA'
				if mark in record.annotations['comment']:
					dna_counter += 1

				rc = record.seq[-lt_win:-1].reverse_complement()
				score = aligner.score(str(record.seq[0:rt_win]), str(rc))
				# print(score)

				if score >= ir_min_score and mark in record.annotations['comment']:
					tp_counter += 1

				if score >= ir_min_score and not mark in record.annotations['comment']:
					fp_counter += 1

			sensitivity = 100 * tp_counter / dna_counter
			if tp_counter > 0 or fp_counter > 0:
				precision   = 100 * tp_counter / (tp_counter + fp_counter)
			else:
				precision = 0	

			if sensitivity > 0 or precision > 0:
				f_measure   = 2 * (sensitivity * precision) / (sensitivity + precision)	
			else:
				f_measure = 0
			
			print('Left window:', lt_win, 'Right window:', rt_win, 'Min score:', ir_min_score, 'Sensitivity:', sensitivity, 'Precision:', precision, 'F-measure:', f_measure)


def test_poly_a():

	for poly_a_min_score in range(5,11):
		for poly_a_win in range(5,31,5):

			pstv_counter = 0
			tp_counter   = 0
			fp_counter   = 0	

			for record in SeqIO.parse("/Users/zakarota/Tools/repeatMasker/Libraries/RepeatMaskerLib.embl","embl"):
				mark1 = 'SINE'
				mark2 = 'LINE'
				if mark1 in record.annotations['comment'] or mark2 in record.annotations['comment']:
					pstv_counter += 1

				aligner = Align.PairwiseAligner()
				aligner.mode = 'local'
				aligner.open_gap_score = -1.0
				aligner.extend_gap_score = -1.0
				aligner.mismatch = -1

				score = aligner.score(record.seq[-poly_a_win:-1], poly_a_win * "a")
				if score < poly_a_min_score:
					# Search for poly-A tail in the negative strand
					score = aligner.score(str(record.seq[0:poly_a_win]), poly_a_win * "t")	

				if score >= poly_a_min_score and (mark1 in record.annotations['comment'] or mark2 in record.annotations['comment']):
					tp_counter += 1

				if score >= poly_a_min_score and not (mark1 in record.annotations['comment'] or mark2 in record.annotations['comment']):
					fp_counter += 1

			sensitivity = 100 * tp_counter / pstv_counter
			if tp_counter > 0 or fp_counter > 0:
				precision   = 100 * tp_counter / (tp_counter + fp_counter)
			else:
				precision = 0	

			if sensitivity > 0 or precision > 0:
				f_measure   = 2 * (sensitivity * precision) / (sensitivity + precision)	
			else:
				f_measure = 0

			print('Window:', poly_a_win, 'Min score:', poly_a_min_score, 'Sensitivity:', sensitivity, 'Precision:', precision, 'F-measure:', f_measure)


def has_CTRR_pstv(s):
	s = s.lower()
	if ('ctaa' in s) or ('ctgg' in s) or ('ctag' in s) or ('ctga' in s):
		return True
	else:
		return False

def has_CTRR_ngtv(s):
	s = s.lower()
	if ('ttag' in s) or ('ccag' in s) or ('ctag' in s) or ('tcag' in s):
		return True
	else:
		return False

				
def test_hlitron():
	helitron_win = 50
	helitron_min_score = 20

	pstv_counter = 0
	tp_counter   = 0
	fp_counter   = 0	

	aligner = Align.PairwiseAligner()
	aligner.mode = 'local'
	aligner.open_gap_score = -1.0
	aligner.extend_gap_score = -1.0
	aligner.mismatch = -1

	for record in SeqIO.parse("/Users/zakarota/Tools/repeatMasker/Libraries/RepeatMaskerLib.embl","embl"):
		mark = 'Helitron'

		if len(record.seq) < 3000:
			continue

		if mark in record.annotations['comment']:
			pstv_counter += 1

		# Skip if it is LTR
		if 'LTR' in record.annotations['comment']:
			continue

		# Skip if it is a DNA transposon	
		rc = record.seq[-dna_win:-1].reverse_complement()
		score = aligner.score(str(record.seq[0:dna_win]), str(rc))
		if score >= dna_min_score:
			continue;

		# Skip if it has Poly-A tail
		score = aligner.score(record.seq[-poly_a_win:-1], poly_a_win * "a")
		if score < poly_a_min_score:
			score = aligner.score(str(record.seq[0:poly_a_win]), poly_a_win * "t")	
		if(score >= poly_a_min_score):
			continue

		# # Search for the palindrome in the positive strand
		# pstv_strand = record.seq[-helitron_win:-1]
		# score = aligner.score(str(pstv_strand), str(pstv_strand.reverse_complement()))
		# if score < helitron_min_score:
		# 	# Search for the palindrom in the negative strand
		# 	ngtv_strand = record.seq[0:helitron_win]
		# 	score = aligner.score(str(ngtv_strand), str(ngtv_strand.reverse_complement()))
		palindrome = find_longest_palindrome(str(record.seq[-helitron_win:-1]))
		if palindrome == ():
			score = 0
			ctrr  = False
		else:	
			score = palindrome[1] - palindrome[0] + 1
			ctrr  = has_CTRR_pstv(record.seq[-(helitron_win - palindrome[1] - 1):-1])
		
		if score < helitron_min_score:
			palindrome = find_longest_palindrome(str(record.seq[0:helitron_win]))
			if palindrome == ():
				score = 0
				ctrr  = False
			else:
				score = palindrome[1] - palindrome[0] + 1
				ctrr = has_CTRR_ngtv(record.seq[0:palindrome[0]])			

		if score >= helitron_min_score and mark in record.annotations['comment']:
			tp_counter += 1

		if score >= helitron_min_score and not (mark in record.annotations['comment']):
			fp_counter += 1
			print(record.id)

	sensitivity = 100 * tp_counter / pstv_counter
	if tp_counter > 0 or fp_counter > 0:
		precision   = 100 * tp_counter / (tp_counter + fp_counter)
	else:
		precision = 0	

	if sensitivity > 0 or precision > 0:
		f_measure   = 2 * (sensitivity * precision) / (sensitivity + precision)	
	else:
		f_measure = 0

	print('Window:', helitron_win, 'Min score:', helitron_min_score, 'Sensitivity:', sensitivity, 'Precision:', precision, 'F-measure:', f_measure)

def is_complement(a, b):
	result = False
	if (a == 'a' and b == 't') or (a == 't' and b == 'a') or (a == 'c' and b == 'g') or (a == 'g' and b == 'c'):
		result = True
	return result	

def find_longest_palindrome(s):
	s = s.lower()

	max_len = -1
	max_cor = ()
	for i in range(0, len(s)-1):
		if is_complement(s[i], s[i+1]):
			lt = i
			rt = i + 1

			if i > 0:
				for j in reversed(range(0, i)):
					j_prime = i+(i-j)+1
					if  j_prime < len(s) and is_complement(s[j], s[j_prime]):
						lt = j
						rt = j_prime
					else:
						break

			if rt - lt + 1 > max_len:
				max_len = rt - lt + 1
				max_cor = (lt, rt)

	# if max_cor == ():
	# 	return ''
	# else:			
	# 	# return s[max_cor[0]:max_cor[1]+1]				
	return max_cor

# Start the pipeline
format_input()
run_Red()
run_LtrDetector()
run_Look4TRs()
group_based_on_structure()
run_meshclust()