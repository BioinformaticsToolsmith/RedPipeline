#
# Author: Hani Z. Girgis, PhD
# The Bioinformatics Toolsmith Laboratory
# 5/27/2019
#

import re
import sys

class Member():
	"""An instance of this class keeps information about one member"""
	def __init__(self, chr='', start=0, end=0, group='', file_nick_name=''):
		self.chr   = chr
		self.start = start
		self.end   = end
		self.group = group
		self.file_nick_name = file_nick_name

	def __str__(self):
		return self.chr + ':' + str(self.start) + '-' + str(self.end) + ':' + self.group	

class Cluster():
	"""An instance of this class keeps information about one cluster"""
	def __init__(self, name='', center=None, member_list=[], note=''):
		self.name        = name
		self.center      = center
		self.member_list = member_list
		self.note        = note
		self.lable       = 'UNDEFINED'

	def check_center(self):
		"""The center is the median member based on lenght if MeShClust does not assign a center"""
		if self.center == None:
			# Sort members according to length
			self.member_list.sort(key = lambda i:i.end - i.start)
			mid_index = round( (len(self.member_list)-1) / 2.0 )
			if mid_index >= 0:
				self.center = self.member_list[mid_index]
			else:
				raise Exception('Something went wrong while assigning a missing center.')
			print('Assinged a missing center.')

	# Sets the sequence associated of the center of this cluster		
	def set_center_seq(self, seq):
		self.seq = seq

	def merge(self, other):
		self.member_list.extend(other.member_list)

	def assign_label(self):
		self.group_dict = {}
		for member in self.member_list:
			if member.group in self.group_dict.keys():
				self.group_dict[member.group] += 1.0
			else:
				self.group_dict[member.group] = 1.0

		total = 0.0	
		for count in self.group_dict.values():
			total += count

		self.label   = ''
		highest = -1.0
		for group in self.group_dict.keys():
			self.group_dict[group] = 100 * self.group_dict[group]/total
			if self.group_dict[group] > highest:
				self.label   = group
				highest = self.group_dict[group] 	

	# Make a label of all families along with their percentage
	def get_detailed_label(self):	
		sorted_groups = sorted(self.group_dict.keys(), key= lambda gp:self.group_dict[gp])
		label_str = ''
		for group in reversed(sorted_groups):
			label_str += group + ':' + str( round(self.group_dict[group],2) ) + '_' 
		return label_str[0:-1]

	def get_label(self):
		return self.label

	def set_name(self, name):
		self.name = name	
		
	def __str__(self):
		s = self.name + '\t' + self.note + '\t' + str(self.center) + '\n'
		for item in self.member_list:
			s += '\t' + str(item) + '\n'
		return s

class ClusterList():
	"""An instance of this class handles MeShClust output"""
	def __init__(self, file_name, threshold=10, note=''):
		self.file_name = file_name
		self.threshold = threshold
		self.cluster_list = []
		self.note = note
		self.parse()

	# >Cluster 0
	# 0       50nt, >chr1:2895421-2895471... *
	# >Cluster 1
	# 0       50nt, >chr5:25848461-25848511... *
	# >Cluster 2
	# 0       50nt, >chr2:10299809-10299859... *
	# 1       50nt, >chr4:7762286-7762336... 
	# 2       50nt, >chr5:8415559-8415609... 
	# 3       50nt, >chr5:14210973-14211023... 
	# >Cluster 3
	# 0       50nt, >chr1:4692758-4692808... *
	# >Cluster 4
	# 0       50nt, >chr2:19256993-19257043... *
	# >Cluster 5
	# 0       50nt, >chr1:5156886-5156936... *
	# >Cluster 6
	# 0       50nt, >chr3:1931906-1931956... *
	def parse(self):
		file = open(self.file_name, 'r')
		cluster = Cluster()
		is_first = True
		for line in file:
			if line[0] == '>':
				if is_first:
					is_first = False
				else:
					if len(cluster.member_list) >= self.threshold:
						cluster.check_center()
						self.cluster_list.append(cluster)
					cluster = Cluster()

				cluster.name = line[1:]
				cluster.member_list = []
				cluster.note = self.note
			else:
				match = re.search(r'>(.+):(\d+)\-(\d+):(.+)\.\.\.', line)
				if match:
					# member = (match.group(1), int(match.group(2)), int(match.group(3)), match.group(4))
					# cluster.member_list.append(member)
					member = Member(match.group(1), int(match.group(2)), int(match.group(3)), match.group(4))
					cluster.member_list.append(member)
					if '*' in line:
						cluster.center = member

		# Handle last cluster
		if len(cluster.member_list) >= self.threshold:
			cluster.check_center()
			self.cluster_list.append(cluster)
		file.close()

	def print(self, file=sys.stdout):
		for cluster in self.cluster_list:
			file.write(cluster.name + '\t' + cluster.note + '\t' + str(cluster.center) + '\n')
			for item in cluster.member_list:
				file.write('\t' + str(item) + '\n')
				