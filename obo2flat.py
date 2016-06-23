#! /usr/bin/env python
# -*- coding: UTF-8 -*-

class extract_info(object):
	'''extract the specified information from the obo aanotation file.\n'''
	def __init__(self, options):
		self.obo_info = {}
		self.default_relation = [[], [], [], []]
		self.default_extract  = [[], [], [], [], [], [], 0, [], {}, []]
		self.obo_output_info = options.output + '/' + 'obo_info' + '/'
		self.obo_output_path = options.output + '/' + 'obo_path' + '/'
		self.info = ["alt_id:", "name:", "namespace:", "def:", "subset:", "is_a:",\
			     	"is_obsolete:", "replaced_by:", "relationship:", "consider:"]
		self.relation = ["regulates", "positively_regulates", "negatively_regulates", "part_of"]

def obo2info(options):
	'''The 'options' is the input parameter by your command, which contained OBO 
	annotation file absolutly path, output path and ect. '''
	fp_obo = open(options.input_obo, 'r')
	obo_all_info = fp_obo.read()
	fp_obo.close()
	carr = obo_all_info.split("[Term]\n")
	leng_carr = len(carr)
	info = extract_info(options)
	len_extract = len(info.info)

	for line in xrange(1,leng_carr-1):
		iterm_id = "GO:" + carr[line].split("GO:")[1].split()[0]
		term = carr[line].split("\n")
		if iterm_id in info.obo_info:
			sys.stderr.write("[Error] Duplicate go id.\n")
		else:
			info.obo_info[iterm_id] = copy.deepcopy(info.default_extract)
			info.obo_info[iterm_id][len_extract-2] = copy.deepcopy(info.default_relation)
		for line1 in term[1:]:
			if line1.startswith(info.info[0]):
				goid = 'GO:' + line1.split('GO:')[1].split()[0]
				info.obo_info[iterm_id][0].append(goid)
			if line1.startswith(info.info[1]):
				go_name = line1.split(info.info[1])[-1].strip()
				info.obo_info[iterm_id][1] = go_name
			if line1.startswith(info.info[2]):
				go_namespace = line1.split(info.info[2])[-1].strip()
				info.obo_info[iterm_id][2] = go_namespace
			if line1.startswith(info.info[3]):
				go_def = line1.split(info.info[3])[-1].strip()
				info.obo_info[iterm_id][3] = go_def
			if line1.startswith(info.info[4]):
				go_subset = line1.split(info.info[4])[-1].strip()
				info.obo_info[iterm_id][4].append(go_subset)
			if line1.startswith(info.info[5]):
				go_isa = 'GO:' + line1.split('GO:')[1].split()[0]
				info.obo_info[iterm_id][5].append(go_isa)
			if line1.startswith(info.info[6]):
				info.obo_info[iterm_id][6] = 1
			if line1.startswith(info.info[7]):
				go_replace = 'GO:' + line1.split('GO:')[1].split()[0]
				info.obo_info[iterm_id][7].append(go_replace)
			if line1.startswith(info.info[8]):
				go_relation = line1.split(info.info[8])[-1].split('!')[0].strip()
				for line2 in info.relation:
					if go_relation.startswith(line2):
						index = info.relation.index(line2)
						go_sub_relation = go_relation.split(line2)[1].strip()
						info.obo_info[iterm_id][8][index].append(go_sub_relation)
			if line1.startswith(info.info[9]):
				go_consider = line1.split(info.info[9])[-1].strip()
				info.obo_info[iterm_id][9].append(go_consider)
	return info.obo_info
	

def output_obo_info(options, obo_info):
	'''Output useful info in OBO info, object represented the obo infomation which had been extracted, 
	obo_output_directory is the output directory of the useful information which mentioned above'''
	info = extract_info(options)
	if os.path.exists(info.obo_output_info):
		pass
	else:
		os.mkdir(info.obo_output_info)
	
	fp_alt = file(info.obo_output_info  + 'GO_alt_id' + '.alt_id', 'w')
	fp_isa = file(info.obo_output_info  + 'GO_isa' + '.isa', 'w')
	fp_part_of = file(info.obo_output_info  + 'GO_part_of' + '.part_of', 'w')
	fp_regulates = file(info.obo_output_info  + 'GO_regulates' + '.regulates', 'w')
	fp_pos_reg = file(info.obo_output_info  + 'GO_positive_regulates' + '.positive_regulates', 'w')
	fp_neg_reg = file(info.obo_output_info  + 'GO_negative_regulates' + '.negative_regulates', 'w')
	fp_def = file(info.obo_output_info  + 'GO_def' + '.def', 'w')
	fp_consider = file(info.obo_output_info  + 'GO_consider' + '.consider', 'w')
	fp_replaced = file(info.obo_output_info  + 'GO_replace' + '.replace', 'w')
	fp_bp = file(info.obo_output_info  + 'GO_bp' + '.bp', 'w')
	fp_mf = file(info.obo_output_info  + 'GO_mf' + '.mf', 'w')
	fp_cc = file(info.obo_output_info  + 'GO_cc' + '.cc', 'w')
	
	fp_alt.write('#GO_id\tAlt_id\n')
	fp_isa.write('#GO_id\tIsa\n')
	fp_part_of.write('#GO_id\tPart_of\n')
	fp_regulates.write('#GO_id\tRegulates\n')
	fp_pos_reg.write('#GO_id\tPositive_regulates\n')
	fp_neg_reg.write('#GO_id\tNegative_regulates\n')
	fp_def.write('#GO_id\tName\tNamespace\tdef\n')
	fp_consider.write('#GO_id\tConsider\n')
	fp_replaced.write('#GO_id\tReplaced\n')
	fp_bp.write("#GO_id\tName\tNamespace\tdef\n")
	fp_mf.write("#GO_id\tName\tNamespace\tdef\n")
	fp_cc.write("#GO_id\tName\tNamespace\tdef\n")
	
	len_info = len(info.default_extract)
	len_relation = len(info.relation)
	slim = {}
	for iterm in obo_info:
		for line in xrange(len_info):
			if line == 0:
				if obo_info[iterm][line]:
					fp_alt.write(iterm + '\t')
					for iterm1 in obo_info[iterm][line]:
						if iterm1 not in obo_info[iterm][line][-1]:
							fp_alt.write(iterm1 + '|')
						else:
							fp_alt.write(iterm1 + '\n')
			if line == 3:
				fp_def.write(iterm + '\t' + '\t'.join([obo_info[iterm][line-2], obo_info[iterm][line-1], obo_info[iterm][line]]) + '\n')
				if obo_info[iterm][line-1] in 'biological_process':
					fp_bp.write(iterm + '\t' + '\t'.join([obo_info[iterm][line-2], obo_info[iterm][line-1], obo_info[iterm][line]]) + '\n')
				elif obo_info[iterm][line-1] in 'molecular_function':
					fp_mf.write(iterm + '\t' + '\t'.join([obo_info[iterm][line-2], obo_info[iterm][line-1], obo_info[iterm][line]]) + '\n')
				elif obo_info[iterm][line-1] in 'cellular_component':
					fp_cc.write(iterm + '\t' + '\t'.join([obo_info[iterm][line-2], obo_info[iterm][line-1], obo_info[iterm][line]]) + '\n')
			if line == 4:
				arrtr = obo_info[iterm][line]
				for go_slim in arrtr:
					if go_slim in slim:
						slim[go_slim].append(iterm)
					else:
						slim[go_slim] = [iterm,]
			if line == 5:
				if obo_info[iterm][line]:
					fp_isa.write(iterm + '\t')
					for iterm1 in obo_info[iterm][line]:
						if iterm1 not in obo_info[iterm][line][-1]:
							 fp_isa.write(iterm1 + '|')
						else:
							 fp_isa.write(iterm1 + '\n')
			if line == 6:
				if obo_info[iterm][line] == 1:
					if obo_info[iterm][-1]:
						fp_consider.write(iterm + '\t')
						for iterm1 in obo_info[iterm][-1]:
							if iterm1 not in obo_info[iterm][-1][-1]:
								fp_consider.write(iterm1 + '|')
							else:
								fp_consider.write(iterm1 + '\n')
			if line == 7:
				if obo_info[iterm][line]:
					fp_replaced.write(iterm + '\t')
					for iterm1 in obo_info[iterm][line]:
						if iterm1 not in obo_info[iterm][line][-1]:
							fp_replaced.write(iterm1 + '|')
						else:
							fp_replaced.write(iterm1 + '\n')
			if line == 8:
				for j in xrange(len_relation):
					if obo_info[iterm][8][j]:
						if j == 0:
							fp_regulates.write(iterm + '\t')
							for iterm1 in obo_info[iterm][8][j]:
								if iterm1 not in obo_info[iterm][8][j][-1]:
									fp_regulates.write(iterm1 + '|')
								else:
									fp_regulates.write(iterm1 + '\n')
						elif j == 1:
							fp_pos_reg.write(iterm + '\t')
							for iterm1 in obo_info[iterm][8][j]:
								if iterm1 not in obo_info[iterm][8][j][-1]:
									fp_pos_reg.write(iterm1 + '|')
								else:
									fp_pos_reg.write(iterm1 + '\n')
						elif j == 2:
							fp_neg_reg.write(iterm + '\t')
							for iterm1 in obo_info[iterm][8][j]:
								if iterm1 not in obo_info[iterm][8][j][-1]:
									fp_neg_reg.write(iterm1 + '|')
								else:
									fp_neg_reg.write(iterm1 + '\n')
						elif j == 3:
							fp_part_of.write(iterm + '\t')
							for iterm1 in obo_info[iterm][8][j]:
								if iterm1 not in obo_info[iterm][8][j][-1]:
									fp_neg_reg.write(iterm1 + '|')
								else:
									fp_neg_reg.write(iterm1 + '\n')
	for subset in slim:
		fp = file(info.obo_output_info + "obo" + "_" + subset + "." + "slim", 'w')
		fp.write("\n".join(slim[subset]))
		fp.close()
	
	fp_alt.close()
	fp_isa.close()
	fp_part_of.close()
	fp_regulates.close()
	fp_pos_reg.close()
	fp_neg_reg.close()
	fp_def.close()
	fp_consider.close()
	fp_replaced.close()
	fp_bp.close()
	fp_cc.close()
	fp_mf.close()

class define_info(object):
	'''Define valiables for searching topology constructure and paths, and add the flag of each edge'''
	def __init__(self, options):
		self.obo_toplogy = {}
		self.relation = {}
		self.go_class = {}
		self.flag = {}
		self.flat = {}
		obo = extract_info(options)
		self.relations = ["is_a", "regulates", "positively_regulates", "negatively_regulates", "part_of"]
		self.label = ["i", "r", "p", "n", "o"]
		self.save_path = obo.obo_output_path
		self.default = obo.default_extract
		self.go_class["GO:0008150"] = "biological_process"
		self.go_class["GO:0003674"] = "molecular_function"
		self.go_class["GO:0005575"] = "cellular_component"

def add_label(options, obo_info):
	'''Save the relationship of each go semantic and add the label for recognition.\n'''
	def_info = define_info(options)	

	for iterm in obo_info:
		for index in [5, 8]:
			if index == 5 and obo_info[iterm][index]:
				def_info.relation[iterm] = copy.deepcopy([])
				def_info.flag[iterm] = copy.deepcopy([])
				def_info.relation[iterm] = def_info.relation[iterm] + obo_info[iterm][index]
				def_info.flag[iterm] = def_info.flag[iterm] + list(def_info.label[0]*len(obo_info[iterm][index]))
			elif index == 8 and obo_info[iterm][index]:
				if obo_info[iterm][5]:
					pass
				else:
					def_info.relation[iterm] = copy.deepcopy([])
					def_info.flag[iterm] = copy.deepcopy([])
				'''The integer 4 represents the length of ["regulates", "positively_regulates", "negatively_regulates", "part_of"]'''
				for i in xrange(4):
					if obo_info[iterm][index][i]:
						def_info.relation[iterm] = def_info.relation[iterm] + obo_info[iterm][index][i]
						def_info.flag[iterm] = def_info.flag[iterm] + list(def_info.label[i + 1]*len(obo_info[iterm][index][i]))
	return def_info.relation, def_info.flag

class goflat(object):
	def __init__(self, obo_relation, obo_flag, go_id):
		self.h = obo_relation
		self.f = obo_flag
		self.goid = go_id
		self.goarray = []
		self.flagarray = []
		self.flat = []
		try:
			self.newarray = self.h[go_id]
			self.newflag = self.f[go_id]
			index = 0
			for line in self.h[go_id]:
				self.goarray.append(tuple([go_id, line]))
				#index = self.h[go_id].index(line)
				self.flagarray.append(self.f[go_id][index])		
				index = index + 1
				if line not in self.flat:
					self.flat.append(line)
		except:
			self.newarray = []
	def increase(self):
		tmplist = []
		for go_id in self.newarray:
			if go_id in self.h:
				tmplist = tmplist + self.h[go_id]
				tmplist = list(set(tmplist))
				index = 0
				for line in self.h[go_id]:
					self.goarray.append(tuple([go_id, line]))
					#index = self.h[go_id].index(line)
					self.flagarray.append(self.f[go_id][index])
					index = index + 1
					if line not in self.flat:
						self.flat.append(line)
		if tmplist:
			self.newarray = tmplist[:]
			self.increase()
		else:
			return 0

def find_all_paths(graph, start, end, path = []):
	path = path + [start]
	if start == end:
		return path
	if not graph.has_key(start):
		return []
	paths = []
	for node in graph[start]:
		if node not in path:
			newpaths = find_all_paths(graph, node, end, path)
			for newpath in newpaths:
				paths.append(newpath)
	return paths

def find_a_path(graph, start, end, path = []):
	path = path + [start]
	if start == end:
		return path
	if not graph.has_key(start):
		return None
	for node in graph[start]:
		if node not in path:
			newpath = find_a_path(graph, start, end, path)
			if newpath:
				return newpath
	return None

def find_shortest_path(graph, start, end, path = []):
	path = path + [start]
	if start == end:
		return path
	if graph.has_key(start):
		return None
	shortest_path = None
	for node in graph[start]:
		if node not in path:
			newpath = find_shortest_path(graph, node, end, path = [])
			if newpath:
				if not shortest or len(newpath) < len(shortest):
					shortest_path = newpath
	return shortest_path

def search_topology(obo_info, obo_relation, obo_flag, options):
	'''Search the topology structure of go semantic and save into a specified path.\n'''
	sys.stderr.write('[INFO] Generate Topology.\n')
	def_info = define_info(options)
	if os.path.exists(def_info.save_path):
		pass
	else:
		os.mkdir(def_info.save_path)
	
	fp_obo = file(def_info.save_path + "/obo_topology.txt", 'w')
	fp_obo.write("#go_id\tgo_name\tsearch_flag\tsearch_path\n")
	for go_id in obo_relation:
		search_path = goflat(obo_relation, obo_flag, go_id)
		search_path.increase()
		def_info.relation[go_id] = search_path.goarray
		def_info.flat[go_id] = search_path.flat
		def_info.flag[go_id] = search_path.flagarray
		fp_obo.write(go_id + '\t')
		fp_obo.write(''.join(obo_info[go_id][1]) + "\t")
		fp_obo.write(''.join(def_info.flag[go_id]) + '\t')
		tmp_path_leng = len(def_info.relation[go_id])
		if tmp_path_leng:
			for i in xrange(tmp_path_leng):
				if i < tmp_path_leng-1:
					fp_obo.write('->'.join(def_info.relation[go_id][i]) + ';')
				else:
					fp_obo.write('->'.join(def_info.relation[go_id][i]) + '\n')
		else:
			fp_obo.write('\t\n')
	fp_obo.close()
	return def_info.relation, def_info.flag, def_info.flat
		
def search_pathes(obo_info, relation, flag, flat, options):		
	fp = []
	def_info = define_info(options)
	for i in xrange(1, options.go_entire_path + 1):
		fp.append(file(def_info.save_path + 'obo_flat_%s.txt'%(str(i)), 'w'))
	fp_path = file(def_info.save_path + 'obo_flat_entire_path.txt', 'w')
	for go_id in relation:
		graph = {}
		tmp_path = list(set(relation[go_id]))
		if tmp_path:
			for i in xrange(len(tmp_path)):
				if tmp_path[i][0] not in graph.keys():
					graph[tmp_path[i][0]] = []
					graph[tmp_path[i][0]].append(tmp_path[i][1])
				else:
					graph[tmp_path[i][0]].append(tmp_path[i][1])
		index = [0] * (options.go_entire_path + 1)
		for go_class in def_info.go_class.keys():
			if go_class not in flat[go_id]: continue
			if not graph: break
			obo_path = []
			paths = find_all_paths(graph, go_id, go_class)
			for iterm in paths:
				if iterm in go_id:
					tmp_iterm = go_id
				elif iterm in go_class:
					tmp_iterm = tmp_iterm + '->' + iterm
					obo_path.append(tmp_iterm)
				else:
					tmp_iterm = tmp_iterm + '->' + iterm
			for line in obo_path:
				tmp_path = line.split('->')
				len_path = len(tmp_path)
				if index[-1] == 0:
					fp_path.write(go_id + '\t')
					fp_path.write(obo_info[go_id][1] + '\t')
					fp_path.write(obo_info[go_id][2] + '\t')
					index[-1] = index[-1] + 1
				fp_path.write(line + '|')
				if len_path <= options.go_entire_path + 1 and len_path >= 2:
					if index[len_path - 2] == 0:
						fp[len_path - 2].write(go_id + '\t')
						fp[len_path - 2].write(obo_info[go_id][1] + '\t')
						fp[len_path - 2].write(obo_info[go_id][2] + '\t')
						index[len_path - 2] = index[len_path - 2] + 1
					fp[len_path - 2].write(line + '|')
		for number in xrange(options.go_entire_path + 1):
			if index[number] > 0 and number < options.go_entire_path:
				fp[number].write('\n')
			elif number > 0 and number == options.go_entire_path:
				fp_path.write('\n')
	for i in xrange(options.go_entire_path):
		fp[i].close()
	fp_path.close()

def __main():
	usages = "USAGE: obo2flat [options]"
	'''Extracting the obo informatiom from OBO annotation format file, then generated the entire pathes and topology structure.'''
	optParser = OptionParser(usages)

	optParser_obb = OptionGroup(optParser, 'Obbligato options')
	optParser_obb.add_option('--run-obo2flat', dest = 'obo2flat', help = 'Run option for parsing obo annotation file', action = 'store_true')
	
	Parameter = OptionGroup(optParser, 'Obbligato parameter')
	Parameter.add_option('--input', dest = 'input_obo', help = 'Input obo annotation file', metavar = 'STR', type = 'string')
	Parameter.add_option('--output', dest = 'output', help = 'Output directory', metavar = 'DIR', type = 'string', default = '../obo')

	optional_Parameter = OptionGroup(optParser, 'Optional parameter')
	optional_Parameter.add_option('--path', dest = 'go_entire_path', help = 'Generate entire pathes and specified the number of pile pathes', metavar = 'INT', type = 'int')
	optParser.add_option_group(optParser_obb)
	optParser.add_option_group(Parameter)
	optParser.add_option_group(optional_Parameter)

	(options, argvs) = optParser.parse_args()
	if options.obo2flat:
		pass
	else:
		optParser.print_help()
		sys.stderr.write("[Error] Run option for obo annotion file has not been specfied.\n")
		exit(1)
	if options.input_obo == None:
		optParser.print_help()
		sys.stderr.write("[Error] Please input obo annotation file.\n")
		exit(1)
	else:
		if os.path.isfile(options.input_obo):
			pass
		else:
			optParser.print_help()
			sys.stderr.write("[Error] Please input correct obo annotation file path.\n")
			exit(1)
	return options

def call_all_functions():
	options = __main()
	sys.stderr.write('[INFO] Task Beginning.\n')
	obo_info = obo2info(options)
	output_obo_info(options, obo_info)
	obo_relation, obo_flag = add_label(options, obo_info)
	relation, flag, flat = search_topology(obo_info, obo_relation, obo_flag, options)	
	if options.go_entire_path != None:
		sys.stderr.write('[INFO] Output Semantic Pathes.\n')
		search_pathes(obo_info, relation, flag, flat, options)

if __name__ == '__main__':
	Packages = ['import sys', 'import os', 'import time',\
				'import copy', 'import string', 'import pdb',\
				'import networkx as nx', 'import numpy as np',\
				'import matplotlib.pyplot as plt', 'import pygraphviz as pgv',\
				'from pygraphviz import *', 'import pydot', 'import math',\
				'import scipy', 'from optparse import OptionParser, OptionGroup']
	for module in Packages:
		try:
			exec(module)
		except:
			print("[Error] No module named %s in the system, program exit.\n"%(module.split()[1]))
			exit(1)
	del module, Packages
	start_time = time.time()
	call_all_functions()
	final_time = time.time() - start_time
	sys.stderr.write('[CONSUME TIME] %s seconds.\n'%(str(final_time)))
	sys.stderr.write('[INFO] Task Done.\n')
	

	

