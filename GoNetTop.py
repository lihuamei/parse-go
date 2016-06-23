#!/usr/bin/env python

class obo_file_path():
	def __init__(self):
		self.go_class = {}
		self.godef = ""
		self.topology = ""
		self.altgo = ""
		self.repgo = ""
		self.divide = 0
		self.gene = ""
		self.qvalue = 0.05
		self.maxnum = 0
		self.output = ""
		self.go_class = ["GO:0008150", "GO:0003674", "GO:0005575"]
		self.label = ["i", "r", "p", "n", "o"]
		self.colors = ['blue', 'violet', 'green', 'red', 'orange']

def read_go_siglist(sig_gene_path):
	sig_golist = {}
	fp_go = file(sig_gene_path, 'r')
	for line in fp_go:
		if line.startswith("#") or line == '\n':
			continue
		else:
			goid = line.split("\t")[0].strip()
			sig_golist[goid] = [line.split("\t")[5].strip(), line.split("\t")[1].strip(), line.split("\t")[2].strip()]
	fp_go.close()
	return sig_golist
	
def select_go_siglist(go_siglist, qvalue = 0.05, maxnum = 0, qsite = 0):
	sig_go = {}
	if not maxnum:
		for iterm in go_siglist:
			if float(go_siglist[iterm][qsite]) <= qvalue:
				sig_go[iterm] = go_siglist[iterm]
	else:
		qvalue = []
		for iterm in go_siglist:
			if go_siglist[iterm][qsite]:
				qvalue.append(float(go_siglist[iterm][qsite]))
		qvalue_temp = [(v,i) for i,v in enumerate(qvalue)]
		qvalue_temp.sort()
		qvalue_sorted, index = zip(*qvalue_temp)
		cnt = 0
		for i in index:
			sig_go[go_siglist.keys()[i]] = go_siglist[go_siglist.keys()[i]]
			cnt += 1
			if cnt == maxnum:
				break
	return sig_go

def read_alt_rep(opt):
	alt_goid = {"current_id" : [] , "alt_id" : []}
	rep_goid = {"original_id" : [] , "replaced_id" : []}
	fp_alt = file(opt.altgo, 'r')
	fp_rep = file(opt.repgo, 'r')
	while True:
		fp_alt_line = fp_alt.readline()
		if not fp_alt_line:
			break
		else:
			if fp_alt_line.startswith('#'): continue
			current_id = fp_alt_line.split('\t')[0]
			alt_id = fp_alt_line.split('\t')[1].strip().split('|')
			for line in alt_id:
				alt_goid[alt_goid.keys()[0]].append(current_id)
				alt_goid[alt_goid.keys()[1]].append(line)
	del current_id, line
	while True:
		fp_rep_line = fp_rep.readline()
		if not fp_rep_line:
			break
		else:
			if fp_rep_line.startswith('#'): continue
			original_id = fp_rep_line.split('\t')[0]
			rep_id = fp_rep_line.split()[1].strip()
			rep_goid[rep_goid.keys()[0]].append(original_id)
			rep_goid[rep_goid.keys()[1]].append(rep_id)
	fp_alt.close(); fp_rep.close()
	return alt_goid, rep_goid

def check_alt_rep(sig_go, alt_goid, rep_goid):		
	tmp_sig_go = copy.deepcopy(sig_go)
	id_save = []
	for iterm in tmp_sig_go:
		if iterm in alt_goid[alt_goid.keys()[1]]:
			tmp_index = alt_goid[alt_goid.keys()[1]].index(iterm)
			tmp_id = alt_goid[alt_goid.keys()[0]][tmp_index]
			sys.stderr.write('[WARNING] %s was alternated by %s.\n'%(iterm, tmp_id))
			if tmp_id not in sig_go:
				sig_go[tmp_id] = tmp_sig_go[iterm]
				sig_go[tmp_id][1] = None
				sig_go[tmp_id][2] = None
			if iterm not in id_save:
				id_save.append(iterm)
		if iterm in rep_goid[rep_goid.keys()[0]]:
			tmp_index = rep_goid[rep_goid.keys()[0]].index(iterm)
			tmp_id = rep_goid[rep_goid.keys()[1]][tmp_index]
			rep_go = tmp_id.split('|')
			sys.stderr.write('[WARNING] %s was replaced by %s.\n'%(tmp_id, iterm))
			for rep_id in rep_go:
				if rep_id not in sig_go:
					sig_go[rep_id] = tmp_sig_go[iterm]
					sig_go[rep_id][1] = None
					sig_go[rep_id][2] = None
			if iterm not in id_save:
				id_save.append(iterm)
	if id_save:
		for iterm in id_save:
			del sig_go[iterm]
	return sig_go, id_save

def read_gofmt(fileName):
	sig_go = {}
	with open(fileName, 'r') as fp:
		for line in fp:
			if line.startswith('#'): continue
			else:
				go_attr = line.strip()
				if go_attr in sig_go:
					pass
				else:
					sig_go[go_attr] = ['0.0', '', '']
	return sig_go

def def_annotation(sig_go, defPath):
	with open(defPath, 'r') as fp:
		for line in fp:
			if line.startswith('#'): continue
			else:
				goid, des, go_class, other = line.split('\t', 3)
				if goid in sig_go:
					if sig_go[goid][1]:pass
					else:
						sig_go[goid][1] = des
						sig_go[goid][2] = go_class
				else:
					continue
	return sig_go

def read_golist(alt_goid, rep_goid, go_single_list):
	go_path = {}
	if go_single_list:
		if go_single_list.startswith('GO:'):
			go_path[go_single_list] = None
	elif os.path.isfile(go_single_list):
		fp_go = file(go_sing_list, 'r')
		go_read = fp_go.read()
		for iterm in go_read:
			go_path[iterm.split('\n')[0].strip()] = None
	
	tmp_go_path = copy.deepcopy(go_path)
	id_save = []
	for iterm in tmp_go_path:
		if iterm in alt_goid[alt_goid.keys()[1]]:
			tmp_idex = alt_goid[alt_goid.keys()[1]].index(iterm)
			tmp_id = alt_goid[alt_goid.keys()[0]][tmp_index]
			sys.stderr.write('[Warning] %s was alternated by %s.\n'%(iterm, tmp_id))
			if tmp_id not in go_path:
				go_path[tmp_id] = None
			if iterm not in id_save:
				id_save.append(iterm)
		if iterm in rep_goid[rep_goid.keys()[0]]:
			tmp_idex = rep_goid[rep_goid.keys()[0]].index(iterm)
			tmp_id = rep_goid[rep_goid.keys()[1]][tmp_index]
			rep_go = tmp_id.split('|')
			sys.stderr.write('[Warning] %s was replaced by %s.\n'%(tmp_id, iterm))
			for rep_id in rep_go:
				if tmp_id not in go_path:
					go_path[tmp_id] = go_path[iterm]
				if iterm not in id_save:
					id_save.append(iterm)
	if id_save:
		for iterm in id_save:
			del go_path[iterm]
	return go_path, id_save

def read_topo_flag(sig_go, id_save, opt, flag = 1):
	obo_path = {}; obo_flag = {};  go_add = []
	fp_path = file(opt.topology, 'r')
	path_str = fp_path.read().split("\n")
	fp_path.close()
	
	flag_id = 1; len_sig = len(sig_go)
	for iterm in path_str:
		if iterm.startswith("#"):
			continue
		else:
			tmp_iterm = iterm.split('\t')
			go_id = tmp_iterm[0]
			if not go_id:
				continue
			else:
				if go_id in sig_go:
					go_add = go_add + [go_id]
					if tmp_iterm[2]:
						obo_path[go_id] = copy.deepcopy([[], []])
						obo_flag[go_id] = copy.deepcopy([])
						obo_path[go_id][1] = tmp_iterm[1]
						if go_id in id_save:
							sig_go[go_id][flag] = tmp_iterm[1]
						for i in xrange(len(tmp_iterm[2])):
							obo_flag[go_id].append(tmp_iterm[2][i])
						tmp_iterm_path = tmp_iterm[3].split(";")
						for i in xrange(len(tmp_iterm_path)):
							obo_path[go_id][0].append(tuple(tmp_iterm_path[i].split("->")))
							go_add = list(set(go_add + tmp_iterm_path[i].split("->")))
					else:
						obo_path[go_id] = copy.deepcopy([[], []])
						obo_flag[go_id] = copy.deepcopy([])
						obo_path[go_id][1] = tmp_iterm[1]
					if flag_id == len_sig:
						break
					else:
						flag_id = flag_id + 1
	return obo_path, obo_flag, sig_go, go_add

def read_go_def(go_add, path):
	go_def = {}; flag_id = 1
	fp_def = file(path, 'r')
	def_read = fp_def.read()
	fp_def.close()
	def_read = def_read.split('\n')
	len_go_add = len(go_add)

	for iterm in def_read:
		if iterm.startswith('#'):
			continue
		else:
			iterm = iterm.split('\t')
			if iterm[0] in go_add:
				go_def[iterm[0]] = iterm[1]
				if flag_id == len_go_add:
					break
				else:
					flag_id = flag_id + 1
	return go_def

def color_pvalue(sig_go):
	qvalue = []; nodes_col = {}; flag_id = 0
	for go_id in sig_go:
		qvalue.append(float(sig_go[go_id][0]))
	np_min =  np.min(np.min(qvalue))
	np_max = np.max(np.max(qvalue))
	qvalue = ((qvalue - np_min) + 1)/(np_max - np_min + 1).tolist()
	for go_id in sig_go:
		nodes_col[go_id] = "#FFFF00%2x"%(int(qvalue[flag_id]*128)+64)
		flag_id = flag_id + 1
	return nodes_col

def network_pgv(sig_go, go_def, obo_path, obo_flag, nodes_colors, opt):
	if opt.divide == 0:
		G_pgv = AGraph(directed=True, rankdir = 'BT', )
		G_pgv.node_attr['style']='filled'
		G_pgv.edge_attr['style']='setlinewidth(2.5)'
		for go_id in sig_go:
			edges = obo_path[go_id][0]
			leng_edges = len(edges)
			for line in xrange(leng_edges):
				tmp = []
				for k in xrange(len(edges[line])):
					tmp_str = go_def[edges[line][k]]
					if len(tmp_str) <= 50:
						tmp_label = edges[line][k] + "\\n" +  tmp_str
					else:
						str_split = tmp_str.split(" ")
						if len(str_split) % 2 == 0:
							str_1 = " ".join(str_split[0 : int(float(len(str_split)) / 2)])
							str_2 = " ".join(str_split[int(float(len(str_split)) / 2) : len(str_split)])
							tmp_str = str_1 + "\\n" + str_2
						else:
							str_1 = " ".join(str_split[0 : int(math.ceil(float(len(str_split)) / 2))])
							str_2 = " ".join(str_split[int(math.ceil(float(len(str_split)) / 2)) : ])
							tmp_str = str_1 + "\\n" + str_2
					tmp_label = edges[line][k] + "\\n" +  tmp_str
					tmp.append(edges[line][k])
					if edges[line][k] in sig_go and edges[line][k] not in opt.go_class:
						tmp_color = nodes_colors[edges[line][k]]
						G_pgv.add_node(edges[line][k], shape = "rectangle", label = tmp_label, color = 'red')
						node_pgv = G_pgv.get_node(edges[line][k]); node_pgv.attr['fillcolor']= tmp_color 
					elif edges[line][k] in opt.go_class:
						G_pgv.add_node(edges[line][k], shape = "rectangle", label = tmp_label)
						node_pgv = G_pgv.get_node(edges[line][k]); node_pgv.attr['fillcolor']= 'gray'
					else:
						G_pgv.add_node(edges[line][k], shape = "rectangle",  label = tmp_label)
						node_pgv = G_pgv.get_node(edges[line][k]); node_pgv.attr['fillcolor']= 'white'
				for j in xrange(len(opt.label)):
					if obo_flag[go_id][line][0] in opt.label[j]:
						G_pgv.add_edge(tmp[0], tmp[1], color = opt.colors[j])
						edge_pgv = G_pgv.get_edge(tmp[0], tmp[1]); edge_pgv.attr['width'] = 2.0

		G_pgv.layout(prog = 'dot')
		G_pgv.draw(opt.output + 'GoTopNet_fdr_biological.process_celluar.component_molucular.function.png')
		G_pgv.draw(opt.output + 'GoTopNet_fdr_biological.process_celluar.component_molucular.function.svg')
	elif opt.divide > 0:
		G_pgv = {"biological_process": '', "molecular_function": '', "cellular_component": '' }
		for go_id in sig_go:
			for i in xrange(len(G_pgv.keys())):
				if sig_go[go_id][2] not in G_pgv.keys()[i]: continue
				if sig_go[go_id][2] in G_pgv.keys()[i]:
					if len(G_pgv[G_pgv.keys()[i]]) == 0:
						G_pgv[G_pgv.keys()[i]] = AGraph(directed=True, rankdir = 'BT')
						G_pgv[G_pgv.keys()[i]].node_attr['style']='filled'
						G_pgv[G_pgv.keys()[i]].edge_attr['style']='setlinewidth(2.5)'
					else:
						pass
					edges = obo_path[go_id][0]
					leng_edges = len(edges)
					for line in xrange(leng_edges):
						tmp = []
						for k in xrange(len(edges[line])):
							tmp_str = go_def[edges[line][k]]
							if len(tmp_str) <= 50:
								tmp_label = edges[line][k] + "\\n" +  tmp_str
							else:
								str_split = tmp_str.split(" ")
								if len(str_split) % 2 == 0:
									str_1 = " ".join(str_split[0 : int(float(len(str_split)) / 2)])
									str_2 = " ".join(str_split[int(float(len(str_split)) / 2) : len(str_split)])
									tmp_str = str_1 + "\\n" + str_2
								else:
									str_1 = " ".join(str_split[0 : int(math.ceil(float(len(str_split)) / 2))])
									str_2 = " ".join(str_split[int(math.ceil(float(len(str_split)) / 2)) : ])
									tmp_str = str_1 + "\\n" + str_2
							tmp_label = edges[line][k] + "\\n" +  tmp_str
							tmp.append(edges[line][k])
							if edges[line][k] not in opt.go_class:
								if edges[line][k] in sig_go:
									tmp_color = nodes_colors[edges[line][k]]
									G_pgv[G_pgv.keys()[i]].add_node(edges[line][k], shape = "rectangle", label = tmp_label)
									node_pgv = G_pgv[G_pgv.keys()[i]].get_node(edges[line][k]); node_pgv.attr['fillcolor']= tmp_color
								else:
									G_pgv[G_pgv.keys()[i]].add_node(edges[line][k], shape = "rectangle",  label = tmp_label)
									node_pgv = G_pgv[G_pgv.keys()[i]].get_node(edges[line][k]); node_pgv.attr['fillcolor']= 'white'
							else:
								G_pgv[G_pgv.keys()[i]].add_node(edges[line][k], shape = "rectangle", label = tmp_label)
								node_pgv = G_pgv[G_pgv.keys()[i]].get_node(edges[line][k]); node_pgv.attr['fillcolor']= 'gray'
						for j in xrange(len(opt.label)):
							if obo_flag[go_id][line][0] in opt.label[j]:
								G_pgv[G_pgv.keys()[i]].add_edges_from([(tmp[0], tmp[1])], color = opt.colors[j])
		for iterm in G_pgv.keys():
			try:
				G_pgv[iterm].layout(prog = 'dot')
				G_pgv[iterm].draw(opt.output + "/" + 'GoTopNet_fdr_' + iterm + '.png')
				G_pgv[iterm].draw(opt.output + "/" + 'GoTopNet_fdr_' + iterm + '.svg')
			except:
				sys.stderr.write('[WARNING] %s class has no nodes.\n'%iterm)
	
def __main():
	usage = "USAGE: %prog [options]"
	description = 'Contact: Li Huamei <li_hua_mei@163.com>'
	parser = OptionParser(usage, version = '%prog 1.0', description = description)
	OBB_group = OptionGroup(parser, 'Obligatory options')
	OBB_group.add_option('--run-stat', dest = 'stat', help = 'Starting run sigList program', action = 'store_true')
	OBB_group.add_option('--run-go', dest = 'run_go', help = 'Starting run GoList program', action = 'store_true')
	OBB_group.add_option('--gene-sig', dest = 'gene', help = 'Significance gene list file', metavar = 'STR', type = 'string')
	OBB_group.add_option('--go-list', dest = 'golist', help = 'GO list file', metavar = 'STR', type = 'string')
	OBB_group.add_option("--godef",dest="godef",help="GO items def file",metavar="FILE",type="string",default='/data/database/GO/20150811/GO_def.txt')
	OBB_group.add_option("--topology",dest="topology",help="Topology of semantic",metavar="FILE",type="string",default='/home/lihuamei/pipline/rnaseq/go_info/obo_topology.txt')
	OBB_group.add_option("--altgo",dest="altgo",help="Alternative of semantic",metavar="FILE",type="string",default='/home/lihuamei/pipline/rnaseq/go_info/GO_alt_id.alt_id')
	OBB_group.add_option("--repgo",dest="repgo",help="Replaced by semantic",metavar="FILE",type="string",default='/home/lihuamei/pipline/rnaseq/go_info/GO_replace.replace')

	OBB_group.add_option('--divide', dest = 'divide', help = 'Seprate the significance semantic network as bp, mf and cc or not, combineed bp, cc and mf process in a total network is default', metavar = 'INT', type = 'int', default = 0)
	OBB_group.add_option('--qvalue', dest = 'qvalue', help = 'Set the cutoff value of semantics', metavar = 'FLOAT', type = 'float', default = 0.05)
	OBB_group.add_option('--output', dest = 'output', help = 'Specify the output directory', metavar = 'DIR', type = 'string', default = './')

	Optional_group = OptionGroup(parser, 'Optional commands')
	Optional_group.add_option('--gpath', dest = 'path', help = 'Generate specified semantic paths', metavar = 'FILE or Single GO ID', type = 'string')
	Optional_group.add_option('--gpnet', dest = 'network', help = 'Generate specified semantic network', metavar = 'FILE or Single GO ID', type = 'string')
	Optional_group.add_option('--spath', dest = 'shortest', help = 'Search out the shorest path, semantic id as input', metavar = 'FILE or Single GO ID', type = 'string')
	Optional_group.add_option('--lpath', dest = 'longest', help = 'Search out the longest path, semantic id as input', metavar = 'FILE or Single GO ID', type = 'string')
	Optional_group.add_option('--maxnum', dest = 'max', help = 'Utilize the number of top significant semantic', metavar = 'INT', type = 'int')
	Optional_group.add_option('--piles', dest = 'piles', help = 'Product the number of piles for specifying semantic', metavar = 'FILE or Single GO ID', type = 'string')
	Optional_group.add_option('--pilenum', dest = 'pilenum', help = 'Specifying the number of piles, default = 4', metavar = 'INT', type = 'int')
	
	parser.add_option_group(OBB_group)
	parser.add_option_group(Optional_group)

	(options, argv) =  parser.parse_args()
	opt = obo_file_path()
	opt.gene = options.gene
	opt.divide = options.divide
	opt.qvalue = options.qvalue
	opt.maxnum = options.max
	opt.godef = options.godef
	opt.topology = options.topology
	opt.altgo = options.altgo
	opt.repgo = options.repgo
	opt.golist = options.golist
	opt.output = options.output
	if options.stat:
		opt.run = 1
		if os.path.isfile(opt.godef) and os.path.isfile(opt.altgo) and os.path.isfile(opt.repgo):
			pass
		else:
			parser.print_help()
			exit(1)
	else:
		opt.run = 0
		if not options.run_go:
			parser.print_help()
			exit(1)
		else:
			if not os.path.isfile(options.golist):
				parser.print_help()
				exit(1)

	return opt
def call_function_sig(opt):
	alt_goid, rep_goid = read_alt_rep(opt)
	sig_golist = read_go_siglist(opt.gene)
	sig_go = select_go_siglist(sig_golist, qvalue = opt.qvalue, maxnum = opt.maxnum, qsite = 0)
	if len(sig_go) == 0:
		sys.stderr.write('[WARNING] NO Significance Semantic.\n')
		exit(0)
	sig_go, id_save = check_alt_rep(sig_go, alt_goid, rep_goid)
	obo_path, obo_flag, sig_go, go_add = read_topo_flag(sig_go, id_save, opt, flag = 1)
	go_def = read_go_def(go_add, opt.godef)
	nodes_colors = color_pvalue(sig_go)
	network_pgv(sig_go, go_def, obo_path, obo_flag, nodes_colors, opt)

def call_function_go(opt):
	alt_goid, rep_goid = read_alt_rep(opt)
	sig_go = read_gofmt(opt.golist)
	if len(sig_go) == 0:
		sys.stderr.write('[WARNING] NO Significance Semantic.\n')
		exit(0)
	sig_go = def_annotation(sig_go, opt.godef)
	sig_go, id_save = check_alt_rep(sig_go, alt_goid, rep_goid)
	obo_path, obo_flag, sig_go, go_add = read_topo_flag(sig_go, id_save, opt, flag = 1)
	go_def = read_go_def(go_add, opt.godef)
	nodes_colors = color_pvalue(sig_go)
	network_pgv(sig_go, go_def, obo_path, obo_flag, nodes_colors, opt)


if __name__ == '__main__':
	Packages = ['from __future__ import division', 'import sys',\
			'import os', 'import time', 'import copy', 'import string',\
			'import pdb','import networkx as nx', 'import numpy as np',\
			'import matplotlib.pyplot as plt', 'import pygraphviz as pgv',\
			'from pygraphviz import *', 'import pydot', 'import math',\
			'import scipy', 'from optparse import OptionParser, OptionGroup']
	print ('[INFO] Test Modules.')
	for module in Packages:
		try:
			exec(module)
		except:
			 print("[Error] No module named %s in the system, program exit.\n"%(module.split()[1]))
			 exit(1)
	del module, Packages
	options = __main()
	start_time = time.time()
	if options.run:
		call_function_sig(options)
	else:
		call_function_go(options)
	finish_time = time.time() - start_time
	sys.stderr.write('[CONSUME TIME] Consume %s seconds.\n'%str(finish_time))
	sys.stderr.write('[INFO] Task Done.\n')

