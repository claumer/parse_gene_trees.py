# parse_gene_trees.py v1.0
#
# Requires the Environment for Tree Exploration v2 (ETE2), argparse 1.1+.
#
# The purpose of this script is to parse a large set of gene trees with missing taxa and,
# for each rooted node in a single user-defined input tree, annotate the proportion of genes
# in the large set that are a.) potentially decisive for that node, that is, containing at 
# least one member of each descendant branch plus two distinct outgroups, and 
# b.) the number within the above set that are actually congruent with that node. 
#
# Usage is as follows (by default prints to stdout):
# python parse_gene_trees.py -i [newick-formatted input gene tree containing all taxa, whose nodes are to be annotated] \
# -t [large set of gene trees with missing taxa, newick-formatted, one per line] \ 
# -o [String representing the name of some outgroup with which to root the input gene tree]

from ete2 import Tree
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i')
parser.add_argument('-t')
parser.add_argument('-o')
args = parser.parse_args()

def load_PGTs():
	'''
	Opens a handle and uses the ETE2 tree function to create a list of 
	tree objects, and returns this list.
	'''
	newick_handle = open(args.t, 'r')
	trees = []
	for line in newick_handle:
		t = Tree(line.rstrip('\n'))
		trees.append(t)
	return trees
	
def nodes_to_test():
	'''
	Opens up some individual tree whose nodes we want to annotate with GSFs,
	and returns a list of the nodes in that tree.
	'''
	tree = Tree(args.i)
	nodelist = []
	for node in tree.traverse('preorder'):
		if len(node.get_leaf_names()) > 1: 
			nodelist.append(node)
	return nodelist
			
def find_splits(tree, node):
	'''
	For an input node on the tree-of-interest, returns a three-item list: 
	the taxa on one side of a split, those on the other side of it, 
	and the taxa outside the node. Each item is itself a set.
	'''
	one_side = set(node.get_descendants()[0].get_leaf_names())
	other_side = set(node.get_descendants()[0].get_sisters()[0].get_leaf_names())
	#other_side = set(node.get_descendants()[1].get_leaf_names())
	outgroups = []
	taxa = [taxon.get_leaf_names()[0] for taxon in tree]
	for taxon in taxa:
		if taxon not in one_side and taxon not in other_side:
			outgroups.append(taxon)
	outgroups = set(outgroups)
	assert set(taxa) == one_side | other_side | outgroups
	return (one_side, other_side, outgroups)
	
def test_to_include(tree, node):
	'''
	Given an input partial gene tree and an input node, tests if that tree contains at least
	one taxon on both sides of the split defined by that node, and one taxon outside the node
	of interest. Trees that meet this definition are able to non-trivially support the
	node of interest.
	'''
	tree_of_interest = Tree(args.i)
	tree_of_interest.set_outgroup(args.o)
	splits = find_splits(tree_of_interest,node)	
	for split in splits[:2]: #this for statement checks whether there is at least one taxon in both branches past the node of interest
		missing = 0
		for taxon in split:
			if taxon not in tree:
				missing += 1
		if missing == len(split):
			return False
	missing = 0
	for taxon in splits[2]: # this for statement checks whether there are at least two taxa outside the node of interest present in the PGT
		if taxon not in tree:
			missing += 1
	if len(splits[2]) - missing < 2:
		return False
	return True

def node_dict():
	'''
	Returns a dictionary with keys consisting of all nodes in the input tree, and
	values consisting of all trees that could possibly support that node (that is,
	all trees with one or more taxa on either side of the split).
	'''	
	node_dict = {}
	nodes = nodes_to_test()
	trees = load_PGTs()
	for node in nodes:
		node_dict[node] = []
		for tree in trees:
			if test_to_include(tree,node):
				node_dict[node].append(tree)
	return node_dict
	
def is_tree_congruent(tree, splits):
	taxa = [taxon.get_leaf_names()[0] for taxon in tree]
	one_side, other_side = [], []
	for taxon in splits[0]:
		if taxon in taxa:
			one_side.append(taxon)
	for taxon in splits[1]:
		if taxon in taxa:
			other_side.append(taxon)
	one_side, other_side = set(one_side), set(other_side)
	both = one_side | other_side
	def get_mrca(taxgroup):
		if len(taxgroup) == 1:
			return set(taxgroup)
		else:
			return set(tree.get_common_ancestor(taxgroup).get_leaf_names())
	if both != get_mrca(both):
		return False
	return True

node_dict = node_dict()
tree_of_interest = Tree(args.i)
tree_of_interest.set_outgroup(args.o)
for node in node_dict:
	trees = node_dict[node]
	print node
	splits = find_splits(tree_of_interest, node)
	count = 0
	for tree in trees:
		if is_tree_congruent(tree, splits):
			count += 1
	print count, len(trees)
