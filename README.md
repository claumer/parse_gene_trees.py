Parse_gene_trees.py v1.0

Requires the Environment for Tree Exploration v2 (ETE2), argparse 1.1+.

The purpose of this script is to parse a large set of gene trees with missing taxa and,
for each rooted node in a single user-defined input tree, annotate the proportion of genes
in the large set that are a.) potentially decisive for that node, that is, containing at 
least one member of each descendant branch plus two distinct outgroups, and 
b.) the number within the above set that are actually congruent with that node. 

Usage is as follows (by default prints to stdout):

python parse_gene_trees.py -i [newick-formatted input gene tree containing all taxa, whose nodes are to be annotated] \
-t [large set of gene trees with missing taxa, newick-formatted, one per line] \ 
-o [String representing the name of some outgroup with which to root the input gene tree]
