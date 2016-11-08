from suffix_tree import SuffixTree
from suffix_array import suffix_array


text = "panamabanana$"
suffix_tree = SuffixTree( text )

suffix_tree.print_edges(1)
#suffix_tree.tree_to_file("st")
#suffix_array(suffix_tree)