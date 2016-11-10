
def suffix_array(suffix_tree):
	root = suffix_tree.get_root_node()
	cur_node = root

	print root
	for edge_char, edge_start in cur_node.get_edges().iteritems():
	    edge = suffix_tree.get_tree_edge( cur_node, edge_char )
	    edge_idx = cur_node.get_edge( edge_char )
	    
	    print "edge char:", edge_char, "(", edge_start, ")"
	    print "\tchild start:", edge.get_start()
	    print "\tchild end:", edge.get_end()
	    print "\tchild substring:", suffix_tree.edge_string( edge_idx )


	edge_chars = [i[0] for i in root.get_edges().items()]
	edge_chars.sort()

	print "-------------------------------------"
	print root.get_edges()
	print edge_chars