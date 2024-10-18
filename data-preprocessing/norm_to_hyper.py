normal_name = "graph.txt"
bipartite_name = "bipartite_hyper.txt"

# Use a set to track unique edges and avoid duplicates
unique_edges = set()

with open(normal_name, 'r') as normal_graph, open(bipartite_name, 'w') as bi_graph:
    idx = 0
    for line in normal_graph:
        # Split the line by tab and strip any trailing newlines
        nodes = line.strip().split(" ")
        #nodes[1] = nodes[1][0:-1]
        nodes.sort()
        # Ensure we have exactly two nodes
        if len(nodes) != 2:
            continue
        if (nodes[0] == nodes[1]):
            continue
        
        # Use frozenset to handle unordered pairs (an edge is the same regardless of direction)
        edge = frozenset(nodes)

        # Avoid self-loops (where both nodes are the same)
        if len(edge) == 1:
            continue
        
        # Add the edge to the set if it's not already present
        if edge not in unique_edges:
            unique_edges.add(edge)
            
            # Write the bipartite edges
            node_list = list(edge)
            bi_graph.write(f"{node_list[0]} {node_list[1]}\n")
            bi_graph.write(f"{node_list[1]} {node_list[0]}\n")
            idx += 1

print("converted")
