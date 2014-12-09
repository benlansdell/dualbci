#
# Python Script to generate complete graph
#

NUM_NODES = 20

print "Creating complete graph with " + str(NUM_NODES) + " nodes..."

new_network = cyAppAdapter.getCyNetworkFactory().createNetwork()
new_network.getRow(new_network).set("name", "Complete Graph Created by Python Script")
cyAppAdapter.getCyNetworkManager().addNetwork(new_network);

# Add nodes
nodes = [];
for i in range (NUM_NODES):
    node_name = "Node " + str(i)
    node = new_network.addNode()
    new_network.getRow(node).set("name", node_name)
    nodes.append(node)

# Add edges
edge_count = 0;
for source in nodes:
    for target in nodes:
        if new_network.containsEdge(source, target) is False \
                and new_network.containsEdge(target, source) is False \
                and source is not target:
            edge = new_network.addEdge(source, target, True)
            edge_count = edge_count + 1
            new_network.getRow(edge).set("name", "Edge " + str(edge_count))
            new_network.getRow(edge).set("interaction", "interacts_with")