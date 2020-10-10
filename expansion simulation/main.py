## real stimulation
# we still use rectangular grid
# a single cell in a "cell"

import networkx as nx

m = 5  # number of rows of hexagons in the lattice
n = 5  # number of columns of hexagons in the lattice
isPBC = False # if True, use periodic boundary conditions
#
### build graph using networkx
tri_graph = nx.generators.lattice.triangular_lattice_graph(m, n, periodic=isPBC)
# label graph nodes by consecutive integers
tri_graph = nx.convert_node_labels_to_integers(tri_graph)
# set number of lattice sites
N = tri_graph.number_of_nodes()
print('constructed triagonal lattice with {0:d} sites.\n'.format(N))
# visualise graph
pos = nx.spring_layout(tri_graph, seed=42, iterations=100)
nx.draw(tri_graph, pos=pos, with_labels=True)
plt.show(block=False)

#%%
import turtle as t
t.setup(width=600, height=600)
t.left(90)

for time in range(3):
    t.fillcolor("black")
    t.begin_fill()
    for i in range(6):
        t.forward(10)
        t.right(60)
    t.end_fill()
    t.forward(30)

ts = t.getscreen()

#%%
deltas = [[1, -1], [0, 1, [-1, 0], 1], -1, [1, 0]]


class HexGrid():
    def __init__(self, radius):
        self.radius = radius
        self.tiles = {(0, 0): "X"}
        for r in range(radius):
            a = 0
            b = -r
            c = +r
            for j in range(6):
                num_of_hexas_in_edge = r
                for i in range(num_of_hexas_in_edge):
                    a = a + deltas[j][0]
                    b = b + deltas[j][1]
                    c = c + deltas[j][2]
                    self.tiles[a, b, c] = "X"

    def show(self):
        l = []
        for y in range(20):
            l.append([])
            for x in range(60):
                l[y].append(".")
        for (a, c), tile in self.tiles.iteritems():
            l[self.radius - 1 - b][a - c + (2 * (self.radius - 1))] = self.tiles[a, c]
        mapString = ""
        for y in range(len(l)):
            for x in range(len(l[y])):
                mapString += l[y][x]
            mapString += "\n"
        print(mapString)

