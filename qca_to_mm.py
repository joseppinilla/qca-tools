import os
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt

from parse_qca import parse_qca_file

from embedding_methods.utilities.graph_mmio import write_networkx

class QCANetwork(nx.Graph):

    def __init__(self, qca_file=None, full_adj=True):

        self.full_adj = full_adj

        # read qca file
        cells, spacing, J_mtx = parse_qca_file(qca_file)
        self.cells = cells
        self.spacing = spacing

        nx.Graph.__init__(self, J_mtx)

        self.pos = {}
        for cell in self.cells:
            num = cell['num']
            self.pos[num] = cell['x'], cell['y']

        # create graph
        n, m = J_mtx.shape

        # qca adjacency
        self.qca_adj = {i: J_mtx[i].nonzero()[0].tolist() for i in range(J_mtx.shape[0])}

        # cell functions
        self.fixed = set()
        self.normal = set()
        self.outputs = set()
        self.drivers = set()
        for cell in self.cells:
            num = cell['num']
            func = cell['cf']
            if func == 'QCAD_CELL_FIXED':
                self.fixed.add(num)
            elif func == 'QCAD_CELL_NORMAL':
                self.normal.add(num)
            elif func == 'QCAD_CELL_INPUT':
                self.drivers.add(num)
            elif func == 'QCAD_CELL_OUTPUT':
                self.outputs.add(num)
            else:
                raise ValueError('Cell %s function %s is not valid.' % (cell,func))

        # set up input polarization vector
        P = np.asmatrix(np.zeros([n, 1], dtype=float))

        # fixed cell contributions
        for num in self.fixed:
            P[num] = cells[num]['pol']


        # driver cell contributions
        for num in self.drivers:
            cname = cells[num]['name']
            pol = input('Input polarization for %s(%s):' % (cname,num))
            P[num] = pol

        # get h coefficients
        h_vect = np.round(np.asmatrix(J_mtx)*P, 8)

        active_cells, reduced_qca_adj = self.get_reduced_qca_adj()

        h = {c:float(h_vect[c]) for c in active_cells}
        self.h = h
        for x, val in h.items():
            self.add_edge(x,x,weight=val)

        J = {}
        for u,v in self.edges():
            if u in reduced_qca_adj and v in reduced_qca_adj:
                J[(u,v)] = J_mtx[u,v]
        self.J = J

    def qca_layout(self):
        return self.pos


    def get_reduced_qca_adj(self):
        '''Get a reduced form of qca_adj only for non-driver/fixed cells'''

        # check function for membership in drivers or fixed
        check = lambda cell: cell not in self.drivers.union(self.fixed)

        reduced_adj = {c1: [c2 for c2 in self.qca_adj[c1] if check(c2)]
            for c1 in self.qca_adj if check(c1)}

        return sorted(reduced_adj), reduced_adj

if __name__ == "__main__":

    bench_dir = './benchmarks/'
    fn = 'SRFlipFlop.qca'
    base, ext = os.path.splitext(fn)
    if not os.path.exists(bench_dir):
        os.makedirs(bench_dir)

    G = QCANetwork(qca_file=bench_dir+fn)

    pos = G.qca_layout()

    nx.draw(G, pos=pos, with_labels=True)
    plt.gca().invert_yaxis()
    plt.show()

    write_networkx(G, pos=pos, mtx_name=base+'_G', mm_dir=bench_dir)
