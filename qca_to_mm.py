%matplotlib tk

import os
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt

from parse_qca import parse_qca_file

from embedding_methods.utilities.graph_mmio import write_networkx


ROUND_VAL = 8

#class QCANetworkX(QCANetwork, nx.Graph):

class QCANetwork(nx.Graph):

    def __init__(self, qca_file=None, full_adj=True, pols={}):

        # properties of network
        self.full_adj = full_adj
        self.pols = pols

        # read qca file
        cells, spacing, J_mtx = parse_qca_file(qca_file)
        self.cells = cells
        self.spacing = spacing

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

        n, m = J_mtx.shape

        # set up input polarization vector
        P = np.asmatrix(np.zeros([n, 1], dtype=float))

        # fixed cell contributions
        for num in self.fixed:
            P[num] = cells[num]['pol']

        # driver cell contributions
        for num in self.drivers:
            cname = cells[num]['name']
            pol = pols.setdefault(cname, input('Input polarization for %s(%s):' % (cname,num)))
            P[num] = pol

        active_cells = self.normal.union(self.outputs)
        self.active_cells = active_cells

        # get h coefficients
        h_vect = np.round(np.asmatrix(J_mtx)*P, ROUND_VAL)
        h = {c:float(h_vect[c]) for c in active_cells}
        self.h = h

        # get J coefficients.
        J = {}
        for u in active_cells:
            for v in self.qca_adj[u]:
                if (v in active_cells) and (u < v):
                    J[(u,v)] = np.round(J_mtx[u,v], ROUND_VAL)
        self.J = J

        # Create Graph from Ising model
        nx.Graph.__init__(self)
        for v, val in h.items():
            self.add_edge(v, v, weight=val)
        for edge, val in J.items():
            self.add_edge(*edge, weight=val)

        # Store positions with driver and fixed cells removed
        self.pos = {}
        for num in active_cells:
            cell = cells[num]
            self.pos[num] = cell['x'], cell['y']

    def qca_layout(self):
        return self.pos


if __name__ == "__main__":

    bench_dir = './benchmarks/'
    #fn = 'SRFlipFlop.qca'
    #fn = 'mux2to1.qca'
    #fn = 'NOT.qca'
    #fn = 'AND4.qca'
    fn = 'half_adder.qca'

    base, ext = os.path.splitext(fn)
    if not os.path.exists(bench_dir):
        os.makedirs(bench_dir)

    G = QCANetwork(qca_file=bench_dir+fn)

    pos = G.qca_layout()

    nx.draw(G, pos=pos, labels=G.h)
    _ = nx.draw_networkx_edge_labels(G,pos=pos, edge_labels=G.J)
    plt.gca().invert_yaxis()

    pols = G.pols

    comments = "Source: %s\nPolarizations: %s" % (fn, pols)

    write_networkx(G, pos=pos, mtx_name=base, mm_dir=bench_dir, comment=comments)
