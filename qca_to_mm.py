import os
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt

from parse_qca import parse_qca_file

from embedding_methods.utilities.graph_mmio import write_networkx


ROUND_VAL = 8

#class QCANetworkX(QCANetwork, nx.Graph):

class QCANetwork(nx.Graph):

    def __init__(self, qca_file=None, full_adj=True, pols={}, ancilla=True):

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

        # driver cell polarizations
        for num in self.drivers:
            cname = cells[num]['name']
            if cname not in pols:
                pols[cname] = float(input('Input polarization for %s(%s):' % (cname,num)))

        if not ancilla:
            # set up input polarization vector
            P = np.asmatrix(np.zeros([n, 1], dtype=float))

            # fixed cell contributions
            for num in self.fixed:
                P[num] = cells[num]['pol']

            # driver cell contributions
            for num in self.drivers:
                cname = cells[num]['name']
                P[num] = pols[cname]

            active_cells = self.normal.union(self.outputs)

            # Only for active cells
            # get h coefficients
            h_vect = np.round(np.asmatrix(J_mtx)*P, ROUND_VAL)
            h = {c:float(h_vect[c]) for c in active_cells}

            # get J coefficients.
            J = {}
            for u in active_cells:
                for v in self.qca_adj[u]:
                    if (v in active_cells) and (u < v):
                        J[(u,v)] = np.round(J_mtx[u,v], ROUND_VAL)

            # Create Graph from Ising model
            nx.Graph.__init__(self)
            for v, val in h.items():
                self.add_edge(v, v, weight=val)
            for edge, val in J.items():
                self.add_edge(*edge, weight=val)

        else:
            # set values of ancillas in matrix
            for num in self.fixed:
                pol = cells[num]['pol']
                J_mtx[num,num] = -pol
            for num in self.drivers:
                cname = cells[num]['name']
                pol = pols[cname]
                J_mtx[num,num] = -pol

            # Create Graph
            nx.Graph.__init__(self, J_mtx)
            # Create Ising model from Graph
            h = {u:data['weight'] for u,v,data in self.edges(data=True) if u==v}
            J = {(u,v):data['weight'] for u,v,data in self.edges(data=True) if u!=v}

            active_cells = set(self.nodes())

        self.h = h
        self.J = J
        self.active_cells = active_cells

        # cell positions
        self.pos = {}
        for num in active_cells:
            cell = cells[num]
            self.pos[num] = cell['x'], cell['y']

    def qca_layout(self):
        return self.pos
