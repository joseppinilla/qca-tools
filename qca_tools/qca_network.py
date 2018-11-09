import os
import numpy as np
import networkx as nx

import matplotlib.pyplot as plt

from qca_tools.parse_qca import parse_qca_file

from embedding_methods.utilities.graph_mmio import write_networkx
from dimod import BinaryQuadraticModel, Vartype

ROUND_VAL = 8

class QCANetwork(BinaryQuadraticModel):

    def __init__(self, qca_file=None, full_adj=True, pols={}, ancilla=True):

        # init empty BQM
        BinaryQuadraticModel.__init__(self,{},{},0.0,Vartype.SPIN)

        # properties of QCANetwork
        self.full_adj = full_adj
        self.pols = pols

        # properties from file
        self.pos = {}
        self.cells = {}
        self.spacing = None
        self.qca_adj = {}
        self.fixed = set()
        self.normal = set()
        self.outputs = set()
        self.drivers = set()
        self.active_cells = set()

        # empty network
        if qca_file is None: return

        # read qca file
        cells, spacing, J_mtx = parse_qca_file(qca_file)

        # qca adjacency
        qca_adj = {i: J_mtx[i].nonzero()[0].tolist() for i in range(J_mtx.shape[0])}

        # cell functions
        for cell in cells:
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

        # driver cell polarizations
        for num in self.drivers:
            cname = cells[num]['name']
            if cname not in pols:
                pols[cname] = float(input('Input polarization for %s(%s):' % (cname,num)))

        n, m = J_mtx.shape
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

            # only for active cells
            self.active_cells.update(self.normal, self.outputs)
            # get h coefficients
            h_vect = np.round(np.asmatrix(J_mtx)*P, ROUND_VAL)
            for v in self.active_cells:
                bias = h_vect[v]
                self.add_variable(v, bias)

            # get J coefficients.
            J = {}
            for u in self.active_cells:
                for v in qca_adj[u]:
                    if (v in self.active_cells) and (u < v):
                        bias = np.round(J_mtx[u,v], ROUND_VAL)
                        self.add_interaction(u,v,bias)

        else:
            # set values of ancillas in matrix
            for num in self.fixed:
                pol = cells[num]['pol']
                J_mtx[num,num] = -pol
            for num in self.drivers:
                cname = cells[num]['name']
                pol = pols[cname]
                J_mtx[num,num] = -pol

            # fill BQM
            variable_order = list(range(n))
            for (row, col), bias in np.ndenumerate(J_mtx):
                if row == col:
                    self.add_variable(variable_order[row], bias)
                elif bias:
                    self.add_interaction(variable_order[row], variable_order[col], bias)

            self.active_cells.update(self.fixed, self.normal, self.outputs, self.drivers)

        # cell positions
        for name in self.active_cells:
            cell = cells[name]
            self.pos[name] = cell['x'], cell['y']

        self.spacing = spacing
        self.qca_adj = qca_adj
        self.cells = cells

    def qca_layout(self):
        return self.pos

class QCANetworkX(nx.Graph):
    def __init__(self, qca_file=None, full_adj=True, pols={}, ancilla=True):
        QCA = QCANetwork(qca_file, full_adj, pols, ancilla)
        self.QCA =  QCA
        self.pos = QCA.qca_layout()

        nx.Graph.__init__(self)
        for v, val in QCA.linear.items():
            self.add_edge(v, v, weight=val)
        for edge, val in QCA.quadratic.items():
            self.add_edge(*edge, weight=val)


    def draw_qca(self, with_biases=False, with_weights=False):

        QCA = self.QCA
        pos = self.pos

        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("Matplotlib required for draw_qca()")
        except RuntimeError:
            print("Matplotlib unable to open display")
            raise

        cf = plt.gcf()
        ax = cf.gca()

        node_color = []
        for node in self.nodes:
            if node in QCA.drivers:
                node_color.append('blue')
            elif node in QCA.outputs:
                node_color.append('yellow')
            elif node in QCA.fixed:
                node_color.append('orange')
            elif node in QCA.normal:
                node_color.append('green')

        weights = nx.get_edge_attributes(self,'weight')
        edge_color = list(weights.values())
        edge_vmin = min(QCA.quadratic.values())
        edge_vmax = max(QCA.quadratic.values())

        params = {  'pos': pos,
                    'node_shape': '$\u2683$', #TODO: if node is rotated?
                    'with_labels': True,
                    'node_color': node_color,
                    'node_size': QCA.spacing*100,
                    'width': 4,
                    'edge_color': edge_color,
                    'edge_cmap': plt.cm.RdBu,
                    'edge_vmin': edge_vmin,
                    'edge_vmax': edge_vmax
                     }

        if with_biases:
            params['labels'] = {v:weights[v,v] for v in self}

        nx.draw_networkx(self, **params)


        if with_weights:
            edge_labels = {(u,v):weights[u,v] for u,v in weights if u!=v}
            _ = nx.draw_networkx_edge_labels(self, pos=pos, edge_labels=edge_labels)

        ax.set_axis_off()
        ax.invert_yaxis()

    @classmethod
    def from_qca_network(cls, QCA):

        self.QCA = QCA


if __name__ == "__main__":
    QCA = QCANetworkX('../examples/benchmarks/NOT_FT.qca', ancilla=True)
    QCA.draw_qca()
    plt.show()
