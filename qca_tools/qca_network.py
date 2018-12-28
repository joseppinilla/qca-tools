import os
import numpy as np
import networkx as nx
import matplotlib.markers as markers
import matplotlib.pyplot as plt

from qca_tools.parse_qca import parse_qca_file

from embedding_methods.utilities.graph_mmio import write_networkx
from dimod import BinaryQuadraticModel, Vartype

ROUND_VAL = 8
R_MAX = 2.1

__all__ = ['QCANetwork', 'QCANetworkX']

class QCANetwork(BinaryQuadraticModel):
    """ Load from a *.qca file and generate the adjacency matrix.
    """

    def __init__(self, qca_file=None,
                full_adj=True,
                pols={},
                ancilla=True,
                r_max=R_MAX):

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
        cells, spacing, J_mtx = parse_qca_file(qca_file, r_max=r_max)

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
    """ QCANetwork to NetworkX
    """

    def __init__(self, QCA):
        self.QCA =  QCA
        self.pos = QCA.qca_layout()

        nx.Graph.__init__(self)
        for v, val in QCA.linear.items():
            self.add_edge(v, v, weight=val)
        for edge, val in QCA.quadratic.items():
            self.add_edge(*edge, weight=val)

    def draw_qca(self, with_cmap=False, with_labels=False,
                with_biases=False, with_weights=False,
                scale=2.0, edge_width=2.0):
        """ Draw the QCA layout with heatcolored edges between adjacent cells

        Note: node_size is constructed from the cell's size and scale
        Default scale: nanometers
        Default cell size: 18x18(nm)

        Args (optional):

            with_cmap (bool, default=False): Add colormap legend for edge weights.

            with_biases (bool, default=False): Overlay bias numbers on top of
                nodes. Not recommended for good visibility.

            with_weights (bool, default=False): Overlay weight numbers on top of
                edges. Not recommended for good visibility.

            scale (float, default=2.0): A value of 2 is used to scale the layout
                size so that cells are distinguishable in a *.png exported file.

            edge_width (float, default=2.0): Configurable width of edges between
                adjacent cells.

        """

        # Plot setup
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            raise ImportError("Matplotlib required for draw_qca()")
        except RuntimeError:
            print("Matplotlib unable to open display")
            raise
        fig, ax = plt.subplots()
        ax.autoscale()

        # Dimensions setup
        QCA = self.QCA
        pos = {k:(scale*x, scale*y) for k,(x,y) in QCA.pos.items()}
        spacing = scale*QCA.spacing
        node_size = (spacing*72./fig.dpi)**2
        node_shape = markers.MarkerStyle(marker='$\u2683$')

        # Set drawing parameters
        deg0_cells = [cell['num'] for cell in QCA.cells if not cell['rot']]
        deg45_cells = [cell['num'] for cell in QCA.cells if cell['rot']]
        weights = nx.get_edge_attributes(self,'weight')

        deg0_color = []
        for cell in deg0_cells:
            if cell in QCA.drivers:
                deg0_color.append('blue')
            elif cell in QCA.outputs:
                deg0_color.append('peru')
            elif cell in QCA.fixed:
                deg0_color.append('orange')
            elif cell in QCA.normal:
                deg0_color.append('green')


        deg45_color = []
        for cell in deg45_cells:
            if cell in QCA.drivers:
                deg45_color.append('blue')
            elif cell in QCA.outputs:
                deg45_color.append('peru')
            elif cell in QCA.fixed:
                deg45_color.append('orange')
            elif cell in QCA.normal:
                deg45_color.append('green')

        # Draw nodes
        node_params = { 'pos': pos,
                        'with_labels': with_labels,
                        'node_size': node_size}

        nx.draw_networkx_nodes(self, nodelist=deg0_cells,
                                node_color=deg0_color,
                                node_shape=node_shape,
                                **node_params)
        node_shape._transform = node_shape.get_transform().rotate_deg(45)
        nx.draw_networkx_nodes(self, nodelist=deg45_cells,
                                node_color=deg45_color,
                                node_shape=node_shape,
                                **node_params)

        # Node labels
        if with_labels:
            if with_biases:
                node_labels = {v:weights[v,v] for v in self}
            else:
                node_labels = {v:v for v in self}
            nx.draw_networkx_labels(self, pos=pos, labels=node_labels)

        # Use scaled dim of drawing to set figure size
        xmin,xmax = plt.xlim()
        ymin,ymax = plt.ylim()
        xdim = (xmax-xmin+spacing)/fig.dpi
        ydim = (ymax-ymin+spacing)/fig.dpi
        fig.set_size_inches(xdim, ydim)

        # Draw edges
        edge_color = list(weights.values())
        edge_vmin = min(QCA.quadratic.values())
        edge_vmax = max(QCA.quadratic.values())
        edge_params = { 'pos': pos,
                        'width': edge_width,
                        'edge_color': edge_color,
                        'edge_cmap': plt.cm.RdBu,
                        'edge_vmin': edge_vmin,
                        'edge_vmax': edge_vmax,
                        'alpha': 0.5
                        }
        lc = nx.draw_networkx_edges(self, **edge_params)

        # Use NetworkX LineCollection to generate a colormap legend
        if with_cmap:
            plt.colorbar(lc)

        # Weigh labels
        if with_weights:
            edge_labels = {(u,v):weights[u,v] for u,v in weights if u!=v}
            _ = nx.draw_networkx_edge_labels(self,
                                        pos=pos, edge_labels=edge_labels)

        ax.set_axis_off()
        ax.invert_yaxis()

    @classmethod
    def from_qca_file(cls, qca_file=None, full_adj=True, pols={}, ancilla=True,
                    r_max=R_MAX):

        QCA = QCANetwork(qca_file, full_adj, pols, ancilla, r_max)
        return cls(QCA)

if __name__ == "__main__":

    dir= '../examples/benchmarks/'

    benchmarks = []
    benchmarks.append('4BITACCUM')
    benchmarks.append('4BITMUX')
    benchmarks.append('SERADD')
    benchmarks.append('LOOPMEMCELL')
    benchmarks.append('FULLADD')
    benchmarks.append('XOR')
    benchmarks.append('SELOSC')
    benchmarks.append('SRFlipFlop')

    for name in benchmarks:
        filepath = os.path.join(dir, name+'.qca')
        mm_dir = os.path.join(dir, name)
        R_MAX = 2.5
        QCAnx = QCANetworkX(filepath, ancilla=True, r_max=R_MAX)
        QCAnx.draw_qca()

        comments = "Source: %s\nPolarizations: %s" % (name, QCAnx.QCA.pols)
        write_networkx(QCAnx, pos=QCAnx.pos, mtx_name=name, mm_dir=mm_dir, comment=comments)
