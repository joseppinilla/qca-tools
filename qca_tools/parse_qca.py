#!/usr/bin/env python

#---------------------------------------------------------
# Name: parse_qca.py
# Purpose: Parsing functions for QCADesigner files
# Author: Jacob Retallick
# Created: 2015.10.22
# Last Modified: 2015.10.22
#---------------------------------------------------------

# NOTE
# the original parse script no longer seems to work (change in networkx?)
# for the purposes of the embedder, we don't need to consider clock zones so
# I have simplified the parseing script to remove that functionality.

import re

import networkx as nx
import numpy as np
#from auxil import getEk, CELL_FUNCTIONS, CELL_MODES
from itertools import combinations

from functools import reduce

## general global parameters

CELL_FUNCTIONS = {'QCAD_CELL_NORMAL': 0,
                  'QCAD_CELL_INPUT': 1,
                  'QCAD_CELL_OUTPUT': 2,
                  'QCAD_CELL_FIXED': 3}

CELL_MODES = {'QCAD_CELL_MODE_NORMAL': 0,
              'QCAD_CELL_MODE_CROSSOVER': 1,
              'QCAD_CELL_MODE_VERTICAL': 2,
              'QCAD_CELL_MODE_CLUSTER': 3}

R_MAX = 2.1         # max cell-cell interaction range (rel to grid spacing)
EK_THRESH = 1e-3    # threshold for included Ek, relative to max(abs(Ek))
X_ROUND = 4         # places to round to when deciding if cell is rotated

# physical parameters
eps0 = 8.85412e-12  # permittivity of free space
epsr = 12.          # relative permittivity
q0 = 1.602e-19      # elementary charge

def getEk(c1, c2, DR=2):
    '''Compute the kink energy using the qdot positions (in nm) of two cells.
    If the cell displacement is greater than DR return False'''

    # check cell-cell range
    dx = c1['x']-c2['x']
    dy = c1['y']-c2['y']
    if dx*dx+dy*dy > DR*DR:
        return False

    qdots_1 = c1['qdots']
    qdots_2 = c2['qdots']

    # compute displacements

    x1 = [qd['x'] for qd in qdots_1]
    y1 = [qd['y'] for qd in qdots_1]

    x2 = [qd['x'] for qd in qdots_2]
    y2 = [qd['y'] for qd in qdots_2]

    X1 = np.array([x1, y1]).T.reshape([4, 1, 2])
    X2 = np.array([x2, y2]).T.reshape([1, 4, 2])

    R = np.sqrt(np.sum(pow(X1-X2, 2), axis=2))

    if np.min(R) == 0:
        print ('qdot overlap detected')
        return 0.

    # QCADesigner orders qdots either CW or CCW so same and diff configurations
    # are always alternating indices.

    Q = np.array([1, -1, 1, -1])    # template for charge arrangement

    Q = np.outer(Q, Q)

    Ek = -1e9*q0*np.sum((Q/R))/(8*np.pi*eps0*epsr)

    return Ek

### FILE PROCESSING

def build_hierarchy(fn):
    '''Build a dict hierarchy containing all objects, their parameters, and
    childen.'''

    fp = open(fn, 'r')

    linemap = lambda s: s.replace(',', '.')

    # general re expression. may need to change if future format changes
    re_start = re.compile('^\[.+\]$')
    re_term = re.compile('^\[#.+\]$')

    hier = {'label': 'Hierarchy', 'children': [], 'vars': {}}

    key_stack = ['Hierarchy']  # stack of active keys, pop of top of stack
    dict_stack = [hier]   # stack of corresponding dict objects.

    line_cnt = 0
    for line in fp:
        line = linemap(line)
        line_cnt += 1
        line = line.strip()  # remove endline and possible whitespace

        # must check object termination first
        if re_term.match(line):
            key = line[2:-1]
            if key_stack[-1] == key:
                d = dict_stack.pop()
                key_stack.pop()
                try:
                    dict_stack[-1]['children'].append(d)
                except:
                    print('Somehow over-popped dict_stack...')
                    return None
            else:
                print('Start-end mismatch in line {0}'.format(line_cnt))
                return None

        # for a new object, create a new dict template
        elif re_start.match(line):
            key = line[1:-1]
            key_stack.append(key)
            d = {'label': key, 'children': [], 'vars': {}}
            dict_stack.append(d)

        # otherwise check for new variable to add to most recent dict
        else:
            if '=' in line:
                var, val = line.split('=')
                dict_stack[-1]['vars'][var] = val
    fp.close()

    return hier


def proc_hierarchy(hier):
    '''Process the extracted data hierarchy to extract useful information. In
    the current information, we are interested in the overall cell grid spacing
    (for deciding on the range of included cell) and the properties of each
    cell in the circuit'''

    cells = []
    spacing = None

    # hierarchy should only have two children: VERSION and TYPE:DESIGN. The
    # former might be useful in later implentations for selecting formatting
    # options but for now all we care about are the DESIGN objects

    hier = [child for child in hier['children']
            if child['label'] == 'TYPE:DESIGN'][0]

    # for now assert that there can be only one cell layer, no vertical x-over
    layers = [child for child in hier['children']
              if child['label'] == 'TYPE:QCADLayer']

    # isolate cell layers
    cell_layers = [layer for layer in layers if layer['vars']['type'] == '1']

    # merge cell layers, will lead to qdot conflict if vertical x-over
    cell_dicts = [layer['children'] for layer in cell_layers]
    cell_dicts = reduce(lambda x, y: x+y, cell_dicts)

    # get grid-spacing (average cell bounding box)
    cx = float(cell_dicts[0]['vars']['cell_options.cxCell'])
    cy = float(cell_dicts[0]['vars']['cell_options.cyCell'])

    spacing = np.sqrt(cx*cy)
    # create cell objects
    cells = []

    for cd in cell_dicts:
        cell = {}

        # cell type
        cell['cf'] = cd['vars']['cell_function']
        cell['cm'] = cd['vars']['cell_options.mode']
        cell['clk'] = int(cd['vars']['cell_options.clock'])

        # just for show sol
        cell['cx'] = float(cd['vars']['cell_options.cxCell'])
        cell['cy'] = float(cd['vars']['cell_options.cyCell'])

        # position, first child will be the QCADesignObject
        design_object = cd['children'][0]
        cell['x'] = float(design_object['vars']['x'])
        cell['y'] = float(design_object['vars']['y'])

        # quantum dots
        qdot_dicts = [child for child in cd['children']
                      if child['label'] == 'TYPE:CELL_DOT']

        qdots = []
        for d in qdot_dicts:
            dot = {}
            dot['x'] = float(d['vars']['x'])
            dot['y'] = float(d['vars']['y'])
            dot['c'] = float(d['vars']['charge'])
            qdots.append(dot)
        cell['qdots'] = qdots

        # cell name label
        for child in cd['children']:
            if child['label'] == 'TYPE:QCADLabel':
                cell['name'] = child['vars']['psz']

        # determine if cell is rotated, will have three x values
        x = set([round(dt['x'], X_ROUND) for dt in qdots])
        if len(x) == 3:
            cell['rot'] = True
        elif len(x) == 2:
            cell['rot'] = False
        else:
            print('Could not decide cell rotation')
            cell['rot'] = False

        # keep track of polarization if cell is fixed: don't rely on labels
        if cell['cf'] == 'QCAD_CELL_FIXED':
            pol = qdots[0]['c']+qdots[2]['c']-qdots[1]['c']-qdots[3]['c']
            pol /= qdots[0]['c']+qdots[2]['c']+qdots[1]['c']+qdots[3]['c']
            cell['pol'] = pol

        cells.append(cell)

    return cells, spacing


## CIRCUIT PROCESSING

def build_J(cells, spacing, r_max=R_MAX):
    '''Build the J matrix for the given circuit. Restricts the interaction
    distance to r_max but does not apply any adjacency contraints'''

    N = len(cells)

    # contruct connectivvity matrix
    J = np.zeros([N, N], dtype=float)
    DR = r_max*spacing
    for i,j in combinations(range(N), 2):
        Ek = getEk(cells[i], cells[j], DR=DR)
        if Ek:
            J[i,j] = J[j,i] = Ek

    # remove very weak interactions
    J = -J*(np.abs(J) >= np.max(np.abs(J)*EK_THRESH))

    return J

def reorder_cells(cells, J):
    '''Renumber cells by position rather than the default QCADesigner placement
    order. Cells ordered by the tuple (zone, y, x)'''

    keys = {}

    # assign sortable tuples for each cell
    for ind, cell in enumerate(cells):
        keys[ind] = (cell['y'], cell['x'])

    order = list(zip(*sorted([(keys[i], i) for i in keys])))[1]

    # relabel cells and reorder the J matrix
    cells = [cells[i] for i in order]
    J = J[order, :][:, order]

    for i in range(len(cells)):
        cells[i]['num'] = i
        cells[i]['number'] = i

    return cells, J


## MAIN FUNCTION

def parse_qca_file(fn, r_max=R_MAX, verbose=False):
    '''Parse a QCADesigner file to extract cell properties. Returns an ordered
    list of cells, the QCADesigner grid spacing in nm, a list structure of the
    indices of each clock zone (propogating from inputs), and a coupling matrix
    J which contains the Ek values for cells within a radius of R_MAX times the
    grid spacing'''

    # build data hierarchy
    hier = build_hierarchy(fn)

    # extract useful information from data hierarchy
    cells, spacing = proc_hierarchy(hier)

    if verbose:
        print('Parsed QCA file...')

    # construct J matrix
    J = build_J(cells, spacing, r_max=r_max)

    # reorder cells by zone and position
    cells, J = reorder_cells(cells, J)

    return cells, spacing, J

if __name__ == "__main__":
    bench_dir = './benchmarks/'
    fn = 'S_R_Flip_Flop.qca'

    parse_qca_file(bench_dir+fn)
