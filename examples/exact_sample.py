import os
import pickle
import itertools
import networkx as nx
import matplotlib.pyplot as plt

from embedding_methods.utilities.graph_mmio import read_networkx
from embedding_methods.utilities.graph_mmio import write_networkx

from qca_tools.qca_network import QCANetworkX

from dimod.reference.samplers.exact_solver import ExactSolver
from dimod.reference.samplers.simulated_annealing import SimulatedAnnealingSampler

bench_dir = './benchmarks/'
#solver = 'dwave'
#solver = 'exact'
solver = 'sa'


def sample(mm_name, dir):

    Sg = read_networkx(mm_name,mm_dir=dir)
    pos = Sg.graph['pos']

    plt.figure(1)
    plt.clf()
    nx.draw(Sg, pos=pos, with_labels=True, node_shape='s')
    plt.gca().invert_yaxis()
    plt.savefig(dir + 'problem.png')

    h = {}
    J = {}
    for u,v,data in Sg.edges(data=True):
        if u==v:
            h[u] = data['weight']
        else:
            J[(u,v)] = data['weight']

    if solver=='sa':
        sampler = SimulatedAnnealingSampler()
    elif solver=='exact':
        sampler = ExactSolver()
    elif solver=='dwave':
        #TODO: EmbeddingComposite
        sampler = DWaveSampler()

    if solver=='exact':
        response = sampler.sample_ising(h,J)
    else:
        response = sampler.sample_ising(h,J,num_reads=1)


    with open(dir + 'problem.pkl','wb') as fp:
        pickle.dump(Sg, fp)

    with open(dir + 'response.pkl','wb') as fp:
        pickle.dump(response, fp)

    # energies = [datum.energy for datum in response.data()]
    # plt.figure(2)
    # _ = plt.hist(energies, bins=100)
    # plt.savefig('response.png')

    datum = next(response.data())
    sample = datum.sample

    plt.figure(3)
    plt.clf()
    nx.draw(Sg, pos=pos, labels=sample, with_labels=True, node_shape='s')
    _ = nx.draw_networkx_edge_labels(Sg,pos=pos,edge_labels=J)
    plt.gca().invert_yaxis()
    plt.savefig(dir + 'ground_state.png')

if __name__ == "__main__":

    EXHAUSTIVE = False

    benchmarks = []
    benchmarks.append('NOT_FT')
    benchmarks.append('MAJ5_A')
    benchmarks.append('MAJ5_B')
    benchmarks.append('MUX')
    benchmarks.append('COPLANARX')

    # Vector inputs
    vectors = { 'NOT_FT':       [{'A':-1}],
                'MAJ5_A':       [{'A':-1, 'B':1, 'C':1, 'D':-1, 'E':1}],
                'MAJ5_B':       [{'A':-1, 'B':1, 'C':1, 'D':-1, 'E':1}],
                'MUX':          [{'S':-1,'A':1,'B':-1}],
                'COPLANARX':    [{'A':-1,'X':1}]
            }

    for name in benchmarks:
        if EXHAUSTIVE:
            vector0 = vectors[name][0]
            vector_size = len(vector0)
            inputs = list(vector0.keys())

            combs = itertools.product([-1,1], repeat=vector_size)
            pols = [  dict(zip(inputs,comb)) for comb in combs ]
        else:
            pols = vectors[name]

        # run
        for i, pol in enumerate(pols):
            dir = bench_dir + name + solver + str(i) + '/'
            if not os.path.exists(dir):
                os.makedirs(dir)

            G = QCANetworkX(bench_dir+name+'.qca', pols=pol)
            G.draw_qca()
            comments = "Source: %s\nPolarizations: %s" % (name, pol)

            write_networkx(G, pos=G.pos, mtx_name=name, mm_dir=dir, comment=comments)

            sample(name, dir)
