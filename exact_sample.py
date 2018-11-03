import os
import pickle
import networkx as nx
import matplotlib.pyplot as plt

from embedding_methods.utilities.graph_mmio import read_networkx
from embedding_methods.utilities.graph_mmio import write_networkx

from qca_to_mm import QCANetwork

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
    nx.draw(Sg, pos=pos, with_labels=True)
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
        response = sampler.sample_ising(h,J,num_reads=100)


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
    nx.draw(Sg, pos=pos, labels=sample, with_labels=True)
    _ = nx.draw_networkx_edge_labels(Sg,pos=pos,edge_labels=J)
    plt.gca().invert_yaxis()
    plt.savefig(dir + 'ground_state.png')

if __name__ == "__main__":

    #benchmarks = ['NOT']
    #benchmarks = ['NOT', 'mux2to1']
    #benchmarks = ['NOT', 'mux2to1', 'SRFlipFlop']
    benchmarks = ['NOT', 'mux2to1', 'SRFlipFlop', 'AND4', 'half_adder']

    pols = {    'NOT':[{'A':-1}, {'A':1}],
                'mux2to1':[{'I0':-1,'I1':1,'S':-1}, {'I0':-1,'I1':1,'S':1}, {'I0':1,'I1':-1,'S':-1}],
                'SRFlipFlop':[{'S':1,'R':-1}, {'S':1,'R':1}],
                'AND':[{'A':-1, 'B':-1, 'C':1, 'D':-1}, {'A':-1, 'B':1, 'C':1, 'D':-1}, {'A':-1, 'B':-1, 'C':-1, 'D':-1}],
                'half_adder':[{'A':1, 'B':1}, {'A':1, 'B':-1}]
            }

    for name in benchmarks:
        for i, pol in enumerate(pols[name]):
            dir = bench_dir + name + solver + str(i) + '/'
            if not os.path.exists(dir):
                os.makedirs(dir)

            G = QCANetwork(qca_file=bench_dir+name+'.qca', pols=pol)
            pos = G.qca_layout()
            comments = "Source: %s\nPolarizations: %s" % (name, pol)

            write_networkx(G, pos=pos, mtx_name=name, mm_dir=dir, comment=comments)

            sample(name, dir)
