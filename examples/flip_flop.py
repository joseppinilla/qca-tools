import networkx as nx
import matplotlib.pyplot as plt
from embedding_methods.utilities.graph_mmio import read_networkx

from dimod.reference.samplers.exact_solver import ExactSolver
from dimod.reference.samplers.simulated_annealing import SimulatedAnnealingSampler

bench_dir = '../benchmarks/'
mm_name = 'SRFlipFlop'

Sg = read_networkx(mm_name,mm_dir=bench_dir)
plt.figure(1)
nx.draw(Sg, pos=Sg.graph['pos'], with_labels=True)
plt.show()

h = {}
J = {}
for u,v,data in Sg.edges(data=True):
    if u==v:
        h[u] = data['weight']
    else:
        J[(u,v)] = data['weight']

#sampler = SimulatedAnnealingSampler()
sampler = ExactSolver()

#response = sampler.sample_ising(h,J,num_reads=200)
response = sampler.sample_ising(h,J)

import pickle
with open('problem.pkl','wb') as fp:
    pickle.dump(Sg, fp)

with open('response.pkl','wb') as fp:
    pickle.dump(response, fp)

energies = [datum.energy for datum in response.data()]
plt.figure(2)
_ = plt.hist(energies, bins=100)
plt.show()
