import networkx as nx
import matplotlib.pyplot as plt
from embedding_methods.utilities.graph_mmio import read_networkx

from dimod.reference.samplers.exact_solver import ExactSolver
from dimod.reference.samplers.simulated_annealing import SimulatedAnnealingSampler

%matplotlib tk

bench_dir = '../../benchmarks/'
mm_name = 'AND4'

Sg = read_networkx(mm_name,mm_dir=bench_dir)
pos = Sg.graph['pos']
plt.figure(1)
nx.draw(Sg, pos=pos, with_labels=True)
plt.gca().invert_yaxis()
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

#response = sampler.sample_ising(h,J,num_reads=500)
response = sampler.sample_ising(h,J)

import pickle
with open('problem.pkl','wb') as fp:
    pickle.dump(Sg, fp)

with open('response.pkl','wb') as fp:
    pickle.dump(response, fp)

energies = [datum.energy for datum in response.data()]
plt.figure(2)
_ = plt.hist(energies, bins=100)
plt.savefig('response.png')

with open('response.pkl','rb') as fp:
    response = pickle.load(fp)

data10 = []
for i in range(10):
    data10.append( next(response.data()) )

sample = data10[0].sample

plt.figure(3)
nx.draw(Sg, pos=pos, labels=sample, with_labels=True)
_ = nx.draw_networkx_edge_labels(Sg,pos=pos,edge_labels=J)
plt.show()


from math import exp

beta = 1.0
Z=0.0
for datum in response.data():
    Ei = datum.energy
    Z += exp(-beta*Ei)

print(Z)

pdf = []
for datum in response.data():
    Ei = datum.energy
    pdf.append( (Ei**2)*exp(-beta*Ei) )

plt.plot(pdf)
plt.show
