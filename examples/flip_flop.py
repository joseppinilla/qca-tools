import networkx as nx
import matplotlib.pyplot as plt
from embedding_methods.utilities.graph_mmio import read_networkx

bench_dir = '../benchmarks/'
mm_name = 'SRFlipFlop'

Gnx = read_networkx(mm_name,mm_dir=bench_dir)
plt.figure(1)
nx.draw(Gnx, pos=Gnx.graph['pos'])
plt.show()
