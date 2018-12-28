""" Parse a *.qca file from the provided benchmarks
and draw using the provided methods.
"""
import os
import matplotlib.pyplot as plt
from qca_tools.qca_network import QCANetworkX

# Set-up file path
dir= './benchmarks/'
qca_name = 'XOR'
qca_filename = qca_name + '.qca'
filepath = os.path.join(dir, qca_filename)

# Parse
R_MAX = 2.5
QCAnx = QCANetworkX.from_qca_file(filepath, ancilla=True, r_max=R_MAX)

# Draw
QCAnx.draw_qca(with_labels=True)
plt.title(qca_name)
plt.savefig(os.path.join(dir, qca_name, qca_name + '.svg'))
