""" Example usage of the D-Wave composites to obtain the exact
    solution for a QCA input network file.
"""

from qca_tools import QCAComposite
from dimod.reference.samplers import ExactSolver

bench_dir = '../examples/benchmarks/'
qca_name = 'NOT_FT.qca'
qca_filename = bench_dir + qca_name

# Create sampler
child_sampler = ExactSolver()
qca_sampler = QCAComposite(child_sampler, qca_filename)

# Collect samples
response = qca_sampler.sample()

# First Energy level is the ground state, only because it's the exact solution
ground_state, energy, _ = response.first

# Plot Energies
energies = [datum.energy for datum in response.data()]
_ = plt.hist(energies, bins=100)
plt.show()
