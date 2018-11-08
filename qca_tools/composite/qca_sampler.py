import dimod
from os.path import isfile
from qca_tools import QCANetwork

class QCAComposite(QCANetwork, dimod.ComposedSampler):

    def __init__(self, child_sampler, qca_filename, pols={}, use_ancilla=True):
        if not isfile(qca_filename):
            raise ValueError("QCA input file not found")
        self._qca_filename = qca_filename
        self._use_ancilla = ancilla
        self._children = [child_sampler]

        QCANetwork.__init__(self)


    def _get_outputs(self):
        outputs = self.outputs


    def sample(self, **parameters):
        child = self.child
        bqm = self._bqm

        response = child.sample(bqm)
