import dimod
from os.path import isfile
from qca_tools.qca_network import QCANetwork

class QCAComposite(dimod.ComposedSampler):

    def __init__(self, child_sampler, qca_filename, pols={}, use_ancilla=True):
        if not isfile(qca_filename):
            raise ValueError("QCA input file not found")
        self._qca_filename = qca_filename
        self._use_ancilla = use_ancilla
        self._children = [child_sampler]
        self._QCA = QCANetwork(qca_filename)

    @property
    def children(self):
        """list: Children property inherited from :class:`dimod.Composite` class.

        For an instantiated composed sampler, contains the single wrapped structured sampler.

        .. _configuration: http://dwave-cloud-client.readthedocs.io/en/latest/#module-dwave.cloud.config

        """
        return self._children

    @property
    def parameters(self):
        """dict[str, list]: Parameters in the form of a dict.

        For an instantiated composed sampler, keys are the keyword parameters accepted by the child sampler.

        .. _configuration: http://dwave-cloud-client.readthedocs.io/en/latest/#module-dwave.cloud.config

        """
        # does not add or remove any parameters
        param = self.child.parameters.copy()
        return param

    @property
    def properties(self):
        """dict: Properties in the form of a dict.

        For an instantiated composed sampler, contains one key :code:`'child_properties'` that
        has a copy of the child sampler's properties.

        .. _configuration: http://dwave-cloud-client.readthedocs.io/en/latest/#module-dwave.cloud.config

        """
        properties = {'child_properties': self.child.properties.copy()}
        properties['qca_filename'] = self._qca_filename
        properties['use_ancilla'] = self._use_ancilla
        return properties

    def get_outputs(self):
        outputs = self.outputs

    def get_qca_network(self):
        return self._QCA

    def sample(self, **parameters):
        child = self.child
        QCA = self._QCA
        response = child.sample(QCA)
        return response
