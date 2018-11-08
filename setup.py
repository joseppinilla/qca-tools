#!/usr/bin/env python

from setuptools import setup

packages = ['qca_tools',
			'qca_tools.composite',
			'qca_tools.drawing',
		]

install_requires = ['networkx>=2.0,<3.0',
                    'decorator>=4.1.0,<5.0.0',
                    'dimod>=0.6.8,<0.8.0']


setup(name='qca_tools',
      version='0.0.1',
      description='QCA Tools',
      long_description="Collection of tools for parsing and simulation of QCA circuits.",
      author='Jose Pinilla',
      author_email='jpinilla@ece.ubc.ca',
      url='https://github.com/joseppinilla/qca_tools',
      packages=packages,
      platforms='any',
      install_requires=install_requires,
      license='MIT'
     )
