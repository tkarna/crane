#!/usr/bin/env python

from distutils.core import setup

setup(name='crane',
      version='0.0.1',
      description='Pre- and post processing tools for ocean modeling products',
      author='Tuomas Karna',
      author_email='tuomas.karna@gmail.com',
      url='https://bitbucket.org/tkarna/crane',
      packages=['crane', 'crane.data', 'crane.plotting',
		'crane.files'],
     )

