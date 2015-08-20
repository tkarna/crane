#!/usr/bin/env python
from distutils.core import setup

# install all excecutable files under ./scripts to PREFIX/bin/
import glob
import os
scripts = []
for f in glob.glob('./scripts/*'):
    if os.path.isfile(f) and os.access(f, os.X_OK):
        scripts.append(f)

setup(name='crane',
      version='0.0.1',
      description='Pre- and post processing tools for ocean modeling products',
      author='Tuomas Karna',
      author_email='tuomas.karna@gmail.com',
      url='https://bitbucket.org/tkarna/crane',
      packages=['crane', 'crane.data', 'crane.plotting',
                'crane.files'],
      scripts=scripts,
      )

