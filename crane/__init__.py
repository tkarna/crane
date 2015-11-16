import sys
import os
import glob
import datetime
import numpy as np  # TODO import as numpy
import scipy
import time as timeMod

import matplotlib
# change backend if no display present (e.g. on supercomputers)
if 'TACC_SYSTEM' in os.environ and os.environ['TACC_SYSTEM'] == 'stampede':
    matplotlib.use('Agg', warn=False)
elif 'stccmop' in os.environ.get('HOSTNAME', []):
    matplotlib.use('Agg', warn=False)
elif 'edison' in os.environ.get('HOSTNAME', []):
    matplotlib.use('Agg', warn=False)
elif 'hopper' in os.environ.get('HOSTNAME', []):
    matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt

# bring useful functions/variables to highest level in namespace
from utility import *
from physicalVariableDefs import *
