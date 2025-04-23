# This script takes command line arguments for the following parameters:
# Filename (-f), random seed (-s), number of iterations (-i).

import argparse
import sys
sys.path.insert(1, "../code/experiment_code")
sys.path.insert(1, "../code/utils")

from experiment_utils import *
from celldata import *
from effects import *
from n_effective import *
from experiment_dumps import *

# import experiment code (needs the above imports)
from FDR_experimentV1 import *

# parse command line arguments
parser = argparse.ArgumentParser(description='Run experiment')
parser.add_argument('-f', '--filename', type=str, help='Input filename with entire path')
# random seed default is 0
parser.add_argument('-s', '--seed', type=int, default=0, help='Random seed')
# number of iterations default is 5000
parser.add_argument('-i', '--iterations', type=int, default=5000, help='Number of iterations')


args = parser.parse_args()

# set random seed
np.random.seed(args.seed)

