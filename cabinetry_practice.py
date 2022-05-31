import cabinetry
cabinetry.set_logging()
import copy
import json
import pathlib

import boost_histogram as bh
import numpy as np
import pyhf
from pyhf.contrib.utils import download

config_ex = cabinetry.configuration.load("config_example.yml")
print(config_ex['General'])





