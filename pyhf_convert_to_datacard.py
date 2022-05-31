import pyhf
import uproot
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import json
import string

with open("C:/Users/peter/iris_project/shape_example_ws.json") as serialized:
    spec = json.load(serialized)
workspace = pyhf.Workspace(spec)
sig = spec['channels'][0]['samples'][0]['data']
bkg = spec['channels'][0]['samples'][1]['data']
bkg_shapesys = spec['channels'][0]['samples'][1]['modifiers'][0]['data']
bkg_histosys_up = spec['channels'][0]['samples'][1]['modifiers'][1]['data']['hi_data']
bkg_histosys_down = spec['channels'][0]['samples'][1]['modifiers'][1]['data']['lo_data']
observation = spec['observations'][0]['data']
with uproot.recreate("C:/Users/peter/iris_project/example_shapes.root") as file:
    file["bin1/background"] = np.histogram(range(3), weights=bkg, bins = range(4))
    file["bin1/signal"] = np.histogram(range(3), weights=sig, bins = range(4))
    file["bin1/background_uncorr"] = np.histogram(range(3), weights=bkg_shapesys, bins = range(4))
    file["bin1/background_alphaUp"] = np.histogram(range(3), weights=bkg_histosys_up, bins = range(4))
    file["bin1/background_alphaDown"] = np.histogram(range(3), weights=bkg_histosys_down, bins = range(4))
    file["bin1/data_obs"] = np.histogram(range(3), weights=observation, bins = range(4))
    
    




