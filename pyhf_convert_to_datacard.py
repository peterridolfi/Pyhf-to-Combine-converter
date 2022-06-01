import pyhf
import uproot
import numpy as np
import hist
from hist import Hist
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import json
import string

with open("C:/Users/peter/iris_project/shape_example_ws.json") as serialized:
    spec = json.load(serialized)
workspace = pyhf.Workspace(spec)
sig = spec['channels'][0]['samples'][0]['data']
bkg = spec['channels'][0]['samples'][1]['data']
bkg_err = spec['channels'][0]['samples'][1]['modifiers'][1]['data']
bkg_histosys_up = spec['channels'][0]['samples'][1]['modifiers'][0]['data']['hi_data']
bkg_histosys_down = spec['channels'][0]['samples'][1]['modifiers'][0]['data']['lo_data']
observation = spec['observations'][0]['data']
h_background = hist.Hist.new.Regular(3, 0, 3).Weight()  # 3 bins
h_background[...] = np.stack([bkg, bkg_err], axis=-1)  # set yields + uncertainties
h_signal = hist.Hist.new.Regular(3, 0, 3).Weight()  # 3 bins
h_signal[...] = np.stack([sig, [0, 0, 0]], axis=-1)
h_background_up = hist.Hist.new.Regular(3, 0, 3).Weight()  # 3 bins
h_background_up[...] = np.stack([bkg_histosys_up, bkg_err], axis=-1)
h_background_down = hist.Hist.new.Regular(3, 0, 3).Weight()  # 3 bins
h_background_down[...] = np.stack([bkg_histosys_down, bkg_err], axis=-1)
h_obs = hist.Hist.new.Regular(3, 0, 3).Weight()  # 3 bins
h_obs[...] = np.stack([observation, [0, 0, 0]], axis=-1)
with uproot.recreate("C:/Users/peter/iris_project/example_shapes.root") as file:
    file["bin1/background"] = h_background
    file["bin1/signal"] = h_signal
    file["bin1/background_alphaUp"] = h_background_up
    file["bin1/background_alphaDown"] = h_background_down
    file["bin1/data_obs"] = h_obs
    
    




