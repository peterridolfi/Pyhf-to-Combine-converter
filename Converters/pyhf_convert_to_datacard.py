from CMSSW_11_2_0.pyhf import parameters
import pyhf
import uproot
import numpy as np
import hist
from hist import Hist
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import json
import string
import HiggsAnalysis.CombinedLimit.DatacardParser as DP
from modulefinder import Module
from optparse import OptionParser
from HiggsAnalysis.CombinedLimit.Datacard import Datacard

parser = OptionParser()
parser.add_option("-o", "--out-file", dest="outfile", default="converted_datacard.txt")
options, args = parser.parse_args()


with open("C:/Users/peter/iris_project/shape_example_ws.json") as serialized:
    spec = json.load(serialized)
workspace = pyhf.Workspace(spec)
model = workspace.model()

channels = model.config.channels
channel_bins = model.config.channel_nbins
samples = model.config.samples
modifiers = model.config.modifiers
lumi = model.config.parameters

DC = Datacard()
DC.bins = channels
DC.processes = samples


def addSamples():
    for idxc, channel in enumerate(channel_bins.keys()):
        if channel_bins[channel] == 1:  ##counting
            for idxs in range(spec["channels"][idxc]["samples"]):
                DC.exp[channel].update(
                    {
                        spec["channels"][idxc]["samples"][idxs]["name"]: 
                        spec["channels"][idxc]["samples"][idxs]["data"]
                    }
                )
        else: ##shapes
            DC.hasShapes = True
            for idxs in range(spec["channels"][idxc]["samples"]):
                mods = [
                    spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["type"]
                    for i in range(spec["channels"][idxc]["samples"][idxs]["modifiers"])
                ]
                if "staterror" in mods:
                    h_data = hist.Hist.new.Regular(3, 0, 3).Weight()  # 3 bins
                    h_data[...] = np.stack([spec["channels"][idxc]["samples"][idxs]["data"], spec["channels"][idxc]["samples"][idxs]["modifiers"][mods.index("staterror")]["data"]], axis=-1)
                elif "shapesys" in mods:
                    h_data = hist.Hist.new.Regular(3, 0, 3).Weight()  # 3 bins
                    h_data[...] = np.stack([spec["channels"][idxc]["samples"][idxs]["data"], spec["channels"][idxc]["samples"][idxs]["modifiers"][mods.index("shapesys")]["data"]], axis=-1)
                else:
                    h_data = hist.Hist.new.Regular(3, 0, 3).Weight()  # 3 bins
                    h_data[...] = np.stack([spec["channels"][idxc]["samples"][idxs]["data"], [0 for i in channel_bins[channel]]], axis=-1)
            with uproot.recreate("C:/Users/peter/iris_project/example_shapes.root") as file:
                file[channel/spec["channels"][idxc]["samples"][idxs]["name"]] = h_data

    
