from operator import indexOf
from unicodedata import name
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
parser.add_option('-o', '--out-datacard', dest = "outdatacard", default = "converted_datacard.txt")
parser.add_option('-s', 'shape-file', dest = 'shapefile', default = 'shapes.root')

with open(args[0]) as serialized:
    spec = json.load(serialized)
workspace = pyhf.Workspace(spec)
model = workspace.model()

channels = model.config.channels
channel_bins = model.config.channel_nbins
samples = model.config.samples
modifiers = model.config.modifiers
lumi = model.config.parameters

DC = Datacard()
DC.processes = samples

uproot.recreate(options.shapefile)

def addChannels():
    DC.bins = channels
    for idxc, channel in enumerate(channel_bins.keys()):
        if channel_bins[channel] == 1: ##single bin
            DC.obs.update({channel: workspace.data(model, include_auxdata=False)[idxc]})
        else:
            DC.obs.update({channel: -1})
            with uproot.update(options.shapefile) as file:
                            h_data = hist.Hist.new.Regular(channel_bins[channel], 0, channel_bins[channel]).Weight() 
                            h_data[...] = np.stack([workspace.data(model, include_auxdata=False)[idxc], [0 for i in channel_bins[channel]]], axis=-1)
                            file[channel/'obs_data'] = h_data
                            


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
                DC.exp[channel].update(
                    {
                        spec["channels"][idxc]["samples"][idxs]["name"]: 
                        -1
                    }
                )
                mods = [
                    spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["type"]
                    for i in range(spec["channels"][idxc]["samples"][idxs]["modifiers"])
                ]
                if "histosys" in mods:
                    DC.shapeMap.update({channel: {spec["channels"][idxc]["samples"][idxs]["name"]: [options.shapeFile, channel/spec["channels"][idxc]["samples"][idxs]["name"], channel/spec["channels"][idxc]["samples"][idxs]["name"]+"_$SYSTEMATIC" ]}})
                else:
                    DC.shapeMap.update({channel: {spec["channels"][idxc]["samples"][idxs]["name"]: [options.shapeFile, channel/spec["channels"][idxc]["samples"][idxs]["name"]]}})
                if "staterror" in mods:
                    h_data = hist.Hist.new.Regular(channel_bins[channel], 0, channel_bins[channel]).Weight()  
                    h_data[...] = np.stack([spec["channels"][idxc]["samples"][idxs]["data"], spec["channels"][idxc]["samples"][idxs]["modifiers"][mods.index("staterror")]["data"]], axis=-1)
                elif "shapesys" in mods:
                    h_data = hist.Hist.new.Regular(channel_bins[channel], 0, channel_bins[channel]).Weight() 
                    h_data[...] = np.stack([spec["channels"][idxc]["samples"][idxs]["data"], spec["channels"][idxc]["samples"][idxs]["modifiers"][mods.index("shapesys")]["data"]], axis=-1)
                else:
                    h_data = hist.Hist.new.Regular(channel_bins[channel], 0, channel_bins[channel]).Weight()  
                    h_data[...] = np.stack([spec["channels"][idxc]["samples"][idxs]["data"], [0 for i in channel_bins[channel]]], axis=-1)
            with uproot.update(options.shapefile) as file:
                file[channel/spec["channels"][idxc]["samples"][idxs]["name"]] = h_data
            
def addMods():
    for modifier in modifiers:
        if "normsys" in modifier[1]:
            for idxc, channel in enumerate(channels):
                DC.systs.append((modifier[0], False, 'lnN', [], {channel : {}}))
                for sample in samples:
                    if sample in [spec["channels"][idxc]["samples"][i]["name"] for i in spec["channels"][idxc]["samples"]]: 
                        for i in DC.systs:
                            i[4][channel].update({sample : 0.0})           
        if "histosys" in modifier[1]:
            for channel in channels:
                DC.systs.append((modifier[0], False, 'shape', [], {channel: {}}))
                for sample in samples:
                    if sample in [spec["channels"][idxc]["samples"][i]["name"] for i in spec["channels"][idxc]["samples"]]: 
                        for i in DC.systs:
                            i[4][channel].update({sample : 0.0})  
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate([spec["channels"][idxc]["samples"][i]["name"] for i in spec["channels"][idxc]["samples"]]):
            mods = []
            for i in range(spec["channels"][idxc]["samples"][idxs]["modifiers"]):
                mods.append([spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["name"], spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["type"], spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["data"]])
            for syst in DC.systs:
                name = syst[0]
                type = syst[2]
                if name in mods: ##if systematic is a modifier for this sample
                    if "lnN" in type:
                        syst[4][channel].update({sample: mods[mods.index(name)][2]["hi"]})
                    if "shape" in type:
                        syst[4][channel].update({sample: 1.0})
                        

                    
                        
            
                

