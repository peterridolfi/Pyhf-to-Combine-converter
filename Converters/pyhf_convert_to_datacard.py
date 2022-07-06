from ctypes import sizeof
from operator import indexOf
from unicodedata import name
from numpy.core.fromnumeric import size
from pyhf import parameters
import pyhf
import uproot
import numpy as np
import hist
from hist import Hist
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import json
import string
from HiggsAnalysis.CombinedLimit.DatacardParser import *
from modulefinder import Module
from optparse import OptionParser
from HiggsAnalysis.CombinedLimit.Datacard import Datacard

parser = OptionParser()
parser.add_option("-o", "--out-file", dest="outfile", default="converted_datacard.txt")
parser.add_option(
    "-O", "--out-datacard", dest="outdatacard", default="converted_datacard.txt"
)
parser.add_option("-s", "--shape-file", dest="shapefile", default="shapes.root")
options, args = parser.parse_args()

with open(args[0]) as serialized:
    spec = json.load(serialized)
workspace = pyhf.Workspace(spec)
model = workspace.model()

channels = model.config.channels
channel_bins = model.config.channel_nbins
samples = model.config.samples
modifiers = model.config.modifiers
lumi = model.config.parameters

##generate list of mods without normfactors, shapesys, staterror
systs = []
for mod in modifiers:
    if 'normfactor' not in mod:
        if 'shapesys' not in mod:
            if 'staterror' not in mod:
                systs.append(mod)


DC = Datacard()


file = uproot.recreate(options.shapefile)


def addChannels():
    DC.bins = channels
    count = 0
    for idxc, channel in enumerate(channel_bins.keys()):
        DC.exp.update({channel : {}})
        if channel_bins[channel] == 1:  ##single bin
            DC.obs.update({channel: workspace.data(model, include_auxdata=False)[idxc]})
        else:
            DC.obs.update({channel: -1})
            h_data = hist.Hist.new.Regular(
                channel_bins[channel], 0, channel_bins[channel]
            ).Weight()
            h_data[...] = np.stack(
                [
                    workspace.data(model, include_auxdata=False)[count:count+channel_bins[channel]],
                    [0 for i in range(channel_bins[channel])],
                ],
                axis=-1,
            )
            count = count + 1
            file[channel+ "/data_obs"] = h_data


def addSamples():
    DC.processes = samples
    for idxc, channel in enumerate(channel_bins.keys()):
        if channel_bins[channel] == 1:  ##counting
            for idxs, sample in enumerate(spec["channels"][idxc]["samples"]):
                DC.exp[channel].update(
                    {
                        sample["name"]: spec[
                            "channels"
                        ][idxc]["samples"][idxs]["data"]
                    }
                )
        else:  ##shapes
            DC.shapeMap.update({channel: {}})
            DC.hasShapes = True
            for idxs, sample in enumerate(spec["channels"][idxc]["samples"]):
                DC.exp[channel].update(
                    {sample["name"]: -1}
                )
                mods = [
                    spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["type"]
                    for i, mod in enumerate(spec["channels"][idxc]["samples"][idxs]["modifiers"])
                ]
                if "histosys" in mods:
                    DC.shapeMap[channel].update(
                        {
                            spec["channels"][idxc]["samples"][idxs]["name"]: [
                                options.shapefile,
                                channel+
                                "/" + spec["channels"][idxc]["samples"][idxs]["name"],
                                channel+
                                "/" + spec["channels"][idxc]["samples"][idxs]["name"]
                                + "_$SYSTEMATIC",
                            ]
                        }    
                    )
                else:
                    DC.shapeMap[channel].update(
                        {
                            spec["channels"][idxc]["samples"][idxs]["name"]: [
                                options.shapeFile,
                                channel+
                                "/" + spec["channels"][idxc]["samples"][idxs]["name"]
                                ]
                        }
                    )
                if "staterror" in mods:
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel], 0, channel_bins[channel]
                    ).Weight()
                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            spec["channels"][idxc]["samples"][idxs]["modifiers"][
                                mods.index("staterror")
                            ]["data"],
                        ],
                        axis=-1,
                    )
                    DC.binParFlags.update({channel : True})
                elif "shapesys" in mods:
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel], 0, channel_bins[channel]
                    ).Weight()
                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            spec["channels"][idxc]["samples"][idxs]["modifiers"][
                                mods.index("shapesys")
                            ]["data"],
                        ],
                        axis=-1,
                    )
                    DC.binParFlags.update({channel : True})
                else:
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel], 0, channel_bins[channel]
                    ).Weight()
                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            [0 for i in range(channel_bins[channel])],
                        ],
                        axis=-1,
                    )
                file[channel+  "/" + spec["channels"][idxc]["samples"][idxs]["name"]] = h_data


def addMods():
    for im, modifier in enumerate(systs):
        if "normsys" in modifier[1]:
            DC.systs.append((modifier[0], False, "lnN", [], {}))
            for idxc, channel in enumerate(channels):
                DC.systs[im][4].update({channel: {}})  
                for sample in samples:
                    if sample in [
                        spec["channels"][idxc]["samples"][i]["name"]
                        for i, s in enumerate(spec["channels"][idxc]["samples"])
                    ]:
                        for i in DC.systs:
                            i[4][channel].update({sample: 0.0})
        if "histosys" in modifier[1]:
            DC.systs.append((modifier[0], False, "shape", [], {}))
            for idxc, channel in enumerate(channels):
                DC.systs[im][4].update({channel: {}}) 
                for sample in samples:
                    if sample in [
                        spec["channels"][idxc]["samples"][i]["name"]
                        for i, s in enumerate(spec["channels"][idxc]["samples"])
                    ]:
                        for i in DC.systs:
                            i[4][channel].update({sample: 0.0})       
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(
            [
                spec["channels"][idxc]["samples"][i]["name"]
                for i, s in enumerate(spec["channels"][idxc]["samples"])
            ]
        ):
            mods = []
            for i, mod in enumerate(spec["channels"][idxc]["samples"][idxs]["modifiers"]):
                mods.append(
                    [
                        spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["name"],
                        spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["type"],
                        spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["data"],
                    ]
                )
            names = []
            for mod in mods:
                names.append(mod[0])
            for syst in DC.systs:
                name = syst[0]
                type = syst[2]
                if name in names:  ##if systematic is a modifier for this sample
                    if "lnN" in type:
                        syst[4][channel].update(
                            {sample: mods[names.index(name)][2]["lo"] + "/" + mods[names.index(name)][2]["high"]} ##asymmetric lnN
                        )
                    if "shape" in type:
                        syst[4][channel].update({sample: 1.0})
                        hi_data = hist.Hist.new.Regular(
                        channel_bins[channel], 0, channel_bins[channel]
                        ).Weight()
                        hi_data[...] = np.stack(
                        [
                            mods[names.index(name)][2]["hi_data"],
                            [0 for i in range(channel_bins[channel])],
                        ],
                        axis=-1,
                        )
                        lo_data = hist.Hist.new.Regular(
                        channel_bins[channel], 0, channel_bins[channel]
                        ).Weight()
                        lo_data[...] = np.stack(
                        [
                            mods[names.index(name)][2]["lo_data"],
                            [0 for i in range(channel_bins[channel])],
                        ],
                        axis=-1,
                        )
                        file[channel+  "/" + spec["channels"][idxc]["samples"][idxs]["name"] + "_" + name + "Up"] = hi_data
                        file[channel+  "/" + spec["channels"][idxc]["samples"][idxs]["name"] + "_" + name + "Down"] = lo_data
def addSignal():
    signalMods = []
    for modifier in modifiers:
        if modifier[0] == workspace.get_measurement()['config']['poi']:
            signalMods.append(modifier[0])
            
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(spec["channels"][idxc]["samples"]):
            for i, mod in enumerate(spec["channels"][idxc]["samples"][idxs]["modifiers"]):
                for sig in signalMods:
                    if sig in mod["name"]:
                       
                        DC.isSignal.update({sample["name"]: True})
                        DC.signals.append(sample["name"])

def addRateParams():   
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(
            (spec["channels"][idxc]["samples"])
        ):
            for i, mod in enumerate(spec["channels"][idxc]["samples"][idxs]["modifiers"]): 
                if "normfactor" in mod["type"] or spec["channels"][idxc]["parameters"] != None:
                    if mod["name"] not in workspace.get_measurement()['config']['poi']:
                        DC.rateParams.update({channel + "AND" + sample["name"]: []})
            for idxl, lumi in spec["channels"][idxc]["parameters"]:
                    DC.rateParams[channel + "AND" + sample["name"]].append([[lumi["name"], 1, 0], ''])
            for i, mod in enumerate(spec["channels"][idxc]["samples"][idxs]["modifiers"]):
                if "normfactor" in mod["type"]:
                    if mod["name"] not in workspace.get_measurement()['config']['poi']:
                        DC.rateParams[channel + "AND" + sample["name"]].append([[mod["name"], 1, 0], ''])
                
            

def writeDataCard(path):
    with open(path, 'w') as f:
        f.write("imax " + str(size(DC.bins)) + "\n") 
        f.write("jmax " + str(size(DC.processes)- size(DC.signals))+ "\n") 
        f.write("kmax " + str(size(DC.systs, 0))+ "\n")  
        if DC.hasShapes:
            for channel in DC.shapeMap.keys():
                for sample in DC.shapeMap[channel].keys():
                    f.write("shape " + sample + " " + channel + " " + DC.shapeMap[channel][sample][0] + " " + DC.shapeMap[channel][sample][1] + " " + DC.shapeMap[channel][sample][2] + "\n")    
        f.write('\n---------------------------------\n') 
        f.write("bin ")
        for bin in DC.bins:
            f.write(bin + " ")
        f.write("\n")
        f.write("observation ")
        for channel in DC.obs.keys():
            f.write(str(DC.obs[channel]) + " ")
        f.write('\n---------------------------------\n')
        f.write("bin     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():
                f.write(channel + "  ")
        f.write("\n")
        f.write("process     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():
                f.write(sample + "  ")
        f.write("\n")
        f.write("process     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():
                if sample in DC.signals:
                    f.write(str(0) + "  ")
                else:
                    f.write(str(DC.processes.index(sample) + 1) + "  ")
        f.write("\n")
        f.write("rate     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():
                f.write(str(DC.exp[channel][sample]) + "  ")
        f.write('\n---------------------------------\n')  
        for syst in DC.systs:
            f.write(syst[0] + "  " + syst[2] + "  ")
            for i, rate in enumerate(syst[4]):
                if rate != 0:
                    f.write(str(rate) + "  ")
                else:
                    f.write('-  ')
            f.write("\n")
        f.write('\n---------------------------------\n')  
        for cAp in DC.rateParams.keys():
            dir = cAp.split("AND")
            for i in range(size(DC.rateParams[cAp], 0)):
                f.write(str(DC.rateParams[cAp][i][0][0]) + " " + "rateParam " + dir[0] + " " + dir[1] + " " + str(DC.rateParams[cAp][i][0][1]))
                f.write("\n")
        f.write('\n---------------------------------\n') 
        for idxc, channel in enumerate(channels):
            if channel in DC.binParFlags.keys():
                if DC.binParFlags[channel] == True: ##double check to be safe
                    includeSignal = False
                    shapesys = False
                    staterror  = False
                    for idxs, sample in enumerate(spec["channels"][idxc]["samples"]):
                        if "shapesys" in [mod["type"] for i, mod in enumerate(sample["modifiers"])]:
                            if sample in DC.signals:
                                includeSignal = True
                            shapesys = True
                        elif "staterror" in [mod["type"] for i, mod in enumerate(sample["modifiers"])]:
                            if sample in DC.signals:
                                includeSignal = True
                            staterror = True
                    if shapesys == True:
                        if includeSignal == True:
                            f.write(channel + " autoMCStats " + str(100000) + " " + str(1) + " " + str(1))
                        else:
                            f.write(channel + " autoMCStats " + str(100000) + " " + str(0) + " " + str(1))
                    if staterror == True:
                        if includeSignal == True:
                            f.write(channel + " autoMCStats " + str(0) + " " + str(1) + " " + str(1))
                        else:
                            f.write(channel + " autoMCStats " +str(0) + " " + str(0) + " " + str(1))
        f.close()


addChannels()
addSamples()
addMods()
addSignal()
addRateParams()
file.close()
writeDataCard(options.outdatacard)











                    



            
                    
                    






