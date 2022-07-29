from ctypes import sizeof
from operator import indexOf
from unicodedata import name
import numpy as np
from numpy.core.fromnumeric import size
from pyhf import parameters
import pyhf
import uproot
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
    DC.bins = [channel["name"] for i, channel in enumerate(spec["channels"])]
    for idxc, channel in enumerate(spec["channels"]):
        DC.exp.update({channel["name"] : {}})
        if channel_bins[channel["name"]] == 1:  ##single bin
            DC.obs.update({channel["name"]: spec["observations"][idxc]["data"]})
            
        else:
            data = sum(spec["observations"][idxc]["data"])
            DC.obs.update({channel["name"]: data})
            h_data = hist.Hist.new.Regular(
                channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
            ).Weight()
            h_data[...] = np.stack(
                [
                    spec["observations"][idxc]["data"],
                    [0 for i in range(channel_bins[channel["name"]])],
                ],
                axis=-1,
            )
            
            file[channel["name"]+ "/data_obs"] = h_data


def addSamples():
    DC.processes = samples
    for idxc, channel in enumerate(spec["channels"]):
        if channel_bins[channel["name"]] == 1:  ##counting
            for idxs, sample in enumerate(channel["samples"]):
                DC.exp[channel["name"]].update(
                    {
                        sample["name"]: sample["data"]
                    }
                )
        else:  ##shapes
            DC.shapeMap.update({channel["name"]: {}})
            DC.shapeMap[channel["name"]].update(
                {
                    'data_obs': [
                        options.shapefile,
                        channel["name"]+"/"+'data_obs'
                    ]
                }
            )
            DC.hasShapes = True
            for idxs, sample in enumerate(channel["samples"]):
                data = sum(sample["data"])
                DC.exp[channel["name"]].update(
                    {sample["name"]: data}
                )
                
                mods = [
                    spec["channels"][idxc]["samples"][idxs]["modifiers"][i]["type"]
                    for i, mod in enumerate(spec["channels"][idxc]["samples"][idxs]["modifiers"])
                ]
                if "histosys" in mods:
                    DC.shapeMap[channel["name"]].update(
                        {
                            spec["channels"][idxc]["samples"][idxs]["name"]: [
                                options.shapefile,
                                channel["name"]+
                                "/" + spec["channels"][idxc]["samples"][idxs]["name"],
                                channel["name"]+
                                "/" + spec["channels"][idxc]["samples"][idxs]["name"]
                                + "_$SYSTEMATIC",
                            ]
                        }    
                    )
                else:
                    DC.shapeMap[channel["name"]].update(
                        {
                            spec["channels"][idxc]["samples"][idxs]["name"]: [
                                options.shapefile,
                                channel["name"]+
                                "/" + spec["channels"][idxc]["samples"][idxs]["name"]
                                ]
                        }
                    )
                if "staterror" in mods:
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
                    ).Weight()
                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            [spec["channels"][idxc]["samples"][idxs]["modifiers"][
                                mods.index("staterror")
                            ]["data"][i]**2 for i in spec["channels"][idxc]["samples"][idxs]["modifiers"][
                                mods.index("staterror")
                            ]["data"]],
                        ],
                        axis=-1,
                    )
                    DC.binParFlags.update({channel["name"] : True})
                elif "shapesys" in mods:
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
                    ).Weight()
                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            [spec["channels"][idxc]["samples"][idxs]["modifiers"][
                                mods.index("shapesys")
                            ]["data"][i]**2 for i in spec["channels"][idxc]["samples"][idxs]["modifiers"][
                                mods.index("shapesys")
                            ]["data"]],
                        ],
                        axis=-1,
                    )
                    DC.binParFlags.update({channel["name"] : True})
                else:
                    h_data = hist.Hist.new.Regular(
                        channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
                    ).Weight()
                    h_data[...] = np.stack(
                        [
                            spec["channels"][idxc]["samples"][idxs]["data"],
                            [0 for i in range(channel_bins[channel["name"]])],
                        ],
                        axis=-1,
                    )
                file[channel["name"]+  "/" + spec["channels"][idxc]["samples"][idxs]["name"]] = h_data


def addMods():
    for im, modifier in enumerate(systs):
        if "normsys" in modifier[1]:
            DC.systs.append((modifier[0], False, "lnN", [], {}))
            for idxc, channel in enumerate(spec["channels"]):
                DC.systs[im][4].update({channel["name"]: {}})  
                for idxs, sample in enumerate(channel["samples"]):
                    for i in DC.systs:
                            i[4][channel["name"]].update({sample["name"]: 0.0})
        if "histosys" in modifier[1]:
            DC.systs.append((modifier[0], False, "shape", [], {}))
            for idxc, channel in enumerate(spec["channels"]):
                DC.systs[im][4].update({channel["name"]: {}}) 
                for idxs, sample in enumerate(channel["samples"]):
                    for i in DC.systs:
                            i[4][channel["name"]].update({sample["name"]: 0.0})      
    for idxc, channel in enumerate(spec["channels"]):
        for idxs, sample in enumerate(channel["samples"]):
            mods = [sample["modifiers"][i] for i, m in enumerate(sample["modifiers"])]
            names = []
            for mod in mods:
                names.append(mod["name"])
            for syst in DC.systs:
                name = syst[0]
                type = syst[2]
                if name in names:  ##if systematic is a modifier for this sample
                    if "lnN" in type:
                        syst[4][channel["name"]].update(
                            {sample["name"]: str(mods[names.index(name)]["data"]["lo"]) + "/" + str(mods[names.index(name)]["data"]["hi"])} ##asymmetric lnN
                        )
                    if "shape" in type:
                        syst[4][channel["name"]].update({sample["name"]: 1.0})
                        hi_data = hist.Hist.new.Regular(
                        channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
                        ).Weight()
                        hi_data[...] = np.stack(
                        [
                            mods[names.index(name)]["data"]["hi_data"],
                            [0 for i in range(channel_bins[channel["name"]])],
                        ],
                        axis=-1,
                        )
                        lo_data = hist.Hist.new.Regular(
                        channel_bins[channel["name"]], 0, channel_bins[channel["name"]]
                        ).Weight()
                        lo_data[...] = np.stack(
                        [
                            mods[names.index(name)]["data"]["lo_data"],
                            [0 for i in range(channel_bins[channel["name"]])],
                        ],
                        axis=-1,
                        )
                        file[channel["name"]+  "/" + spec["channels"][idxc]["samples"][idxs]["name"] + "_" + name + "Up"] = hi_data
                        file[channel["name"]+  "/" + spec["channels"][idxc]["samples"][idxs]["name"] + "_" + name + "Down"] = lo_data
def addSignal():
    measurements = []
    for im, measurement in enumerate(spec["measurements"]):
        measurements.append(measurement["config"]["poi"])
    signalMods = []
    for modifier in modifiers:
        if modifier[0] in measurements:
            signalMods.append(modifier[0])       
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(spec["channels"][idxc]["samples"]):
            for i, mod in enumerate(spec["channels"][idxc]["samples"][idxs]["modifiers"]):
                for sig in signalMods:
                    if sig == mod["name"]:
                        DC.isSignal.update({sample["name"]: True})
                        DC.signals.append(sample["name"])

def addRateParams():
    measurements = []
    for im, measurement in enumerate(spec["measurements"]):
        measurements.append(measurement["config"]["poi"])
    signalMods = []
    for modifier in modifiers:
        if modifier[0] in measurements:
            signalMods.append(modifier[0])
        
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(
            (spec["channels"][idxc]["samples"])
        ):
            isSig = False
            for sig in signalMods:
                if sample["name"] in sig:
                    isSig = True
            if isSig == False:
                for i, mod in enumerate(spec["channels"][idxc]["samples"][idxs]["modifiers"]): ##normfactor or lumi 
                    if "normfactor" in mod["type"] or "lumi" in mod["type"] or "shapefactor" in mod["type"]:
                        DC.rateParams.update({channel + "AND" + sample["name"]: []})
                        DC.rateParams[channel + "AND" + sample["name"]].append([[mod["name"], 1, 0], ''])
                        
            
            
                
            

def writeDataCard(path):
    with open(path, 'w') as f:
        f.write("imax " + str(size(DC.bins)) + "\n") 
        f.write("jmax " + str(size(DC.processes)- size(DC.signals))+ "\n") 
        f.write("kmax " + str(size(DC.systs, 0))+ "\n")  
        if DC.hasShapes:  
            for channel in DC.shapeMap.keys():
                for sample in DC.shapeMap[channel].keys():
                    f.write("shapes " + sample + "  " + channel + "  " + DC.shapeMap[channel][sample][0] + "  " + DC.shapeMap[channel][sample][1])
                    if size(DC.shapeMap[channel][sample]) > 2:
                        f.write("  " + DC.shapeMap[channel][sample][2] + "\n")  
                    else:
                        f.write("\n")
        f.write('\n---------------------------------\n') 
        f.write("bin ")
        for bin in DC.obs.keys():
            f.write(bin + " ")
        f.write("\n")
        f.write("observation ")
        for channel in DC.obs.keys():
            f.write(str(DC.obs[channel]) + " ")
        f.write('\n---------------------------------\n')
        f.write("bin     ")
        for channel in DC.obs.keys():
            for sample in DC.exp[channel].keys():
                f.write(channel + "    ")
        f.write("\n")
        f.write("process     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():
                f.write(sample + "    ")
        f.write("\n")
        f.write("process     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():
                if sample in DC.signals:
                    f.write(str(-1*DC.processes.index(sample)) + "     ")
                else:
                    f.write(str(DC.processes.index(sample)) + "     ")
        f.write("\n")
        f.write("rate     ")
        for channel in DC.bins:
            for sample in DC.exp[channel].keys():
                
                f.write(str(DC.exp[channel][sample]) + "     ")
        f.write('\n---------------------------------\n') 
        for syst in DC.systs:
            f.write(syst[0] + "  " + syst[2] + "  ")
            for bin in syst[4].keys():
                for sample in DC.exp[bin].keys(): 
                    if syst[4][bin][sample] != 0:
                        f.write(str(syst[4][bin][sample]) + "  ")
                    else:
                        f.write('-       ')
                        
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











                    



            
                    
                    






