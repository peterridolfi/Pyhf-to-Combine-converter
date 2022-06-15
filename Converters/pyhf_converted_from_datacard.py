from array import array
from email.policy import default
from modulefinder import Module
import string
from hist import Hist
import pyhf
import uproot
import json
import hist
import numpy
import HiggsAnalysis.CombinedLimit.DatacardParser as DP

import sys
from sys import exit
from optparse import OptionParser


from HiggsAnalysis.CombinedLimit.Datacard import Datacard

parser = OptionParser()
DP.addDatacardParserOptions(parser)
parser.add_option("-O", "--out-file", dest = "outfile", default = "converted_workspace.json" )
options, args = parser.parse_args()


DC = Datacard()
with open(args[0]) as dc_file:
    DC = Datacard()
    DC = DP.parseCard(file=dc_file, options=options)
    
channels = [channel for channel in DC.bins]
observations = [obs for channel, obs in DC.obs.items()]
samples = [sample for sample in DC.processes]
exp_values = DC.exp
sig = DC.isSignal
mods = DC.systs



   

def getShapeFile(shapeMap: dict, channel, sample)->string:
    file = ''
    if not channel in shapeMap.keys():
        if '*' in shapeMap.keys():
            if not sample in shapeMap['*']:
                if '*' in shapeMap['*'].keys():
                    file = shapeMap["*"]["*"][0]
            else:
                file = shapeMap["*"][sample][0]
    else:
        file = shapeMap[channel][sample][0]
    return file

def getHistPath(shapeMap: dict, channel, sample)->string:
    path = ''
    if not channel in shapeMap.keys():
        if '*' in shapeMap.keys():
            if not sample in shapeMap['*']:
                if '*' in shapeMap['*'].keys():
                    path = shapeMap["*"]["*"][1].replace("$CHANNEL", channel).replace("$PROCESS", sample)
            else:
                path = shapeMap["*"][sample][1].replace("$CHANNEL", channel)
    else:
        path = shapeMap[channel][sample][1]
    return path

def getUncertPath(shapeMap: dict, channel, sample, name)->string:
    path = ''
    if not channel in shapeMap.keys():
        if '*' in shapeMap.keys():
            if not sample in shapeMap['*']:
                if '*' in shapeMap['*'].keys():
                    path = shapeMap["*"]["*"][2].replace("$CHANNEL", channel).replace("$PROCESS", sample).replace("$SYSTEMATIC", name)
            else:
                path = shapeMap["*"][sample][2].replace("$CHANNEL", channel).replace("SYSTEMATIC", name)
    else:
        path = shapeMap[channel][sample][2].replace("$SYSTEMATIC", name)
    return path

def getHist(shapeMap: dict, channel, sample):
    file = uproot.open(getShapeFile(shapeMap, channel, sample))
    hist = file[getHistPath(shapeMap, channel, sample)]
    return hist
def getUncertUp(shapeMap: dict, channel, sample, name):
    file = uproot.open(getShapeFile(shapeMap, channel, sample))
    hist = file[getUncertPath(shapeMap, channel, sample, name)+"Up"]
    return hist
def getUncertDown(shapeMap: dict, channel, sample, name):
    file = uproot.open(getShapeFile(shapeMap, channel, sample))
    hist = file[getUncertPath(shapeMap, channel, sample, name)+"Down"]
    return hist


def addChannels(spec: dict):
    if DC.hasShapes:
        for idxc, channel in enumerate(channels):
            spec["channels"].append({"name": channel, "samples": []})
            hist = getHist(DC.shapeMap, channel, "data_obs")
            data = hist.values().tolist()
            spec["observations"].append({"name": channel, "data": data})
    else:
        for idxc, channel in enumerate(channels):
            spec["channels"].append({"name": channel, "samples": []})
            spec["observations"].append({"name": channel, "data": [observations[idxc]]})

def addSamples(spec: dict):
    if DC.hasShapes:
        for idxc, channel in enumerate(channels):
            for idxs, sample in enumerate(samples):
                hist = getHist(DC.shapeMap, channel, sample)
                data = hist.values().tolist()
                if sample in exp_values[channel].keys():
                    spec["channels"][idxc]["samples"].append(
                    {
                        "name": sample,
                        "data": data,
                        "modifiers": [],
                    }
                )
    else:    
        for idxc, channel in enumerate(channels):
            for idxs, sample in enumerate(samples):
                if sample in exp_values[channel].keys():
                    spec["channels"][idxc]["samples"].append(
                        {
                            "name": sample,
                            "data": [exp_values[channel][sample]],
                            "modifiers": [],
                        }
                    )


def addMeasurements(spec: dict):
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(samples):
            if sig[sample] == True:
                if (channel + "AND" + sample) in DC.rateParams.keys():
                    name = DC.rateParams[channel + "AND" + sample][0][0][0]
                    spec["measurements"].append(
                        {"name": "Measurement", "config": {"poi": name, "parameters": []}}
                    )
                else:
                    spec["measurements"].append(
                        {"name": "Measurement", "config": {"poi": "mu", "parameters": []}}
                    )


def addNormFactor(spec: dict):
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(samples):
            if (channel + "AND" + sample) in DC.rateParams.keys():
                    name = DC.rateParams[channel + "AND" + sample][0][0][0]
                    spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                        {"name": name, "type": "normfactor", "data": None}
                    )
            elif sig[sample] == True:
                spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                        {"name": "mu", "type": "normfactor", "data": None}
                    )


def addMods(spec: dict):
    for syst in mods:
        name = syst[0]
        mod_type = syst[2]
        if mod_type == "lnN":  ##normsys
            for idxc, channel in enumerate(channels):
                for idxs, sample in enumerate(exp_values[channel].keys()):
                    if sample in exp_values[channel].keys():
                        if syst[4][channel][sample] != 0:
                            spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                                {
                                    "name": name,
                                    "type": "normsys",
                                    "data": {
                                        "hi": syst[4][channel][sample],
                                        "lo": 1 / (syst[4][channel][sample]),
                                    },
                                }
                            )
        if "shape" in mod_type: ##histosys
            for idxc, channel in enumerate(channels):
                for idxs, sample in enumerate(samples):
                    if syst[4][channel][sample] != 0:    
                        if sample in exp_values[channel].keys():    
                            histUp = getUncertUp(DC.shapeMap, channel, sample, name)
                            histDown = getUncertDown(DC.shapeMap, channel, sample, name)
                            hi_data = histUp.values().tolist()
                            lo_data = histDown.values().tolist()
                            spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                            {
                                "name": name,
                                "type": "histosys",
                                "data": {"hi_data": [count / syst[4][channel][sample] for count in hi_data], "lo_data": [count / syst[4][channel][sample] for count in lo_data]},
                            }
                        )
    for idxc, channel in enumerate(channels):  ##staterror/shapesys

        if channel in DC.binParFlags.keys():
            if DC.binParFlags[channel][0] == 0: #if BB lite enabled               
                for idxs, sample in enumerate(samples):
                    if sample in exp_values[channel].keys():
                        hist = getHist(DC.shapeMap, channel, sample)
                        err = hist.errors().tolist()
                        spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                        {
                            "name": "my_stat_err",
                            "type": "staterror",
                            "data": err
                        }
                    )

            elif DC.binParFlags[channel][0] > 0: ##if user wants total BB (shapesys)               
                for idxs, sample in enumerate(samples):
                    if sample in exp_values[channel].keys():
                        hist = getHist(DC.shapeMap, channel, sample)
                        err = hist.errors().tolist()
                        spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                        {
                            "name": "my_shapesys_" + channel + sample,
                            "type": "shapesys",
                            "data": err
                        }
                    )
                

            
                
def toJSON(spec: dict):
    addChannels(spec)
    addSamples(spec)
    addMeasurements(spec)
    addNormFactor(spec)
    addMods(spec)

def writeFileName(name):
    with open(name, "w") as file:
        file.write(json.dumps(spec, indent=2))
    
        
spec = {"channels": [], "observations": [], "measurements": [], "version": "1.0.0"}
toJSON(spec)

writeFileName(options.outfile)





