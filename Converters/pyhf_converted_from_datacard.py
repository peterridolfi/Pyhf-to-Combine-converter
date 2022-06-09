from email.policy import default
from modulefinder import Module
import string
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
   





def addChannels(spec: dict):
    for idxc, channel in enumerate(channels):
        spec["channels"].append({"name": channel, "samples": []})
        spec["observations"].append({"name": channel, "data": [observations[idxc]]})


def addSamples(spec: dict):
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
    for idxs, sample in enumerate(samples):
        if sig[sample] == True:
            spec["measurements"].append(
                {"name": "Measurement", "config": {"poi": "mu", "parameters": []}}
            )


def addNormFactor(spec: dict):
    for idxc, channel in enumerate(channels):
        for idxs, sample in enumerate(samples):
            if sig[sample] == True:
                spec["channels"][idxc]["samples"][idxs]["modifiers"].append(
                    {"name": "mu", "type": "normfactor", "data": None}
                )


def addMods(spec: dict):
    for syst in mods:
        name = syst[0]
        mod_type = syst[2]
        if mod_type == "lnN":
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


def toJSON(spec: dict):
    addChannels(spec)
    addSamples(spec)
    addMeasurements(spec)
    addNormFactor(spec)
    addMods(spec)
    
        
spec = {"channels": [], "observations": [], "measurements": [], "version": "1.0.0"}


def writeFileName(name):
    with open(name, "w") as file:
        file.write(json.dumps(spec, indent=2))



toJSON(spec)
writeFileName(options.outfile)




