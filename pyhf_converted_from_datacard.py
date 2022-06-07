import CMSSW_11_2_0.pyhf as pyhf

import json
import CMSSW_11_2_0.hist as hist
import CMSSW_11_2_0.python.HiggsAnalysis.CombinedLimit.DatacardParser as DP

from sys import exit
from optparse import OptionParser

from CMSSW_11_2_0.python.HiggsAnalysis.CombinedLimit.Datacard import Datacard
parser = OptionParser()
DP.addDatacardParserOptions(parser)
options,args = parser.parse_args()

dc_file=  open("/home/cmsusr/pyhf_converter/CMSSW_11_2_0/src/HiggsAnalysis/CombinedLimit/data/tutorials/counting/realistic-multi-channel.txt")
DC = Datacard()
DC = DP.parseCard(file=dc_file, options=options)
dc_file.close()

##get arrays of datacard info
channels = [channel for channel in DC.bins]
observations = [obs for channel, obs in DC.obs.items()]
samples = [sample for sample in DC.processes]
exp_values = DC.exp
sig = DC.isSignal
mods = DC.systs


spec = {

  "channels": [
    
  ],
  "observations": [
    
  ],
  "measurements": [
  
  ],
  "version": "1.0.0"
}

for idxc, channel in enumerate(channels):
  spec["channels"].append({"name": channel, "samples": []})
  spec["observations"].append({"name": channel, "data": observations[idxc]})
  for idxs, sample in enumerate(samples):
    if sample in exp_values[channel].keys():
      spec["channels"][idxc]["samples"].append({"name": sample, "data": exp_values[channel][sample], "modifiers": []})    
for idx, sample in enumerate(samples):
  if sig[sample]==True:
    spec["measurements"].append({"name": "Measurement", "config": {"poi": "mu", "parameters": []}})  
    
for idxc, channel in enumerate(channels):
  for idxs, sample in enumerate(samples):
    if sig[sample] == True:
      spec["channels"][idxc]["samples"][idxs]["modifiers"].append({"name": "mu", "type": "normfactor", "data": None})    
for syst in mods:
  name = syst[0]
  mod_type = syst[2]
  if mod_type == "lnN":
    for idxc, channel in enumerate(channels):
      for idxs, sample in enumerate(exp_values[channel].keys()):
        if sample in exp_values[channel].keys():
          if syst[4][channel][sample]!=0:
            spec["channels"][idxc]["samples"][idxs]["modifiers"].append({"name": name, "type": mod_type, "data": [syst[4][channel][sample], 1/(syst[4][channel][sample])]})







print(json.dumps(spec, indent=2))










