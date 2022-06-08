import uproot
import pyhf
import matplotlib.pyplot as plt
import numpy as np
import awkward as ak
import vector

file = uproot.open("C:/Users/peter/CMS_combine/data/tutorials/shapes/simple-shapes-TH1_input.root")

bkg = file["background"].to_numpy()[0].tolist()
bkg_up = file["background_alphaUp"].to_numpy()[0].tolist()
bkg_down = file["background_alphaDown"].to_numpy()[0].tolist()
sig = file["signal;1"].to_numpy()[0].tolist()
sig_up = file["signal_sigmaUp"].to_numpy()[0].tolist()
sig_down = file["signal_sigmaDown"].to_numpy()[0].tolist()
data = file["data_obs"].to_numpy()[0].tolist()


spec = {
    "channels": [
        {
            "name": "bin1",
            "samples": [
                {
                    "name": "signal",
                    "data": sig,
                    "modifiers": [
                        {
                            "name": "mu",
                            "type": "normfactor",
                            "data": None
                        },
                        {
                            "name": "sigma",
                            "type": "histosys",
                            "data": {"hi_data": sig_up, "lo_data": sig_down}
                        },
                        {
                            "name": "lumi",
                            "type": "normsys",
                            "data": {"hi": 1.1, "lo": 0.91}
                        }
                    ]
                },
                {
                    "name": "bkg",
                    "data": bkg,
                    "modifiers": [
                        {
                            "name": "sigma",
                            "type": "histosys",
                            "data": {"hi_data": bkg_up, "lo_data": bkg_down}
                        },
                        {
                            "name": "bg_norm",
                            "type": "normsys",
                            "data": {"hi": 1.3, "lo": 0.77}
                        }
                    ]
                }
                
            ]
        }
    ],
    "observations": [
        {
            "name": "bin1",
            "data": data
        }
    ],
    "measurements": [
    {
      "name": "Measurement",
      "config": {
        "poi": "mu",
        "parameters": []
      }
    }
  ],
  "version": "1.0.0" 
}

workspace = pyhf.Workspace(spec)
model=workspace.model()
observations = workspace.data(model, include_auxdata=True)
poi_values = np.linspace(0, 5, 10)
obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upperlimit(observations, model, poi_values, level = 0.05, return_results= True)
print(obs_limit, exp_limits)
