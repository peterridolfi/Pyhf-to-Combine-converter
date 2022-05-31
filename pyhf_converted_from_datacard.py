import pyhf
import numpy as np
import json

spec = {

  "channels": [
    {
      "name": "e_tau",
      "samples": [
        {
          "name": "higgs",
          "data": [
            0.34
          ],
          "modifiers": [
            {
              "name": "mu",
              "type": "normfactor",
              "data": None
            },
            {
              "name": "lumi",
              "type": "normsys",
              "data": {"hi": 1.11, "lo": 0.91}
            },
            {
              "name": "tauid",
              "type": "normsys",
              "data": {"hi": 1.23, "lo": 0.81}
            },
            {
              "name": "effic",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            }
          ]
        },
        {
          "name": "zTT",
          "data": [
            190
          ],
          "modifiers": [
            {
              "name": "tauid",
              "type": "normsys",
              "data": {"hi": 1.23, "lo": 0.81}
            },
            {
              "name": "ztoll",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            },
            {
              "name": "effic",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            }
          ]
        },
        {
             "name": "qcd",
          "data": [
            327
          ],
          "modifiers": [
            {
              "name": "qcdel",
              "type": "normsys",
              "data": {"hi": 1.20, "lo": 0.83}
            }
          ]
        }        
      ]
    },
    {
      "name": "mu_tau",
      "samples": [
        {
          "name": "higgs",
          "data": [
            0.57
          ],
          "modifiers": [
            {
              "name": "mu",
              "type": "normfactor",
              "data": None
            },
            {
              "name": "lumi",
              "type": "normsys",
              "data": {"hi": 1.11, "lo": 0.91}
            },
            {
              "name": "tauid",
              "type": "normsys",
              "data": {"hi": 1.23, "lo": 0.81}
            },
            {
              "name": "effic",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            }
          ]
        },
        {
          "name": "zTT",
          "data": [
            329
          ],
          "modifiers": [
            {
              "name": "tauid",
              "type": "normsys",
              "data": {"hi": 1.23, "lo": 0.81}
            },
            {
              "name": "ztoll",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            },
            {
              "name": "effic",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            }
          ]
        },
        {
             "name": "qcd",
          "data": [
            259
          ],
          "modifiers": [
            {
              "name": "qcdmu",
              "type": "normsys",
              "data": {"hi": 1.10, "lo": 0.91}
            }
          ]
        }        
      ]
    },
    {
      "name": "e_mu",
      "samples": [
        {
          "name": "higgs",
          "data": [
            0.15
          ],
          "modifiers": [
            {
              "name": "mu",
              "type": "normfactor",
              "data": None
            },
            {
              "name": "lumi",
              "type": "normsys",
              "data": {"hi": 0.91, "lo": 1.1}
            },
            {
              "name": "effic",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            }
          ]
        },
        {
          "name": "zTT",
          "data": [
            88
          ],
          "modifiers": [
            {
              "name": "ztoll",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            },
            {
              "name": "effic",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            }
          ]
        },
        {
             "name": "other",
          "data": [
            14
          ],
          "modifiers": [
            {
              "name": "lumi",
              "type": "normsys",
              "data": {"hi": 0.91, "lo": 1.1}
            },
            {
              "name": "effic",
              "type": "normsys",
              "data": {"hi": 1.04, "lo": 0.96}
            },
            {
              "name": "other",
              "type": "normsys",
              "data": {"hi": 1.1, "lo": 0.91}
            }
          ]
        }        
      ]
    }

    
  ],
  "observations": [
    {
      "name": "e_tau",
      "data": [
            517
      ]
    },
    {
        "name":"mu_tau",
        "data":[
            540
        ]
    },
    {
        "name":"e_mu",
        "data":[
            101
        ]
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
model = workspace.model()

observations = workspace.data(model, include_auxdata = True)
fit = pyhf.infer.mle.fit(data = observations, pdf = model)
print (fit)


