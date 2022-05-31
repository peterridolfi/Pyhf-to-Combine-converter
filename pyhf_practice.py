import pyhf
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import json

workspace = pyhf.Workspace({
  "channels": [
    {
      "name": "singlechannel",
      "samples": [
        {
          "name": "signal",
          "data": [
            5.0, 10.0, 15.0
          ],
          "modifiers": [
            {
              "name": "mu",
              "type": "normfactor",
              "data": None
            }
          ]
        },
        {
          "name": "background",
          "data": [
            50.0, 70.0, 80.0
          ],
          "modifiers": [
            {
              "name": "uncorr_bkguncrt",
              "type": "shapesys",
              "data": [
                5.0, 10.0, 10.0
              ]
            },
            {
                "name": "corr_shape",
                "type": "histosys",
                "data": {"hi_data": [20, 15, 15], "lo_data": [10, 5, 10]}
            }
          ]
        }
      ]
    }
  ],
  "observations": [
    {
      "name": "singlechannel",
      "data": [
        53.0, 65.0, 70.0
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
})

model = workspace.model()
poi_values = np.linspace(1, 5 , 5)
observations = workspace.data(model, include_auxdata = True)




"""
#calculating test statistic by hand
mle_fixed = [pyhf.infer.mle.fixed_poi_fit(poi_value, data=observations, pdf=model, return_fitted_val=True) for poi_value in poi_values]
mle_free = [pyhf.infer.mle.fit(data=observations, pdf=model, return_fitted_val=True) for poi_value in poi_values]
numer = [i[1] for i in mle_fixed]
denom = [i[1]for i in mle_free]
test_stat = [0 for i in numer]
for i in range (len(numer)):
    test_stat[i] = float(numer[i]) / float(denom[i])
plt.plot(poi_values, test_stat)
plt.title("test statistic against mu")
plt.xlabel("mu")
plt.ylabel("test stat")
plt.show()
"""

"""graphing CLs by hand
CLs_set = [pyhf.infer.hypotest(poi_value, observations, model, test_stat = "qtilde", return_expected_set=True) for poi_value in poi_values]
CLs_obs = np.asarray([i[0] for i in CLs_set])
CLs_exp = np.asarray([i[1][2] for i in CLs_set])
alpha_level = [0.05 for i in CLs_set]
plt.plot(poi_values, CLs_obs, label="CLs_obs")
plt.plot(poi_values, CLs_exp, linestyle='dashed', label="CLs_exp")
plt.plot(poi_values, alpha_level, 'r--', label="alpha_level")
plt.title("CLs graph against mu")
plt.xlabel("mu")
plt.ylabel("CLs")
plt.legend()
plt.show()
"""

"""
#graph of p-values with varying mu
CLs_set = [pyhf.infer.hypotest(poi_value, observations, model, test_stat = "qtilde", return_tail_probs = True) for poi_value in poi_values]
p_value_sb = [i[1][0] for i in CLs_set]
p_value_b = [i[1][1] for i in CLs_set]
CLs_value = [0 for i in range (len(p_value_b))]
for i in range (len(p_value_b)):
  CLs_value[i] = float(p_value_sb[i])/(1-float(p_value_b[i]))
alpha = [0.05 for i in CLs_set]
plt.plot(poi_values, p_value_sb, 'g', label= "p_value_sb")
plt.plot(poi_values, p_value_b, 'b', label="p_value_b")
plt.plot(poi_values, CLs_value, 'y', label = "CLs_value")
plt.plot(poi_values, alpha, 'r--', label="alpha_level")
plt.title("p-values against mu")
plt.legend()
plt.xlabel("mu")
plt.ylabel("p-value")
plt.show()
"""
"""
#histogram of data using pyhf API
bins = [idx for idx, i in enumerate(data)]
exp_data = [model.expected_data([poi_value, model.config.suggested_init()[1], model.config.suggested_init()[2], model.config.suggested_init()[3], model.config.suggested_init()[4]], include_auxdata=False) for poi_value in poi_values]
fig, ax = plt.subplots(1, 5)
for i in range(len(obs_data)):
  ax[i].scatter(bins, data, 1.0, label="real data", color = "black" ) 
for i in range(len(exp_data)):
  ax[i].bar(bins, exp_data[i], 1.0, label="exp data", edgecolor = "green", alpha = 0.5 )
  ax[i].legend()
  ax[i].set_ylim(50, 150)
plt.show()
"""

