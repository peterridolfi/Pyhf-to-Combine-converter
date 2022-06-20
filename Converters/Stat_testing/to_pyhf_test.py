
import pyhf
import json
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import uproot


with open("shape_test.json") as file:
    spec = json.load(file)

ws = pyhf.Workspace(spec)
model = ws.model()

fit = uproot.open("fitDiagnosticsTest_shape.root")
total = fit['shapes_prefit/bin1/total'].values().tolist()

plt.hist([i for i in range(10)], bins = 10, range = [0, 10], weights = np.array(total).flatten())
plt.savefig('fitDiagnostics_counting')
'''observations = ws.data(model)
poi_values = np.linspace(0, 10)
obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upperlimit(
    observations, model, poi_values, level=0.05, return_results=True
)
init_pars = model.config.suggested_init()
print(init_pars)
exp = model.expected_actualdata(init_pars)
plt.hist([i for i in range(10)], bins = 10, range = [0,10], weights = exp)
plt.savefig('exp_data_counting')'''


"""fig, ax = plt.subplots()
fig.set_size_inches(10.5, 7)
ax.set_title("Hypothesis Tests")

artists = brazil.plot_results(poi_values, results, ax=ax)
plt.savefig('shapes_brazil_band')
"""





