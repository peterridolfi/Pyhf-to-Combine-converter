
import pyhf
import json
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil


with open("./Stat_testing/shape_test.json") as file:
    spec = json.load(file)

ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
poi_values = np.linspace(0, 10)
obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upperlimit(
    observations, model, poi_values, level=0.05, return_results=True
)


fig, ax = plt.subplots()
fig.set_size_inches(10.5, 7)
ax.set_title("Hypothesis Tests")

artists = brazil.plot_results(poi_values, results, ax=ax)
plt.savefig('shapes_brazil_band')





