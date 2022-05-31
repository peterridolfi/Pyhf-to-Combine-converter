import json
import pyhf
import numpy as np

with open("C:/Users/peter/iris_project/shape_example_ws.json") as serialized:
    spec = json.load(serialized)
workspace = pyhf.Workspace(spec)
model = workspace.model()
observations = workspace.data(model, include_auxdata=True)
poi_values = np.linspace(0, 10, 10)
obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upperlimit(observations, model, poi_values, level = 0.05, return_results= True)
print(obs_limit, exp_limits)