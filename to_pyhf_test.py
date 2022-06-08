import pyhf
import json
import numpy as np


with open("converted_workspace.json") as file:
    spec = json.load(file)

ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
poi_values = np.linspace(1, 10, 1)
obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upperlimit(
    observations, model, poi_values, level=0.05, return_results=True
)
print(obs_limit, exp_limits)
