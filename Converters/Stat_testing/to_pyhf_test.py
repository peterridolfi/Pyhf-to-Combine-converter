
import pyhf
import json
import numpy as np


with open("shape_test.json") as file:
    spec = json.load(file)

ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
poi_values = np.linspace(0, 10)
##obs_limit, exp_limits, (scan, results) = pyhf.infer.intervals.upperlimit(
##    observations, model, poi_values, level=0.05, return_results=True
##)
bestFit = pyhf.infer.mle.fit(data = observations, pdf = model)
print(model.config.par_order)
print(bestFit)
