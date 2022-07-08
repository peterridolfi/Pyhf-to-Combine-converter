
from array import array
from pyhf import workspace
import pyhf
import json
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import uproot
import math






with open("converted_workspace.json") as file:
    spec = json.load(file)
ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
poi_values = np.linspace(0, 3, 50).tolist()
NLL = []
for poi in poi_values:
    NLL.append(-1*(model.logpdf(pars=[poi, 1], data=observations))-35.91)

plt.plot(poi_values, NLL, label = "pyhf NLL")
file = uproot.open('higgsCombineTest.MultiDimFit.mH120.root')
tree = file['limit']
branches = tree.arrays()
plt.plot(np.linspace(0, 3, 50), branches['deltaNLL'][1:], label  = "combine NLL")
file.close()
plt.legend()
plt.savefig('NLLplot')
    





'''
fit = uproot.open("fitDiagnosticsTest.root")
hist = fit["shapes_prefit/b1/total"].values().tolist()
plt.hist(np.linspace(0, 10, 10), bins=10, weights = hist, range=(0, 10), label = "combine yields", alpha = 0.5)

with open('converted_workspace.json') as file:
    spec = json.load(file)
ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
exp = model.expected_actualdata([1, 1.5])
plt.hist(np.linspace(0, 10, 10), bins=10, weights = exp, range=(0, 10), label = "pyhf yields", alpha = 0.5)
plt.legend()
plt.savefig('exp_data_counting')

print(pyhf.infer.mle.fit(data = observations, pdf = model))'''


'''with open("counting_test.json") as file:
    spec = json.load(file)

ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
print(pyhf.infer.mle.fit(data = observations, pdf = model))
print(model.config.par_order)'''


"""fig, ax = plt.subplots()
fig.set_size_inches(10.5, 7)
ax.set_title("Hypothesis Tests")

artists = brazil.plot_results(poi_values, results, ax=ax)
plt.savefig('shapes_brazil_band')
"""





