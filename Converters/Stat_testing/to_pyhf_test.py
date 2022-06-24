
from array import array
from pyhf import workspace
import pyhf
import json
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import uproot
import random




fit = uproot.open("fitDiagnosticsTest.root"), 
b1total = [337.7, 340.3, 343.2, 346.2, 349.5, 353, 356.8, 360.9, 365.2, 369.9]
plt.plot(np.linspace(-2, 1.6, 10).tolist(), b1total, label = "combine yields", alpha = 0.5)

'''
file = uproot.open('higgsCombineTest.MultiDimFit.mH120.root')
tree = file['limit']
branches = tree.arrays()
plt.plot([i for i in range(51)], branches['nll'])
file.close()
    


with open("converted_workspace.json") as file:
    spec = json.load(file)

ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)

poi_values = np.linspace(0, 50, 50).tolist()
NLL = []
for poi in poi_values:
    NLL.append(model.logpdf(pars=[poi, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], data=observations))
plt.plot(poi_values, NLL)


plt.savefig('NLL_plot')'''



with open('converted_workspace.json') as file:
    spec = json.load(file)
ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
pars = np.linspace(-2, 1.6, 10).tolist()
exp = []
for p in pars:
    exp.append(model.expected_actualdata([1.5, p]))
print(np.linspace(-2, 1.6, 10).tolist())
plt.plot(np.linspace(-2, 1.6, 10).tolist(), exp, label = "pyhf yields", alpha = 0.5)
plt.legend()
plt.savefig('exp_data_counting')


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





