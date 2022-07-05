
from array import array
from pyhf import workspace
import pyhf
import json
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import uproot
import random




fit = uproot.open("fitDiagnosticsTest.root")
hist = fit["shapes_fit_s/b1/total"].values().tolist()
plt.hist(np.linspace(0, 10, 10), bins=10, weights = hist, range=(0, 10), label = "combine yields", alpha = 0.5)

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
data = [
            22.100000381469727,
            20.649999618530273,
            18.700000762939453,
            16.950000762939453,
            16.600000381469727,
            14.699999809265137,
            13.399999618530273,
            10.949999809265137,
            10.800000190734863,
            9.449999809265137
          ]
err = [
                1.051189802081432,
                1.016120071645079,
                0.9669539802906859,
                0.9205976319760986,
                0.91104335791443,
                0.8573214099741124,
                0.8185352771872451,
                0.7399324293474372,
                0.7348469228349536,
                0.687386354243376
              ]
percent = [err[i] / data[i] for i, num in enumerate(data)]
fit = [0.410499, -.18183, -.51804, -.23517, -.15169, .16846, -.25026, -.89728, -.21374, .15744, -.19908]

exp = model.expected_actualdata([0.410499]+[1+percent[i]*fit[i+1] for i, num in enumerate(percent)])


plt.hist(np.linspace(0, 10, 10), bins=10, weights = exp, range=(0, 10), label = "pyhf yields", alpha = 0.5)
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





