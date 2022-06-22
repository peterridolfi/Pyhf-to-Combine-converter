
from pyhf import workspace
import pyhf
import json
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import uproot




'''fit = uproot.open("fitDiagnosticsTest.root")
emutotal = fit['shapes_prefit/e_mu/total'].values().tolist()
etautotal = fit['shapes_prefit/e_tau/total'].values().tolist()
mutautotal = fit['shapes_prefit/mu_tau/total'].values().tolist()
plt.hist([i for i in range(3)], bins = 3, range = [0, 3], weights = np.array([emutotal, etautotal, mutautotal]).flatten())
'''

with open("counting_test.json") as file:
    spec = json.load(file)

ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)

poi_values = np.linspace(0, 10, 50).tolist()
NLL = []
for poi in poi_values:
    NLL.append(model.logpdf(pars=[poi, 1, 1, 1, 1, 1, 1, 1], data=observations))
plt.plot(poi_values, NLL)

file = uproot.open('higgsCombineTest.MultiDimFit.mH120.root')

tree = file['limit']
branches  =tree.arrays()
print(poi_values, branches['syst'])

plt.savefig('NLL_plot')

'''observations = ws.data(model)
pars = [ 6.38810694e+00,  4.10779419e-02, -5.09961047e-01,  1.56291969e-02,
  4.19327212e-04,  3.01550046e-01, -2.64537207e-01,  3.68338607e-02]
exp = model.expected_actualdata(pars)
plt.hist([i for i in range(3)], bins = 3, range = [0,3], weights = exp)
plt.savefig('exp_data_counting_varyparams')'''


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





