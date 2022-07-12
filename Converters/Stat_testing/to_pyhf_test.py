
from array import array
import scipy
from pyhf import workspace
import pyhf
import json
import numpy as np
import matplotlib.pyplot as plt
from pyhf.contrib.viz import brazil
import uproot
import math
from pyhf import get_backend
from pyhf.exceptions import UnspecifiedPOI

def twice_nll(pars, data, pdf):
    
    return -2 * pdf.logpdf(pars, data)

def _validate_fit_inputs(init_pars, par_bounds, fixed_params):
    for par_idx, (value, bound) in enumerate(zip(init_pars, par_bounds)):
        if not (bound[0] <= value <= bound[1]):
            raise ValueError(
                f"fit initialization parameter (index: {par_idx}, value: {value}) lies outside of its bounds: {bound}"
                + "\nTo correct this adjust the initialization parameter values in the model spec or those given"
                + "\nas arguments to pyhf.infer.fit. If this value is intended, adjust the range of the parameter"
                + "\nbounds."
            )

##Fitting with specific error tolerance
with open("converted_workspace.json") as file:
    spec = json.load(file)
ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
_, opt = get_backend()
init_pars = model.config.suggested_init()
par_bounds = model.config.suggested_bounds()
fixed_params = model.config.suggested_fixed()
_validate_fit_inputs(init_pars, par_bounds, fixed_params)
    # get fixed vals from the model
fixed_vals = [
    (index, init)
    for index, (init, is_fixed) in enumerate(zip(init_pars, fixed_params))
    if is_fixed
]

print(opt.minimize(
    twice_nll, observations, model, init_pars, par_bounds, fixed_vals, tolerance = 0.00001)
)





    



##NLL PLOTS
'''
with open("converted_workspace.json") as file:
    spec = json.load(file)
ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
poi_values = np.linspace(0, 3, 50).tolist()
nuisances = model.config.suggested_init()[1:]
NLL = []
for poi in poi_values:
    NLL.append(-1*(model.logpdf(pars=[poi] + nuisances, data=observations))-35.91)
##plt.plot(poi_values, NLL, label = "pyhf NLL")

file = uproot.open('higgsCombineTest.MultiDimFit.mH120.root')
tree = file['limit']
branches = tree.arrays()
plt.plot(branches['r'], branches['deltaNLL'], label  = "combine NLL")
file.close()

mu_values = branches['r']
nlls = []
for mu in mu_values:
    nll = 0
    rates_per_bin = model.expected_data([mu] + nuisances, include_auxdata=False)
    for i, bin_rate in enumerate(rates_per_bin):
        # build up negative log likelihood as sum over log Poisson likelihoods per bin
        nll +=  -scipy.stats.poisson.logpmf(observations[i], bin_rate)
    nlls.append(nll)

nlls = nlls - min(nlls)  # offset to set minimum to zero
plt.plot(mu_values, nlls, label = "pyhf_manual_NLL")
plt.xlabel("mu")
plt.ylabel("$\Delta$ NLL")
plt.legend()

plt.savefig('NLLplot')
'''
    




##YIELDS
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

##CLs plot
"""fig, ax = plt.subplots()
fig.set_size_inches(10.5, 7)
ax.set_title("Hypothesis Tests")

artists = brazil.plot_results(poi_values, results, ax=ax)
plt.savefig('shapes_brazil_band')
"""





