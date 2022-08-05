
from array import array
from signal import Sigmasks
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
from operator import indexOf


##Fitting with specific error tolerance
'''
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
'''



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
    


'''
with open("./multi_vs_single_bin/single_bin_1_test.json") as file:
    spec1a = json.load(file)
ws1a = pyhf.Workspace(spec1a)
model1a = ws1a.model()
observations1a = ws1a.data(model1a)
with open("./multi_vs_single_bin/single_bin_2_test.json") as file:
    spec1b = json.load(file)
ws1b = pyhf.Workspace(spec1b)
model1b = ws1b.model()
observations1b = ws1b.data(model1b)



file = uproot.open('./multi_vs_single_bin/higgsCombineTest.MultiDimFit.mH120_1.root')
tree = file['limit']
branches = tree.arrays()
NLL = []
for nll in branches['deltaNLL'][1:]:
    NLL.append(nll - min(branches['deltaNLL']))

file.close()



mu_values = branches['r'][1:]
nlls = []
for mu in mu_values:
    nll = 0
    rates_per_bin = model1a.expected_data([mu, 0], include_auxdata=False)
    for i, bin_rate in enumerate(rates_per_bin):
        # build up negative log likelihood as sum over log Poisson likelihoods per bin
        nll +=  -scipy.stats.poisson.logpmf(observations1a[i], bin_rate)
    nlls.append(nll)

nlls = nlls - min(nlls)  # offset to set minimum to zero
plt.plot(mu_values, nlls, label = "singlebin1 pyhf scan")
plt.plot(mu_values, NLL, label = "singlebin1 combine scan")




file = uproot.open('./multi_vs_single_bin/higgsCombineTest.MultiDimFit.mH120_2.root')
tree = file['limit']
branches = tree.arrays()
NLL = []
for nll in branches['deltaNLL'][1:]:
    NLL.append(nll - min(branches['deltaNLL']))

file.close()



mu_values = branches['r'][1:]
nlls = []
for mu in mu_values:
    nll = 0
    rates_per_bin = model1b.expected_data([mu, 0], include_auxdata=False)
    for i, bin_rate in enumerate(rates_per_bin):
        # build up negative log likelihood as sum over log Poisson likelihoods per bin
        nll +=  -scipy.stats.poisson.logpmf(observations1b[i], bin_rate)
    nlls.append(nll)

nlls = nlls - min(nlls)  # offset to set minimum to zero
plt.plot(mu_values, nlls, label = "singlebin2 pyhf scan")
plt.plot(mu_values, NLL, label = "singlebin2 combine scan")




file = uproot.open('./multi_vs_single_bin/higgsCombineTest.MultiDimFit.mH120_multi.root')
tree = file['limit']
branches = tree.arrays()
NLL = []
for nll in branches['deltaNLL'][1:]:
    NLL.append(nll - min(branches['deltaNLL']))

file.close()



mu_values = branches['r'][1:]
nlls = []
for mu in mu_values:
    nll = 0
    rates_per_bin = model2.expected_data([mu, 0], include_auxdata=False)
    for i, bin_rate in enumerate(rates_per_bin):
        # build up negative log likelihood as sum over log Poisson likelihoods per bin
        nll +=  -scipy.stats.poisson.logpmf(observations2[i], bin_rate)
    nlls.append(nll)

nlls = nlls - min(nlls)  # offset to set minimum to zero
plt.plot(mu_values, nlls, label = "two_bin pyhf scan")
plt.plot(mu_values, NLL, label = "two bin combine scan")
plt.xlabel("mu")
plt.ylabel("$\Delta$ NLL")
plt.legend()
plt.savefig('NLLplot')

print(pyhf.infer.mle.fit(pdf = model2, data = observations2))
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


'''
##NLL plots for fitting, minimization
with open("./verification_model/converted_workspace_floats.json") as file:
    spec2 = json.load(file)
ws2 = pyhf.Workspace(spec2)
model2 = ws2.model()
observations2 = ws2.data(model2)



file = uproot.open('./verification_model/higgsCombineTest.MultiDimFit.mH120_fitscan.root')
tree = file['limit']
branches = tree.arrays()
mu_values = branches['r'][1:]
NLL = []
for nll in branches['deltaNLL'][1:]:
    NLL.append(nll - min(branches['deltaNLL']))
plt.plot(mu_values, NLL, label = "combine scan")


sigma = branches['sigma'][0]
shape = branches['shape'][0]
norm1 = branches['norm1'][0]
shape1 = branches['shape1'][0]
norm2 = branches['norm2'][0]
norm3 = branches['norm3'][0]
shape3 = branches['shape3'][0]

staterr0 = branches['prop_binb1_bin0'][0]
staterr1 = branches['prop_binb1_bin1'][0]
staterr2 = branches['prop_binb1_bin2'][0] 
staterr3 = branches['prop_binb1_bin3'][0]
staterr4 = branches['prop_binb1_bin4'][0]
staterr5 = branches['prop_binb1_bin5'][0]
staterr6 = branches['prop_binb1_bin6'][0]
staterr7 = branches['prop_binb1_bin7'][0]
staterr8 = branches['prop_binb1_bin8'][0]
staterr9 = branches['prop_binb1_bin9'][0]
'''



'''
mu_values = branches['r'][1:]
nlls = []
for mu in mu_values:
    nll = 0
    pyhf.set_backend("numpy", "minuit")
    params, n = pyhf.infer.mle.fit(data = observations2, pdf = model2, init_pars = [mu, 0], fixed_params = [True, False], return_fitted_val = True)
    rates_per_bin = model2.expected_data(params, include_auxdata=False)
    for i, bin_rate in enumerate(rates_per_bin):
        # build up negative log likelihood as sum over log Poisson likelihoods per bin
        nll +=  -scipy.stats.poisson.logpmf(observations2[i], bin_rate)
    nlls.append(nll)
plt.plot(mu_values, nlls, label = "pyhf scan")


'''

'''
NLL = []
for poi in mu_values:
    params, nll = pyhf.infer.mle.fit(data = observations2, pdf = model2, init_pars = [poi, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], fixed_params = [True, False, False, False, False, False, False, False, False, False, False], return_fitted_val = True)
    NLL.append(nll)


  # offset to set minimum to zero
NLL = NLL - min(NLL)
plt.plot(mu_values, [i/2 for i in NLL], label = "pyhf scan")

print(pyhf.infer.mle.fit(data = observations2, pdf = model2))


plt.xlabel("mu")
plt.ylabel("$\Delta$ NLL")
plt.legend()
plt.savefig('NLLplot')

file.close()'''

with open("./big_model/bottom-squarks.json") as file:
    spec = json.load(file)
ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)
for param in pyhf.infer.mle.fit(data = observations, pdf = model):
    print(param)
for par in model.config.par_order:
    print(par)



