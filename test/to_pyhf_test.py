
from array import array
from signal import Sigmasks
from numpy.core.fromnumeric import size
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

##make pull plots
'''
with open("./big_model/bottom-squarks.json") as file:
    spec = json.load(file)
ws = pyhf.Workspace(spec)
model = ws.model()
observations = ws.data(model)

pyhf_pulls = [1.95E-03,1.33E-03,5.49E-02,7.43E-02,6.26E-03,1.33E-03,5.02E-03,3.19E-03,-1.04E-02,-2.17E-04,-1.43E-02,-2.10E-02,5.83E-03,1.33E-03,1.33E-03,2.60E-03,7.66E-03,3.39E-03,1.22E-02,2.18E-03,1.71E-03,3.19E-03,3.64E-03,9.57E-04,-2.62E-03,2.01E-03,3.35E-03,1.15E-03,-1.37E-01,2.27E-01,1.85E-01,2.38E-02,5.60E-02,1.30E-01,-3.84E-02,-1.66E-01,-1.44E-02,-3.53E-02,1.42E-02,-4.68E-02,-6.84E-03,1.89E-02,7.90E-03,-7.21E-03,-4.97E-03,-1.29E-03,1.33E-03,2.62E-03,9.02E-03,3.75E-03]
combine_pulls = [-3.06E-06,2.21E-11,-3.00E-06,0.00E+00,-2.00E-07,2.21E-11,-1.00E-07,-5.00E-06,0.00E+00,-6.00E-05,-3.00E-05,-1.60E-05,1.82E-05,2.21E-11,2.21E-11,-6.00E-08,-4.00E-07,-1.00E-07,-8.00E-07,-5.00E-08,-3.00E-08,-1.00E-07,-1.00E-07,-6.00E-07,-1.00E-05,1.40E-11,4.80E-07,0.00E+00,7.00E-05,-8.00E-05,7.00E-05,-4.00E-05,3.63E-06,6.72E-05,4.57E-05,5.86E-05,1.21E-04,7.24E-06,8.89E-06,4.38E-05,-1.00E-06,-3.00E-05,3.00E-05,-2.70E-05,2.20E-05,0.00E+00,3.50E-11,0.00E+00,-2.00E-06,-4.00E-06]
pyhf_uncertainties = [1,1,1,1,1,1,1,1,1,1,1,1,5.14E-01,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,9.29E-01,1.38E+00,1,9.60E-01,8.42E-01,1.26E+00,1.00E+00,1.00E+00,1.47E+00,9.07E-01,1.00E+00,1.00E+00,7.80E-01,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00]
combine_uncertainties = [0.99,1,0.99,0.99,1,1,0.99,0.99,0.99,1,1,1,5.14E-01,1,1.01,1.01,1.01,1,1,1,1,1,1,1,1.01,1.01,1.01,1,9.29E-01,1.38E+00,1,9.60E-01,8.42E-01,1.26E+00,1.00E+00,1.00E+00,1.47E+00,9.07E-01,1.00E+00,1.00E+00,7.80E-01,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00,1.00E+00]

fig, ax = plt.subplots()
fig.set_size_inches(20, 5)

# set up axes labeling, ranges, etc...

ax.set_xlim(-0.5, len(pyhf_pulls) - 0.5)
ax.set_title("Pull Plot", fontsize=18)


# draw the +/- 2.0 horizontal lines
ax.hlines([-2, 2], -0.5, len(pyhf_pulls) - 0.5, colors="black", linestyles="dotted")
# draw the +/- 1.0 horizontal lines
ax.hlines([-1, 1], -0.5, len(pyhf_pulls) - 0.5, colors="black", linestyles="dashdot")
# draw the +/- 2.0 sigma band
ax.fill_between([-0.5, len(pyhf_pulls) - 0.5], [-2, -2], [2, 2], facecolor="yellow")
# drawe the +/- 1.0 sigma band
ax.fill_between([-0.5, len(pyhf_pulls) - 0.5], [-1, -1], [1, 1], facecolor="green")
# draw a horizontal line at pull=0.0
ax.hlines([0], -0.5, len(pyhf_pulls) - 0.5, colors="black", linestyles="dashed")
# finally draw the pulls
ax.scatter([i + 0.1 for i in range(len(pyhf_pulls))], pyhf_pulls, color="black", label = "pyhf params")
ax.errorbar(
    [i + 0.1 for i in range(len(pyhf_pulls))],
    pyhf_pulls,
    color="black",
    xerr=0,
    yerr=pyhf_uncertainties,
    marker=".",
    fmt="none",
)



# finally draw the pulls
ax.scatter(range(len(combine_pulls)), combine_pulls, color="red", label = "Combine params")
#uncertainties
ax.errorbar(
    range(len(combine_pulls)),
    combine_pulls,
    color="red",
    xerr=0,
    yerr=combine_uncertainties,
    marker=".",
    fmt="none",
)
ax.legend()

plt.savefig('pull_plot')
'''
