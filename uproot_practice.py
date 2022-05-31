import uproot
import matplotlib.pyplot as plt
import numpy as np
import awkward as ak
import vector

file = uproot.open('uproot-tutorial-file.root')

tree = file['Events']
branches  =tree.arrays()
print(branches)
two_muons_mask = branches['nMuon']==2
muon_p4 = vector.zip({'pt': branches['Muon_pt'], 'eta': branches['Muon_eta'], 'phi': branches['Muon_phi'], 'mass': branches['Muon_mass']})
two_muons_p4 = muon_p4[two_muons_mask]

first_muon_p4 = two_muons_p4[:, 0]
second_muon_p4 = two_muons_p4[:, 1]

sum_p4 = first_muon_p4 + second_muon_p4

two_muons_charges = branches['Muon_charge'][two_muons_mask]
opposite_sign_muons_mask = two_muons_charges[:, 0] != two_muons_charges[:, 1]
dimuon_p4 = sum_p4[opposite_sign_muons_mask]

plt.hist(dimuon_p4.mass, bins=40, range=(70, 110))
plt.xlabel('Dimuon invariant mass [GeV]')
plt.ylabel('Number of dimuon events / 1 MeV')
plt.show()


