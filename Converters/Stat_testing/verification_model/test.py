import uproot

file = uproot.open("shapes-test.root")
bkg1 = file["bkg1"]
print(bkg1.values())

file.close()
