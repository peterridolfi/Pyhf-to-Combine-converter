#convert datacard to pyhf
combine-to-pyhf converted_datacard.txt --out-file pyhf_workspace.json
#convert pyhf back to datacard
pyhf-to-combine pyhf_workspace.json --out-datacard combine_datacard.txt --shape-file shapes.root
