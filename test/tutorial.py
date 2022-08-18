import pyhf_combine_converter as converter

##convert to pyhf
converter.pyhf_converted_from_datacard.main("converted_datacard.txt", outfile = "pyhf_workspace.json")
##convert back to datacard
converter.pyhf_convert_to_datacard.main("pyhf_workspace.json", outdatacard = "combine_datacard.txt", shapefile = "shapes.root")

