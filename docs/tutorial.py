from pyhf_combine_converter import pyhf_converted_from_datacard
from pyhf_combine_converter import pyhf_convert_to_datacard

# convert to pyhf
pyhf_converted_from_datacard.pyhf_converted_from_datacard(
    "converted_datacard.txt", outfile="pyhf_workspace.json"
)

# # convert back to datacard
# pyhf_convert_to_datacard.pyhf_convert_to_datacard(
#     "pyhf_workspace.json", outdatacard="combine_datacard.txt", shapefile="shapes.root"
# )
