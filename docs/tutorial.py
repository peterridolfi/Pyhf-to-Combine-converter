from pyhf_combine_converter import pyhf_converted_from_datacard
from pyhf_combine_converter import pyhf_convert_to_datacard
import HiggsAnalysis.CombinedLimit.DatacardParser as DP
from optparse import OptionParser

##convert to pyhf
parser = OptionParser()
DP.addDatacardParserOptions(parser)
options, args = parser.parse_args()
pyhf_converted_from_datacard(input_datacard = "converted_datacard.txt", outfile = "pyhf_workspace.json", options = options)
##convert back to datacard
pyhf_convert_to_datacard(workspace = "pyhf_workspace.json", outdatacard = "combine_datacard.txt", shapefile = "shapes.root", options = None)

