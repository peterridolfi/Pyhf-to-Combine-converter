
import pyhf
import pyhf.readxml
import json

spec = pyhf.readxml.parse(
    "C:\Users\peter\pyhf-tutorial\book\data\multichannel_histfactory\config\example.xml", "C:\Users\peter\pyhf-tutorial\book\data\multichannel_histfactory"
)

workspace = pyhf.Workspace(spec)
print(json.dumps(workspace, indent=2))