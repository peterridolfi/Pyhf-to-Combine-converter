# Convenient access to the version number
from pyhf_combine_converter.version import version as __version__

from pyhf_combine_converter.pyhf_convert_to_datacard import pyhf_convert_to_datacard
from pyhf_combine_converter.pyhf_converted_from_datacard import pyhf_converted_from_datacard

__all__ = ["__version__", "pyhf_convert_to_datacard", "pyhf_converted_from_datacard"]


def __dir__():
    return __all__
