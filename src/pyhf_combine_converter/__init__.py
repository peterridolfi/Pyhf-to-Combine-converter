# FIXME: pyhf import is needed only to avoid UserWarnings
# c.f. https://github.com/peterridolfi/Pyhf-to-Combine-converter/issues/30
import pyhf

from pyhf_combine_converter import combine_to_pyhf, pyhf_to_combine
from pyhf_combine_converter.combine_to_pyhf import pyhf_converted_from_datacard
from pyhf_combine_converter.pyhf_to_combine import pyhf_convert_to_datacard

# Convenient access to the version number
from pyhf_combine_converter.version import version as __version__

__all__ = [
    "__version__",
    "combine_to_pyhf",
    "pyhf_convert_to_datacard",
    "pyhf_converted_from_datacard",
    "pyhf_to_combine",
]


def __dir__():
    return __all__
