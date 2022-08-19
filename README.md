# pyhf to Combine converter

Tool to convert models from pyhf to CMS Combine and vise versa

[![DOI](https://zenodo.org/badge/492820903.svg)](https://zenodo.org/badge/latestdoi/492820903)

[![PyPI version](https://badge.fury.io/py/pyhf-combine-converter.svg)](https://badge.fury.io/py/pyhf-combine-converter)

This tool has been fully verified to produce the same model structure and expected yields, as well as to produce fits that are within 1% of the central value of each other. This has been done through the comparison of NLL plots, pull plots, and raw difference calculations.

## Environment for Combine

In order to use the package, one must first create a docker container that is capable of running Combine. To do so, enter the following commands once you have Docker installed.

```
docker pull pyhf/pyhf-combine-converter:cmssw-11.2.0-python3
docker run --rm -ti -P --device /dev/fuse --cap-add SYS_ADMIN --security-opt apparmor:unconfined -e CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org" pyhf/pyhf-combine-converter:cmssw-11.2.0-python3
```

If the docker container is running, but the Combine commands are not being mounted correctly, run the following commands within the container:

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
```

## Install

`pyhf-combine-converter` can be installed from PyPI. Inside of the Python environment containing Combine run

```
python3 -m pip install pyhf-combine-converter
```

## Use

`pyhf-combine-converter` provides a CLI API for bidirectional conversion between pyhf and Combine.

### Convert from pyhf to Combine

```
pyhf-to-combine $JSON_FILE_NAME --shape-file $SHAPES_FILE_NAME --out-datacard $DATACARD_FILE_NAME
```

### Convert from Combine to pyhf

```
combine-to-pyhf $DATACARD_FILE_NAME --out-file $JSON_FILE_NAME
```

Any questions or issues should be referred to the docs/ folder, in which the translation is put into more detail.

## Acknowledgements

This work was done as part of [Peter Ridolfi's IRIS-HEP Fellow project](https://iris-hep.org/fellows/peterridolfi.html) which was supported by NSF
cooperative agreement [OAC-1836650][NSF-award-1836650].
[![NSF Award Number](https://img.shields.io/badge/NSF-1836650-blue.svg)][NSF-award-1836650]

[NSF-award-1836650]: https://nsf.gov/awardsearch/showAward?AWD_ID=1836650
