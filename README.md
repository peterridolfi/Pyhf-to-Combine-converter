# Pyhf-to-Combine-converter
Tool to convert models from pyhf to CMS Combine and vise versa

This tool has been fully verified to produce the same model structure and expected yields, as well as to produce fits that are within 1% of the central value of each other. This has been done through the comparison of NLL plots, pull plots, and raw difference calculations.

In order to use the package, one must first create a docker container that is capable of running Combine. To do so, enter the following commands once you have docker installed.

```
docker pull pyhf/pyhf-combine-converter:cmssw-11.2.0-python3
docker run --rm -ti -P --device /dev/fuse --cap-add SYS_ADMIN --security-opt apparmor:unconfined -e CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org" pyhf/pyhf-combine-converter:cmssw-11.2.0-python3
```

If the docker container is running, but the Combine commands are not being mounted correctly, run the following commands within the container:

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv
```

After the container has been built, use the following commands to run the tool

For pyhf->combine:

```
python3 pyhf_convert_to_datacard.py $JSON_FILE_NAME --shape-file $SHAPES_FILE_NAME --out-datacard $DATACARD_FILE_NAME
```

For combine->pyhf:

```
python3 pyhf_converted_from_datacard.py $DATACARD_FILE_NAME --out-file $JSON_FILE_NAME


Any questions or issues should be referred to the docs/ folder, in which the translation is put into more detail.