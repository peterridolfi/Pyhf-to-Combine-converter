ARG BASE_IMAGE=pyhf/pyhf-combine-converter:commit-build
FROM ${BASE_IMAGE} as base

ARG CMSSW_VERSION="CMSSW_11_2_0"
ENV CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org"

WORKDIR /home/cmsusr/"${CMSSW_VERSION}"/src
