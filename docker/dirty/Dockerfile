ARG BASE_IMAGE=pyhf/pyhf-combine-converter:commit-build
FROM ${BASE_IMAGE} as base

RUN echo 'printf "\n# Setup the environment with:\n# $ cmsenv\n"' >> /home/cmsusr/.bashrc

ENV CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org"

WORKDIR /home/cmsusr/CMSSW_10_2_13/src
