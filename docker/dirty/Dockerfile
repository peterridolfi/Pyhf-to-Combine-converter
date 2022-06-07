ARG BASE_IMAGE=pyhf/pyhf-combine-converter:commit-build
FROM ${BASE_IMAGE} as base

ARG CMSSW_VERSION="CMSSW_11_2_0"

# RUN echo 'export PATH=/home/cmsusr/.local/bin:${PATH}' >> /home/cmsusr/.bashrc && \
#     echo 'printf "\n# Setup the environment with:\n# $ cmsenv\n"' >> /home/cmsusr/.bashrc
RUN echo 'printf "\n# Setup the environment with:\n# $ cmsenv\n"' >> /home/cmsusr/.bashrc

ENV CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org"

WORKDIR /home/cmsusr/"${CMSSW_VERSION}"/src
