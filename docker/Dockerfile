ARG BASE_IMAGE=aperloff/cms-cvmfs-docker:latest
FROM ${BASE_IMAGE} as base

ENV CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org"

ARG CMSSW_VERSION="CMSSW_11_2_0"
ARG COMBINED_LIMIT_TAG="py3"
ARG COMBINE_HARVESTER_TAG="113x"

RUN . /cvmfs/cms.cern.ch/cmsset_default.sh && \
    cmsrel "${CMSSW_VERSION}" && \
    cd "${CMSSW_VERSION}"/src && \
    cmsenv && \
    git clone \
        --depth 1 \
        https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git \
        --branch "${COMBINED_LIMIT_TAG}" \
        --single-branch \
        HiggsAnalysis/CombinedLimit && \
    git clone \
        --depth 1 \
        https://github.com/nsmith-/CombineHarvester.git \
        --branch "${COMBINE_HARVESTER_TAG}" \
        --single-branch \
        CombineHarvester && \
    scram b --jobs $(($(nproc) - 1)) && \
    cmsenv && \
    python3 -m pip --no-cache-dir --verbose install --upgrade --user pip setuptools wheel && \
    python3 -m pip --no-cache-dir --verbose install --upgrade --user 'pyhf[xmlio,minuit,contrib]'

WORKDIR /home/cmsusr/"${CMSSW_VERSION}"/src
