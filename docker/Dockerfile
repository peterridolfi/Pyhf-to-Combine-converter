ARG BASE_IMAGE=aperloff/cms-cvmfs-docker:latest
FROM ${BASE_IMAGE} as base

ENV CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org"

RUN . /cvmfs/cms.cern.ch/cmsset_default.sh && \
    cmsrel CMSSW_10_2_13 && \
    cd CMSSW_10_2_13/src && \
    cmsenv && \
    git clone \
        --depth 1 \
        https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git \
        --branch master \
        --single-branch \
        HiggsAnalysis/CombinedLimit && \
    git clone \
        --depth 1 \
        https://github.com/cms-analysis/CombineHarvester.git \
        --branch master \
        --single-branch \
        CombineHarvester && \
    scram b --jobs $(($(nproc) - 1)) && \
    cmsenv

WORKDIR /home/cmsusr/CMSSW_10_2_13/src
