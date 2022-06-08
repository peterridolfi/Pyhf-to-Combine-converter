# Pyhf-to-Combine-converter
Tool to convert models from phyf to CMS Combine and vise versa

## Docker image

A base Docker image that contains a Python 3 environment with CMS Combine and pyhf installed is available on [Docker Hub][Docker Hub]

```console
$ docker pull pyhf/pyhf-combine-converter:cmssw-11.2.0-python3
```

### Running the image

As noted in the [using the cms-cvmfs-docker image tutorial][cvmfs-image-tutorial], the images can be run locally using the following run commands (all options shown are required).

#### Linux

On Linux systems the [`--security-opt` option][security-opt-option] is needed

```console
$ docker run \
    --rm \
    -ti \
    -P \
    --device /dev/fuse \
    --cap-add SYS_ADMIN \
    --security-opt apparmor:unconfined \
    -e CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org" \
    pyhf/pyhf-combine-converter:cmssw-11.2.0-python3
```

#### macOS

```console
$ docker run \
    --rm \
    -ti \
    -P \
    --device /dev/fuse \
    --cap-add SYS_ADMIN \
    -e CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org" \
    pyhf/pyhf-combine-converter:cmssw-11.2.0-python3
```

[Docker Hub]: https://hub.docker.com/r/pyhf/pyhf-combine-converter/tags
[security-opt-option]: https://docs.docker.com/engine/reference/commandline/run/#optional-security-options---security-opt
[cvmfs-image-tutorial]: https://awesome-workshop.github.io/docker-singularity-hats/07-cms-cvmfs-docker/index.html
