default: image

all: image

image:
	docker pull aperloff/cms-cvmfs-docker:latest
	docker build . \
		-f docker/Dockerfile \
		--build-arg BASE_IMAGE=aperloff/cms-cvmfs-docker:latest \
		--build-arg CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org" \
		--tag pyhf/pyhf-combine-converter:debug-local

run:
	docker run \
		--rm \
		-it \
		-P \
		--device /dev/fuse \
		--cap-add SYS_ADMIN \
		--security-opt apparmor:unconfined \
		-e CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org" \
		aperloff/cms-cvmfs-docker:latest

image_dirty:
	docker build . \
		-f docker/dirty/Dockerfile \
		--build-arg BASE_IMAGE=pyhf/pyhf-combine-converter:commit-build \
		--build-arg CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org" \
		--tag pyhf/pyhf-combine-converter:latest \
		--tag pyhf/pyhf-combine-converter:cmssw-11.2.0-python3

run_dirty:
	docker run \
		--rm \
		-it \
		-P \
		--device /dev/fuse \
		--cap-add SYS_ADMIN \
		--security-opt apparmor:unconfined \
		-e CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org" \
		pyhf/pyhf-combine-converter:latest
