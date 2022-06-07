default: image

all: image

image:
	docker pull aperloff/cms-cvmfs-docker:latest
	docker build . \
	-f docker/Dockerfile \
	--build-arg BASE_IMAGE=aperloff/cms-cvmfs-docker:latest \
	--tag pyhf/pyhf-combine-converter:debug-local

run:
	docker run \
		-it \
		-P \
		--device /dev/fuse \
		--cap-add SYS_ADMIN \
		--security-opt apparmor:unconfined \
		-e CVMFS_MOUNTS="cms.cern.ch oasis.opensciencegrid.org" \
		pyhf/pyhf-combine-converter:debug-local
