TAG := v3
MANAGER := podman

.PHONY: all clean

all:
	echo "VERSION = \"$(shell git log -1 --format="%h (%ai)")\"" > annovep/_version.py
	${MANAGER} build -t annovep:${TAG} .

clean:
# Needs podman v3.3.1: https://github.com/containers/podman/issues/10832
	${MANAGER} image prune --force
	buildah rm --all
