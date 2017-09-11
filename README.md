# 2D/3D Material Point Method (mpm)
> Cambridge Berkeley - Geomechanics

> Based on [Aspect](https://github.com/geodynamics/aspect)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/cb-geo/mpm/develop/license.md)
[![Developer docs](https://img.shields.io/badge/developer-docs-blue.svg)](http://cb-geo.github.io/mpm)
[![User docs](https://img.shields.io/badge/user-docs-blue.svg)](https://mpm.cb-geo.com/)
[![](https://img.shields.io/github/issues-raw/cb-geo/mpm.svg)](https://github.com/cb-geo/mpm/issues)
[![Project management](https://img.shields.io/badge/projects-view-ff69b4.svg)](https://github.com/cb-geo/mpm/projects/)

## Install dependencies

* Docker image for CB-Geo mpm code [https://hub.docker.com/r/cbgeo/mpm](https://hub.docker.com/r/cbgeo/mpm)

* Instructions for running mpm docker container: [https://github.com/cb-geo/docker-mpm/blob/master/README.md](https://github.com/cb-geo/mpm-container/blob/master/README.md).

### Prerequisite packages
> The following prerequisite packages can be found in the docker image:

* [Dealii](http://dealii.org/)

## Compile and Run
> See https://mpm-doc.cb-geo.com/ for more detailed instructions. 

0. Run `mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release /path/to/CMakeLists.txt`.

1. Run `make clean && make -jN` (where N is the number of cores).
