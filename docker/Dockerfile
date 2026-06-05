#Dockerfile
#Engineer: T Looby
#Date: 03/20/2024
#Description:  Dockerfile for building HEAT docker image
#
# This dockerfile is broken into two stages.  The first stage (HEATbuilder)
# builds OpenFOAM, swak4Foam, M3DC1, MAFOT, and heatFoam.  Order is cache-friendly:
# HEAT is cloned only when building heatFoam, so changing HEAT_REF only rebuilds
# that layer (clone + wmake heatFoam), not the long OpenFOAM/swak4Foam build.
# HEATbuilder can be built once and pushed to Docker Hub; then the final image
# can use "FROM plasmapotential/heatbuilder:<tag>" and skip stage 1 entirely.
#
# The second stage is the final HEAT container.  It only contains
# the packages required to run HEAT, and none of the files used
# to build everything from source.  
#
# Build from the repository root (parent of docker/), so COPY paths are correct:
#   docker build -f docker/Dockerfile -t plasmapotential/heat .
# With a specific HEAT branch/tag for the image:
#   docker build -f docker/Dockerfile --build-arg HEAT_REF=v4.2.6 -t plasmapotential/heat:v4.2.6 .
#
# Cross platform build (e.g. apple silicon):
#   docker buildx build --platform linux/amd64,linux/arm64 -f docker/Dockerfile -t plasmapotential/heat . --push
#
# Build first stage only (HEATbuilder):
#   docker build -f docker/Dockerfile --target HEATbuilder -t plasmapotential/heatbuilder .
#
# after creating an image, tag it using a command like this
#   docker tag plasmapotential/heat:latest plasmapotential/heat:<tag>
#
# upload to dockerhub via a command like this
#   docker push plasmapotential/heat:<tag>
#
# launch bash terminal in container interatively
#   docker-compose run HEAT /bin/bash
#
#======================================================================
#1st build stage.  build all necessary packages from source
#======================================================================
# start from base
FROM ubuntu:24.04 AS heatbuilder
SHELL ["/bin/bash", "-c"] 

LABEL maintainer="tlooby@cfs.energy"

# environment variables
ENV SHELL /bin/bash
ENV runMode docker
ENV DEBIAN_FRONTEND noninteractive
ENV PYTHONPATH /usr/lib/freecad-python3/lib

#update layer
RUN apt-get -yqq update

#python layer
RUN apt-get -yqq install python3 && \
    apt-get -yqq install python3-pip && \
    apt-get -yqq install python3-pkg-resources
#RUN apt-get -yqq install mdsplus-python

# Create venv for PEP 668 (ubuntu 24.04+)
#RUN apt-get -yqq install python3-venv && \
#    python3 -m venv /opt/venv && \
#    /opt/venv/bin/pip install --upgrade pip
## Make sure venv is always active
#ENV VIRTUAL_ENV=/opt/venv
#ENV PATH="$VIRTUAL_ENV/bin:$PATH"

#utility tools/libs layers
#RUN apt-get -yqq libnglib-6.2 
RUN apt-get -yqq install coreutils libcurl4 libunwind-dev nano git git-lfs gfortran iputils-ping iproute2 htop 

#3D plasmas layer
RUN apt-get -yqq install libopenmpi-dev libhdf5-openmpi-dev libnetcdf-dev libnetcdff-dev libblas-dev liblapack-dev

#linking layer
RUN ln -s /usr/lib/x86_64-linux-gnu/libunwind.so /usr/lib/x86_64-linux-gnu/libunwind.so.1

#clone EFIT
#RUN git clone https://github.com/ORNL-Fusion/EFIT.git /root/source/EFIT

WORKDIR /root

#clone openfoam, third party libs
#OpenFOAM version defined by user
ARG OFversion=2306
RUN apt-get -yqq update && \
    apt-get install -yqq python-is-python3 flex bison wget mercurial libscotch-dev libptscotch-dev && \
    wget -P /root/source/openfoam/ https://dl.openfoam.com/source/v${OFversion}/ThirdParty-v${OFversion}.tgz && \
    wget -P /root/source/openfoam/ https://dl.openfoam.com/source/v${OFversion}/OpenFOAM-v${OFversion}.tgz && \
    tar -xzf /root/source/openfoam/ThirdParty-v${OFversion}.tgz -C /root/source/openfoam/ && \
    tar -xzf /root/source/openfoam/OpenFOAM-v${OFversion}.tgz -C /root/source/openfoam/ && \
    mkdir -p /root/source/openfoam/swak4Foam

#build openfoam (without heatFoam first so this layer does not depend on HEAT_REF and stays cached)
ENV WM_THIRD_PARTY_DIR /root/source/openfoam/ThirdParty-v${OFversion}
WORKDIR /root/source/openfoam/OpenFOAM-v${OFversion}
RUN mkdir -p applications/solvers/custom && \
    source /root/source/openfoam/OpenFOAM-v${OFversion}/etc/bashrc && \
    ./Allwmake -j -s -q -l

# clone HEAT and build only heatFoam (this is the only layer that depends on HEAT_REF; OpenFOAM already built above)
ARG HEAT_REF=v4.3
RUN git clone -b ${HEAT_REF} --single-branch https://github.com/plasmapotential/HEAT.git /root/source/HEAT && \
    cp -au /root/source/HEAT/source/heatFoam /root/source/openfoam/OpenFOAM-v${OFversion}/applications/solvers/custom/heatFoam && \
    source /root/source/openfoam/OpenFOAM-v${OFversion}/etc/bashrc && \
    cd /root/source/openfoam/OpenFOAM-v${OFversion}/applications/solvers/custom/heatFoam && \
    wmake

#clone and build swak4foam
WORKDIR /root/source/openfoam/
RUN source /root/source/openfoam/OpenFOAM-v${OFversion}/etc/bashrc && \
    hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam && \
    cd /root/source/openfoam/swak4Foam && \
    hg update develop && \
    ./AllwmakeAll

#clone and build M3DC1 fusion-io libs (buildM3DC1 uses make -k; doc/latex errors are harmless)
RUN git clone https://github.com/nferraro/fusion-io.git /root/source/M3DC1/fusion-io
COPY docker/buildM3DC1 /root/source/M3DC1/fusion-io/buildM3DC1
COPY docker/make.inc.ubuntu_x86 /root/source/M3DC1/fusion-io/install/make.inc.ubuntu_x86
# Upstream util/makefile forces Intel mpiicpc / -qopenmp; Ubuntu image uses GCC + Open MPI only.
# we create a custom makefile that uses GCC + Open MPI.
COPY docker/fusion-io-util-makefile /root/source/M3DC1/fusion-io/util/makefile
# write_neo_input.cpp calls MPI but omits <mpi.h>; mpiicpc often still compiles, Open MPI mpic++ does not.
#this needs to be fixed in fusion-io upstream, but for now we can add the include here.
RUN grep -q '#include <mpi.h>' /root/source/M3DC1/fusion-io/util/write_neo_input.cpp || \
    sed -i '/#include <netcdf.h>/a#include <mpi.h>' /root/source/M3DC1/fusion-io/util/write_neo_input.cpp
RUN mkdir -p /root/source/M3DC1/build && /root/source/M3DC1/fusion-io/buildM3DC1

#install CUDA toolkit (nvcc + cudart-dev) so MAFOT can compile GPU support.
#cudart is linked STATICALLY into the MAFOT binaries (see make.inc.HEAT / makefile),
#so the runtime stage needs no CUDA packages -- only the host driver, injected by
#`docker run --gpus`, and only when -g is actually used.
RUN apt-get update && apt-get install -y --no-install-recommends wget ca-certificates && \
    wget -q https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2404/x86_64/cuda-keyring_1.1-1_all.deb -O /tmp/cuda-keyring.deb && \
    dpkg -i /tmp/cuda-keyring.deb && rm /tmp/cuda-keyring.deb && \
    apt-get update && \
    apt-get install -y --no-install-recommends cuda-nvcc-12-6 cuda-cudart-dev-12-6 && \
    ln -sfn /usr/local/cuda-12.6 /usr/local/cuda && \
    rm -rf /var/lib/apt/lists/*

#clone and build MAFOT
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/hdf5/openmpi
ARG MAFOT_BRANCH_CACHE_BUST=1
RUN git clone -b mafot_gpu --single-branch https://github.com/ORNL-Fusion/MAFOT.git /root/source/MAFOT
COPY docker/buildMAFOT /root/source/MAFOT/buildMAFOT
RUN mkdir -p /root/source/MAFOT/build/bin && \
    mkdir -p /root/source/MAFOT/build/lib && \
    /root/source/MAFOT/buildMAFOT

#clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

#======================================================================
#2nd build stage.  build container without unncesary build artifacts
#to keep final image smaller
#======================================================================
#if user doesnt want to build first stage, they can pull HEATbuilder from
#dockerhub using the following command (need to comment out first stage too):
#FROM plasmapotential/HEATbuilder:<tag> as HEATbuilder

FROM ubuntu:24.04
LABEL maintainer="tlooby@cfs.energy"
ARG DEBIAN_FRONTEND=noninteractive
ENV runMode docker

#update layer
RUN apt-get -yqq update && apt-get -yqq upgrade

#python layer
RUN apt-get -yqq install python3 && \
    apt-get -yqq install python3-pip && \
    apt-get -yqq install python3-pkg-resources && \
    apt-get -yqq install python3-paraview

#RUN apt-get -yqq install mdsplus-python

# freecad layer
#RUN apt-get -yqq install libfreecad-python3-0.20 #22.04
#RUN apt-get -yqq install libfreecad-python3-0.19 #22.10
#RUN apt-get -yqq install software-properties-common && \
#    add-apt-repository ppa:freecad-maintainers/freecad-daily && \
#    apt-get update && \
#    apt-get install -yqq freecad-daily-python3
RUN apt-get -yqq install software-properties-common && \
    add-apt-repository ppa:freecad-maintainers/freecad-stable && \
    apt-get update && \
    apt-get install -yqq freecad-python3

# elmer layer / OF layer
RUN apt-get -yqq update && \
    apt-add-repository ppa:elmer-csc-ubuntu/elmer-csc-ppa && \
    apt-get update && \
    apt-get install -yqq elmerfem-csc libglu1 libscotch-dev libptscotch-dev python-is-python3

#RUN apt-get -yqq libnglib-6.2 
RUN apt-get -yqq update && apt-get -yqq install coreutils libcurl4 libunwind-dev nano git gfortran iputils-ping iproute2 htop git-lfs

#3D plasmas layer
RUN apt-get -yqq install libopenmpi-dev libhdf5-openmpi-dev libnetcdf-dev libnetcdff-dev

#linking layer
RUN ln -s /usr/lib/x86_64-linux-gnu/libunwind.so /usr/lib/x86_64-linux-gnu/libunwind.so.1

#clone EFIT (note branch can be changed)
RUN git clone -b master --single-branch https://github.com/ORNL-Fusion/EFIT.git /root/source/EFIT

#clone HEAT (branch or tag; set via: docker build --build-arg HEAT_REF=v4.2.6 ...)
ARG HEAT_REF=v4.3
# Docker caches the git clone; bump HEAT_CACHE_BUST (e.g. --build-arg HEAT_CACHE_BUST=$(date +%s))
# to force a fresh clone after pushing HEAT changes (e.g. requirements.txt) to the HEAT_REF branch.
ARG HEAT_CACHE_BUST=1
RUN git clone -b ${HEAT_REF} --single-branch https://github.com/plasmapotential/HEAT.git /root/source/HEAT
WORKDIR /root/source/HEAT
RUN     git lfs install && git lfs pull
WORKDIR /root

# fetch HEAT specific python dependencies
# PIP_BREAK_SYSTEM_PACKAGES: Ubuntu 24.04 Python 3.12 marks the env PEP 668
# "externally managed"; HEAT needs its pip deps in the same site-packages as the
# apt-installed freecad/paraview, so install system-wide. --ignore-installed on the
# pip self-upgrade avoids "Cannot uninstall pip ... RECORD file not found".
# scipy: requirements.txt leaves it unpinned so pip defers to the apt python3-scipy
# (built for numpy 1.x); but the pip stack pulls numpy 2.x -> ABI break
# ("numpy.dtype size changed"). Force the pip scipy build so everything is numpy 2.x.
RUN export PIP_BREAK_SYSTEM_PACKAGES=1 && \
    python3 -m pip install --upgrade --ignore-installed pip && \
    python3 -m pip install --ignore-installed blinker && \
    python3 -m pip install -r /root/source/HEAT/docker/requirements.txt && \
    python3 -m pip install --ignore-installed scipy && \
    python3 -m pip install -U --ignore-installed matplotlib

# PyFoam 2023.7 vendors six 1.10.0, whose six.moves importer relies on the find_module
# protocol that Python 3.12 (Ubuntu 24) removed -> "No module named PyFoam.ThirdParty.six.moves".
# Overlay the pip-installed six (>=1.13, has find_spec) onto PyFoam's vendored copy.
RUN python3 -c "import six,os,shutil; d=os.path.dirname(six.__file__); shutil.copy(six.__file__, os.path.join(d,'PyFoam','ThirdParty','six.py'))"


#make build dirs
RUN mkdir -p /root/builds/M3DC1 && \
    mkdir -p /root/builds/MAFOT && \
    mkdir -p /root/builds/openfoam && \
    mkdir -p /root/builds/swak4foam

#copy M3DC1 fusion-io bins / libs
COPY --from=heatbuilder /root/source/M3DC1/build/bin /root/builds/M3DC1/bin
COPY --from=heatbuilder /root/source/M3DC1/build/lib /root/builds/M3DC1/lib
COPY --from=heatbuilder /root/source/M3DC1/build/include /root/builds/M3DC1/include
COPY --from=heatbuilder /root/source/M3DC1/build/modules /root/builds/M3DC1/modules

#copy MAFOT bins / libs
COPY --from=heatbuilder /root/source/MAFOT/build/bin /root/builds/MAFOT/bin
COPY --from=heatbuilder /root/source/MAFOT/build/lib /root/builds/MAFOT/lib

#copy openFOAM bin / libs
ENV OFversion 2306 
COPY --from=heatbuilder /root/source/openfoam/OpenFOAM-v${OFversion}/platforms /root/builds/openfoam/platforms
COPY --from=heatbuilder /root/source/openfoam/OpenFOAM-v${OFversion}/build /root/builds/openfoam/build
COPY --from=heatbuilder /root/source/openfoam/OpenFOAM-v${OFversion}/bin /root/builds/openfoam/bin
COPY --from=heatbuilder /root/source/openfoam/OpenFOAM-v${OFversion}/etc /root/builds/openfoam/etc

#copy swak4foam bins / libs
COPY --from=heatbuilder /root/OpenFOAM/user-v2306/platforms/linux64GccDPInt32Opt/bin /root/builds/swak4Foam/bin
COPY --from=heatbuilder /root/OpenFOAM/user-v2306/platforms/linux64GccDPInt32Opt/lib /root/builds/swak4Foam/lib

#paths for libs
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/qt5
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/netgen
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/builds/MAFOT/lib
#ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/paraview/lib #comment to prevent qt5 conflicts
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/builds/openfoam/lib
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/builds/openfoam/platforms/linux64GccDPInt32Opt/lib
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/builds/openfoam/platforms/linux64GccDPInt32Opt/lib/sys-openmpi
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/builds/openfoam/platforms/linux64GccDPInt32Opt/lib/dummy
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/builds/swak4Foam/lib
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/builds/M3DC1/lib
ENV PYTHONPATH $PYTHONPATH:/root/builds/M3DC1/lib

#set up paths for bins
ENV PATH $PATH:/root/builds/MAFOT/bin
#ENV PATH $PATH:/root/builds/paraview/bin
ENV PATH $PATH:/root/builds/openfoam/bin
ENV PATH $PATH:/root/builds/openfoam/platforms/linux64GccDPInt32Opt/bin
ENV PATH $PATH:/root/builds/swak4foam/bin
ENV PATH $PATH:/root/builds/M3DC1/bin
ENV PYTHONPATH $PYTHONPATH:/usr/lib/freecad-python3/lib

# expose port
EXPOSE 8050

WORKDIR /root

#clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# start app
# start app — ENTRYPOINT is fixed, CMD provides defaults (overridable at runtime)
# GUI mode: docker run plasmapotential/heat (uses defaults below)
# TUI mode: docker run plasmapotential/heat --m t --f /path/to/batchFile.dat
ENTRYPOINT [ "python3", "./source/HEAT/source/launchHEAT.py" ]
CMD [ "--a", "0.0.0.0", "--p", "8050", "--m", "g" ]
