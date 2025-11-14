#Dockerfile
#Engineer: T Looby
#Date: 03/20/2024
#Description:  Dockerfile for building HEAT docker image
#
# This dockerfile is broken into two stages.  The first stage
# clones all the source code repos and builds everything from
# source inside of a container.  This is time consuming, so
# it is made into its own stage called HEATbuilder.  HEATbuilder
# will also be saved on dockerhub from now on, with a version 
# tag.  Users can pull HEATbuilder to save time.
#
# The second stage is the final HEAT container.  It only contains
# the packages required to run HEAT, and none of the files used
# to build everything from source.  
#
# A user can build both stages like this (command also works to only
#   build final stage if HEATbuilder is already complete):
#   docker build -t plasmapotential/heat .
#
# note:  cross platform build (for apple silicon)
#   docker buildx build --platform linux/amd64,linux/arm64 -t plasmapotential/heat -f ./Dockerfile . --push
#
# User can also only build the first stage (HEATbuilder) like this:
#   docker build --target HEATbuilder -t plasmapotential/heatbuilder .
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
FROM ubuntu:22.04 AS heatbuilder
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
    apt-get -yqq install python3-pkg-resources && \
    apt-get -yqq install python3-distutils
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

#clone HEAT
RUN git clone -b v4.3 --single-branch https://github.com/plasmapotential/HEAT.git /root/source/HEAT

#change to root (home) dir
WORKDIR /root

#clone and build M3DC1 fusion-io libs
RUN git clone https://github.com/nferraro/fusion-io.git /root/source/M3DC1/fusion-io
COPY buildM3DC1 /root/source/M3DC1/fusion-io/buildM3DC1
COPY make.inc.ubuntu_x86 /root/source/M3DC1/fusion-io/install/make.inc.ubuntu_x86
RUN mkdir -p /root/source/M3DC1/build && /root/source/M3DC1/fusion-io/buildM3DC1

#clone and build MAFOT
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/hdf5/openmpi
#note branch here can be changed 
RUN git clone -b master --single-branch https://github.com/ORNL-Fusion/MAFOT.git /root/source/MAFOT
COPY buildMAFOT /root/source/MAFOT/buildMAFOT
RUN mkdir -p /root/source/MAFOT/build/bin && \
    mkdir -p /root/source/MAFOT/build/lib && \
    /root/source/MAFOT/buildMAFOT

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

#build openfoam
ENV WM_THIRD_PARTY_DIR /root/source/openfoam/ThirdParty-v${OFversion}
WORKDIR /root/source/openfoam/OpenFOAM-v${OFversion}
RUN mkdir -p /root/source/openfoam/OpenFOAM-v${OFversion}/applications/solvers/custom && \
    cp -au /root/source/HEAT/source/heatFoam /root/source/openfoam/OpenFOAM-v${OFversion}/applications/solvers/custom/heatFoam && \
    source /root/source/openfoam/OpenFOAM-v${OFversion}/etc/bashrc && \
    ./Allwmake -j -s -q -l && \
    cd /root/source/openfoam/OpenFOAM-v${OFversion}/applications/solvers/custom/heatFoam && \
    wmake

#clone and build swak4foam
WORKDIR /root/source/openfoam/
RUN source /root/source/openfoam/OpenFOAM-v${OFversion}/etc/bashrc && \
    hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam && \
    cd /root/source/openfoam/swak4Foam && \
    hg update develop && \
    ./AllwmakeAll

#clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

#======================================================================
#2nd build stage.  build container without unncesary build artifacts
#to keep final image smaller
#======================================================================
#if user doesnt want to build first stage, they can pull HEATbuilder from
#dockerhub using the following command (need to comment out first stage too):
#FROM plasmapotential/HEATbuilder:<tag> as HEATbuilder

FROM ubuntu:22.04
LABEL maintainer="tlooby@cfs.energy"
ARG DEBIAN_FRONTEND=noninteractive
ENV runMode docker

#update layer
RUN apt-get -yqq update && apt-get -yqq upgrade

#python layer
RUN apt-get -yqq install python3 && \
    apt-get -yqq install python3-pip && \
    apt-get -yqq install python3-pkg-resources && \
    apt-get -yqq install python3-distutils && \
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

#clone HEAT (note branch can be changed)
RUN git clone -b v4.3 --single-branch https://github.com/plasmapotential/HEAT.git /root/source/HEAT
WORKDIR /root/source/HEAT
RUN     git lfs install && git lfs pull
WORKDIR /root

# fetch HEAT specific python dependencies
RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install --ignore-installed blinker && \
    python3 -m pip install -r /root/source/HEAT/docker/requirements.txt && \
    python3 -m pip install -U matplotlib 


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
CMD [ "python3", "./source/HEAT/source/launchHEAT.py", "--a", "0.0.0.0", "--p", "8050", "--m", "g" ]
