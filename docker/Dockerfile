#Dockerfile
#Engineer: T Looby
#Date: 01/27/2022
#Description:  Dockerfile for building HEAT docker image
#
#this script should only be run on dev machine
#script is not run directly, but is used when buildDocker bash
#script is run.  buildDocker script calls this file using syntax
#similar to this:
# docker build -t plasmapotential/heat -f ./Dockerfile ./github/source
#
#cross platform build (for apple silicon)
#docker buildx build --platform linux/amd64,linux/arm64 -t plasmapotential/heat -f ./Dockerfile . --push

#after creating image, tag it using the following command
# docker tag plasmapotential/heat:latest plasmapotential/heat:<tag>
#
#upload to dockerhub via the following command
# docker push plasmapotential/heat:<tag>
#

# start from base
FROM ubuntu:22.04
SHELL ["/bin/bash", "-c"] 

LABEL maintainer="tlooby@cfs.energy"

# environment variables
ENV runMode docker
ENV AppDir /root
ENV APPDIR /root
ENV DEBIAN_FRONTEND noninteractive
ENV PYTHONPATH /usr/lib/freecad-python3/lib

#set up paths for libs
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/qt5
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/netgen
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/mafot/lib
#ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/paraview/lib #comment to prevent qt5 conflicts
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/openfoam/lib
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/openfoam/platforms/linux64GccDPInt32Opt/lib
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/openfoam/platforms/linux64GccDPInt32Opt/lib/sys-openmpi
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/swak4foam/lib
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/openfoam/platforms/linux64GccDPInt32Opt/lib/dummy
#ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/open3d/lib:/root/opt/open3d/lib/python_package
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/m3dc1/lib

#set up paths for bins
ENV PATH $PATH:/root/opt/mafot/bin
ENV PATH $PATH:/root/opt/paraview/bin
ENV PATH $PATH:/root/opt/openfoam/bin
ENV PATH $PATH:/root/opt/openfoam/platforms/linux64GccDPInt32Opt/bin
ENV PATH $PATH:/root/opt/swak4foam/bin
#ENV PATH $PATH:/root/opt/open3d/bin
ENV PATH $PATH:/root/opt/m3dc1/bin

#===system-wide deps for HEAT
#update layer
RUN apt-get -yqq update

#python layer
RUN apt-get -yqq install python3 && \
    apt-get -yqq install python3-pip && \
    apt-get -yqq install python3-pkg-resources && \
    apt-get -yqq install python3-distutils
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

# elmer layer
RUN apt-add-repository ppa:elmer-csc-ubuntu/elmer-csc-ppa && \
    apt-get update && \
    apt-get install -yqq elmerfem-csc

#utility tools/libs layers
#RUN apt-get -yqq libnglib-6.2 
RUN apt-get -yqq install coreutils libcurl4 libunwind-dev nano git gfortran iputils-ping iproute2 htop

#3D plasmas layer
RUN apt-get -yqq install libopenmpi-dev libhdf5-openmpi-dev libnetcdf-dev libnetcdff-dev

#linking layer
RUN ln -s /usr/lib/x86_64-linux-gnu/libunwind.so /usr/lib/x86_64-linux-gnu/libunwind.so.1

#clone EFIT
RUN git clone https://github.com/ORNL-Fusion/EFIT.git /root/source/EFIT

#clone HEAT
RUN git clone https://github.com/plasmapotential/HEAT.git /root/source/HEAT

#change to root (home) dir
WORKDIR /root

#clone and build M3DC1 fusion-io libs
RUN git clone https://github.com/nferraro/fusion-io.git /root/source/M3DC1/fusion-io
COPY buildM3DC1 /root/source/M3DC1/fusion-io/buildM3DC1
COPY make.inc.ubuntu_x86 /root/source/M3DC1/fusion-io/install/make.inc.ubuntu_x86
RUN mkdir -p /root/source/M3DC1/build && /root/source/M3DC1/fusion-io/buildM3DC1

#clone and build MAFOT
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/hdf5/openmpi
RUN git clone -b HEATv4.1 --single-branch https://github.com/ORNL-Fusion/MAFOT.git /root/source/MAFOT
COPY buildMAFOT /root/source/mafot/buildMAFOT
RUN /root/source/mafot/buildMAFOT

#clone openfoam
RUN git clone https://develop.openfoam.com/Development/openfoam.git /root/source/openfoam/OpenFOAM-v2312 && \
    git clone https://develop.openfoam.com/Development/ThirdParty-common.git /root/source/openfoam/ThirdParty-v2312 && \
    apt-get install -yqq mercurial && \
    hg clone http://hg.code.sf.net/p/openfoam-extend/swak4Foam /root/source/openfoam/swak4Foam && \
    mkdir /root/source/openfoam/OpenFOAM-v2312/applications/solvers/custom && \
    cp -au /root/source/HEAT/source/heatFoam /root/source/openfoam/OpenFOAM-v2312/applications/solvers/custom/heatFoam

#build openfoam
ENV SHELL /bin/bash
ENV WM_THIRD_PARTY_DIR /root/source/openfoam/ThirdParty-v2312
WORKDIR /root/source/openfoam/OpenFOAM-v2312
RUN apt-get -yqq update && \
    apt-get install -yqq python-is-python3 flex bison && \
    source /root/source/openfoam/OpenFOAM-v2312/etc/bashrc && \
#    foamSystemCheck && \
#    foam 
    ./Allwmake -j -s -q -l 
#    cd /root/source/openfoam/swak4Foam && \
#    hg update develop && \
#    ./Allwmake



#COPY . .
#
# fetch app specific deps
#RUN python3 -m pip install --upgrade pip && \
#    python3 -m pip install --ignore-installed blinker && \
#    python3 -m pip install -r /root/source/HEAT/docker/requirements.txt



#clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/*


# expose port
EXPOSE 8050

WORKDIR /root

# start app
CMD [ "python3", "./source/HEAT/source/launchHEAT.py", "--a", "0.0.0.0", "--p", "8050", "--m", "g" ]
