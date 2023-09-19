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

# start from base
FROM ubuntu:22.04

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
ENV LD_LIBRARY_PATH $LD_LIBRARY_PATH:/root/opt/open3d/lib:/root/opt/open3d/lib/python_package

#set up paths for bins
ENV PATH $PATH:/root/opt/mafot/bin
ENV PATH $PATH:/root/opt/paraview/bin
ENV PATH $PATH:/root/opt/openfoam/bin
ENV PATH $PATH:/root/opt/openfoam/platforms/linux64GccDPInt32Opt/bin
ENV PATH $PATH:/root/opt/swak4foam/bin
ENV PATH $PATH:/root/opt/open3d/bin

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

#utility tools/libs layers
RUN apt-get -yqq install coreutils libnglib-6.2 libcurl4 libunwind-dev nano git gfortran iputils-ping iproute2

#3D plasmas layer
RUN apt-get -yqq install libopenmpi-dev libhdf5-openmpi-dev libnetcdf-dev libnetcdff-dev

#linking layer
RUN ln -s /usr/lib/x86_64-linux-gnu/libunwind.so /usr/lib/x86_64-linux-gnu/libunwind.so.1

# copy context
WORKDIR /root
COPY . .

# fetch app specific deps
RUN python3 -m pip install --upgrade pip && \
    python3 -m pip install --ignore-installed blinker && \
    python3 -m pip install -r requirements.txt

#clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/*


# expose port
EXPOSE 8050

# start app
CMD [ "python3", "./source/HEAT/launchHEAT.py", "--a", "0.0.0.0", "--p", "8050", "--m", "g" ]
