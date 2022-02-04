#Dockerfile
#Engineer: T Looby
#Date: 01/27/2022
#Description:  Dockerfile for building HEAT docker image
#
#this script should only be run on ORNL dev machine
#script is not run directly, but is used when buildDocker bash
#script is run.  buildDocker script calls this file using syntax
#similar to this:
# docker build -t plasmapotential/heat -f ./Dockerfile ./github/source

# start from base
FROM ubuntu:20.04

MAINTAINER Tom Looby <loobytp@ornl.gov>

# environment variables
ENV runMode docker
ENV AppDir /root
ENV APPDIR /root
ENV DEBIAN_FRONTEND noninteractive
ENV PYTHONPATH /usr/lib/freecad-python3/lib
ENV LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/netgen:/root/opt/mafot/lib:/root/opt/paraview/lib:/root/opt/openfoam/lib:/root/opt/openfoam/platforms/linux64GccDPInt32Opt/lib:/root/opt/openfoam/platforms/linux64GccDPInt32Opt/lib/sys-openmpi:/root/opt/swak4foam/lib:/root/opt/openfoam/platforms/linux64GccDPInt32Opt/lib/dummy
ENV PATH $PATH:/root/opt/mafot/bin:/root/opt/paraview/bin:/root/opt/openfoam/bin:/root/opt/openfoam/platforms/linux64GccDPInt32Opt/bin:/root/opt/swak4foam/bin

# install system-wide deps for HEAT
RUN apt-get -yqq update
RUN apt-get -yqq install python3
RUN apt-get -yqq install python3-pip
RUN apt-get -yqq install python3-pkg-resources
RUN apt-get -yqq install python3-distutils
#RUN apt-get -yqq install mdsplus-python
RUN apt-get -yqq install libfreecad-python3-0.18
RUN apt-get -yqq install coreutils
RUN apt-get -yqq install libnglib-6.2.1804
RUN apt-get -yqq install libcurl4
RUN apt-get -yqq install nano
RUN apt-get -yqq install git

# copy context
WORKDIR /root
COPY . .

# fetch app specific deps
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -r requirements.txt

# expose port
EXPOSE 8050

# start app
CMD [ "python3", "./source/HEAT/launchHEAT.py", "--a", "0.0.0.0", "--p", "8050", "--m", "g" ]
