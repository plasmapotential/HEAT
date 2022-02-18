HEAT via the docker container
=============================
This page provides information on downloading and running HEAT from the docker
container.

Installing the docker container
===============================

To run the HEAT docker container, the user will need docker and docker-compose
installed on the local machine.

To download docker and docker-compose
-------------------------------------

First download docker (docker engine)
 - (`<https://docs.docker.com/engine/install/>`_)

Next, set the user up for running docker.  This includes adding the user to the
docker group (example link below for Linux)
 - (`<https://docs.docker.com/engine/install/linux-postinstall/>`_)

Install docker-compose, which is necessary to configure the HEAT environment
  - (`<https://docs.docker.com/compose/install/>`_)

Download HEAT docker container
------------------------------

Now that docker is installed, you will need to pull the HEAT docker container
from docker hub.  The HEAT docker hub page is located here:
 - (`<https://hub.docker.com/r/plasmapotential/heat>`_)

To pull from docker hub, execute the following command::

    docker pull plasmapotential/heat

Download HEAT source code from github
-------------------------------------

To run HEAT using docker, you will need to have several files from the HEAT
github page that set up the docker environment.  The easiest way to download the
HEAT source code is to create a new directory, and pull the source using git::

    cd <sourcePath>
    git clone https://github.com/plasmapotential/HEAT.git

Where <sourcePath> is the path where you want to save HEAT  Once you have the
HEAT source code downloaded, the files you need to run docker are located in the
docker directory, <sourcePath>/docker


Starting HEAT from docker
=========================
In HEAT v2+ here are two ways a user can run HEAT:
 - In an html5 based Graphical User Interface (GUI)
 - From a Terminal User Interface (TUI)

Starting HEAT from both interfaces is covered in the following sections.  For
both modes, the user needs docker, docker-compose, and the HEAT source code
installed.

Start HEAT in GUI mode
------------------------
To start HEAT using the graphical user interface, perform the following steps:

  1) Navigate to the HEAT source code docker directory, <sourcePath>/docker
  2) Once in the docker directory, make sure the last 4 lines appear as follows::

      #run docker compose
      docker-compose up
      #run docker compose interactively (for terminal mode)
      #docker-compose run HEAT /bin/bash
  3) Run docker compose from within the docker directory::
      docker-compose up

Start HEAT in TUI mode
------------------------
To start HEAT using the terminal user interface, perform the following steps:

  1) Navigate to the HEAT source code docker directory, <sourcePath>/docker
  2) Edit the docker-compose.yml recipe file.  Under the volumes section,
     the user can bind directories on their local host machine into the docker
     container.  For each of these lines, the host path and container path are
     in the following format:

        <hostPath>:<containerPath>
     You should not need to edit the <containerPath>, but you will need to edit
     the <hostPath>.  For example, to bind the HEAT source code that you
     downloaded from github at the path <sourcePath> into the container, you
     would have the following line under volumes in the recipe::

          - <sourcePath>:/root/source/HEAT
     You should uncomment the lines that correspond to the local packages that
     you have installed.  The HEAT data directory should always be uncommented
     and binded::

          - ${HOME}/HEAT:/root/HEAT
     For running in terminal mode, you will need to uncomment the line that
     binds your local batchMode directory into the container::

          - <batchModePath>:/root/terminal
     where <batchModePath> is the directory where your batchFile lives.

  3) In the docker directory, make sure the last line appears as follows::

      docker-compose run HEAT /bin/bash
  4) Run docker compose from within the docker directory::

      docker-compose up
  5) Running docker-compose in terminal mode launches a bash terminal inside the
     container.  Once inside the container, navigate to the HEAT source code
     directory::

      cd /root/source/HEAT/
  6) Once in the source directory, launch HEAT using the batchFile.dat that
     was binded into the container in step 2)::

      python3 launchHEAT.py --m t --f /root/terminal/batchFile.dat
