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
