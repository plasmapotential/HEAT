INSTALLATION
============

There are two ways to run HEAT as a third party user.  The first method is via
an appimage, which is basically a Linux executable.  The second method is via
a docker container, which requires docker to be installed.  Both methods provide
the user with the ability to run HEAT, but only the docker container enables
the user to edit the source code and develop new HEAT modules.

If your goal is to get up an running quickly on a Linux machine, then choose
the appImage.  If your goal is to develop and edit the source code, then choose
the docker container.


Installing the appImage
-----------------------

HEAT appImage installation is relatively simple, and consists of downloading the latest HEAT
appImage from github (`<https://github.com/plasmapotential/HEAT>`_).  For
visualizing HEAT output, the user should have a local copy of ParaVIEW (preferably 5.8+).

To download the appImage
^^^^^^^^^^^^^^^^^^^^^^^^
Commands for installing HEAT are given below.  A video installation tutorial is provided below the commands.

In linux terminal run the following command from a directory of your choosing, with the <tag>
modified to reflect the latest version::

    wget https://github.com/plasmapotential/HEAT/releases/download/<tag>/HEAT_AppImage-<tag>-x86_64.AppImage

For example, to download <tag> = v2.0.0, the command would appear as follows:

    wget https://github.com/plasmapotential/HEAT/releases/download/v2.0.0/HEAT_AppImage-v2.0.0-x86_64.AppImage

Alternatively, you may download the latest release directly from github:
`<https://github.com/plasmapotential/HEAT/releases>`_

After downloading, make file executable::

    chmod +x HEAT_AppImage-<tag>-x86_64.AppImage


HEAT Installation Video:

    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/mDui3_z_2oM" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>


To download test HEAT case
^^^^^^^^^^^^^^^^^^^^^^^^^^
After downloading, you can test your HEAT installation by using a test case I
have prepared.  The test case can be downloaded and extracted by using the following commands
(again replacing <tag> with latest version ie vX.X.X)::

    wget https://github.com/plasmapotential/HEAT/releases/download/<tag>/testCase.tar.gz
    tar -xvzf testCase.tar.gz


Test Case Installation Video:

    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/SQXY8lI4s-o" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>




Installing the docker container
-------------------------------

To run the HEAT docker container, the user will need docker and docker-compose
installed on the local machine.

To download docker and docker-compose
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
First download docker (docker engine)
 - (`<https://docs.docker.com/engine/install/>`_)

Next, set the user up for running docker.  This includes adding the user to the
docker group (example link below for Linux)
 - (`<https://docs.docker.com/engine/install/linux-postinstall/>`_)

 Install docker-compose, which is necessary to configure the HEAT environment
  - (`<https://docs.docker.com/compose/install/>`_)

Download HEAT docker container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Now that docker is installed, you will need to pull the HEAT docker container
from docker hub.  The HEAT docker hub page is located here:
 - (`<https://hub.docker.com/r/plasmapotential/heat>`_)

To pull from docker hub, execute the following command
    docker pull plasmapotential/heat
