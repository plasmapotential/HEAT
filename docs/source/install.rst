Installation
############

There are two ways to run HEAT as a third party user.  The first method is via
an appimage, which is basically a Linux executable.  The second method is via
a docker container, which requires docker to be installed.  Both methods provide
the user with the ability to run HEAT, but only the docker container enables
the user to edit the source code and develop new HEAT modules.

If your goal is to get up an running quickly on a Linux machine, then choose
the appImage.  If your goal is to develop and edit the source code, then choose
the docker container.

The links below provide access to tutorials both methods.

.. toctree::
   :maxdepth: 0

   appimage
   docker

To download test HEAT case
**************************
After downloading, you can test your HEAT installation by using a test case I
have prepared.  The test case can be downloaded and extracted by using the following commands
(again replacing <tag> with latest version ie vX.X.X)::

    wget https://github.com/plasmapotential/HEAT/releases/download/<tag>/testCase.tar.gz
    tar -xvzf testCase.tar.gz


Test Case Installation Video:

    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/PiU2yCTtfq0" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>
