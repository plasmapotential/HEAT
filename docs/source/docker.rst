HEAT via the docker container
#############################
This page provides information on downloading and running HEAT from the docker
container.

Installing the docker container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

To pull from docker hub, execute the following command::

    docker pull plasmapotential/heat:<tag>

where <tag> reflects the latest HEAT version (ie v3.0, v4.1, or whatever version you want)

Download HEAT source code from github
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run HEAT using docker, it can be useful to have several files from the HEAT
github page that set up the docker environment.  The easiest way to download the
HEAT source code is to create a new directory, and pull the source using git::

    cd <sourcePath>
    git clone https://github.com/plasmapotential/HEAT.git

Where <sourcePath> is the path where you want to save HEAT  Once you have the
HEAT source code downloaded, the files you need to run docker are located in the
docker directory, <sourcePath>/docker

If you already have the HEAT source code downloaded, then you can pull the latest with:

.. code-block:: bash

    git pull

If you want to force the pull to overwrite your local changes:

.. code-block:: bash

    git reset --hard HEAD


Starting HEAT with docker
^^^^^^^^^^^^^^^^^^^^^^^^^
In HEAT v2+ here are two ways a user can run HEAT:
 - In an html5 based Graphical User Interface (GUI)
 - From a Terminal User Interface (TUI)

Starting HEAT from both interfaces is covered in the following sections.  For
both modes, the user needs docker, docker-compose, and the HEAT source code
installed.  Additionally, the following video provides an introduction to using
HEAT from the docker container:


    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/ygNJRAYitAI" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>


One-time setup: docker/.env
^^^^^^^^^^^^^^^^^^^^^^^^^^^
The docker compose recipe (docker/docker-compose.yml) is parameterized via an
environment file, docker/.env.  Every variable has a sensible default, so this
step is optional, but running it once is recommended::

    ./docker/setup.sh

This generates docker/.env from docker/.env.example and records:

 - ``dockerUID`` / ``dockerGID`` — your user/group IDs, passed into the
   container so files HEAT writes are owned by you (see Permissions below)
 - ``HEAT_DATA_DIR`` — the host directory mounted at /root/HEAT for HEAT
   output (default ``~/HEAT``)
 - ``HEAT_IMAGE_TAG`` — the HEAT image version to run
 - ``HEAT_PORT`` — the host port for the GUI (default 8050)
 - ``HEAT_GPU`` — set to 1 to enable the NVIDIA GPU reservation (requires the
   nvidia container runtime)
 - ``HEAT_RUNS_DIR`` — optional (left empty by setup.sh): a host directory
   mounted at /root/HEAT_runs in every container, so input files can reference
   shared inputs by absolute ``/root/HEAT_runs/...`` paths (see the TUI mode
   section below)

Edit docker/.env at any time to change these.  A variable exported in your
shell overrides .env for one-off runs, e.g.
``HEAT_IMAGE_TAG=v4.2.7 ./run.sh gui``.

Developers: ``./docker/setup.sh --dev`` additionally writes
docker/docker-compose.override.yml, which mounts your HEAT checkout over the
source baked into the image, so local code changes run without rebuilding.
Delete that file to revert to the image's built-in source.

Running HEAT: the run.sh wrapper
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A wrapper script at the repository root, run.sh, is the recommended way to
launch every HEAT mode.  It validates your configuration before starting a
container and prints actionable errors::

    ./run.sh gui                     # GUI on http://localhost:8050 (Ctrl-C to stop)
    ./run.sh gui -d                  # same, detached
    ./run.sh tui path/to/batchFile.dat   # terminal (batch) mode
    ./run.sh shell                   # interactive bash shell in the container
    ./run.sh test optical            # run an integration test case (mirrors CI)
    ./run.sh fix_permissions         # chown the HEAT data dir back to your user
    ./run.sh cleanup                 # remove orphaned/stopped containers
    ./run.sh rm_docker               # remove all project images/containers/volumes
    ./run.sh help                    # full command list

Start HEAT in GUI mode
^^^^^^^^^^^^^^^^^^^^^^
The image uses ENTRYPOINT + CMD so the GUI is the default mode::

    ./run.sh gui

or, with raw compose from the docker directory::

    docker compose up

Then open http://localhost:8050 in a browser.

Start HEAT in TUI mode
^^^^^^^^^^^^^^^^^^^^^^
Point run.sh at your batchFile::

    ./run.sh tui <batchModePath>/batchFile.dat

HEAT runs in terminal mode and the container is removed when the run finishes.
No editing of docker-compose.yml is required.

Alternatively, get an interactive shell and run HEAT manually.  Because the
image uses ENTRYPOINT to run HEAT, the entrypoint must be overridden to get a
shell (``./run.sh shell`` does this for you)::

    docker compose run --rm --entrypoint "" HEAT /bin/bash

then, inside the container::

    cd /root/source/HEAT/
    python3 ./source/launchHEAT.py --m t --f /root/terminal/batchFile.dat

How file paths resolve in TUI mode
""""""""""""""""""""""""""""""""""
The **parent directory of the batchFile** is bind-mounted at /root/terminal
inside the container for that run.  The files named in the batchFile columns
(GEQDSK, CAD, PFC, Input) are written as bare filenames; HEAT resolves them
inside the ``<MachFlag>/`` subdirectory of that folder, where MachFlag is the
batchFile's first column — so they must live at
``<batchFile folder>/<MachFlag>/<file>``.

Paths written *inside* an input file (for example ``radFile``) are
different: HEAT opens them verbatim inside the container, so they must be
written as container paths.  Two styles work:

 - ``/root/terminal/...`` — resolves inside the batchFile's own folder via
   the per-run mount above.  No configuration needed, but the referenced
   file must travel with the run folder.
 - ``/root/HEAT_runs/...`` — resolves inside the stable mount configured by
   ``HEAT_RUNS_DIR`` in docker/.env.  This mount is identical for every run,
   so it suits inputs shared between run folders; if ``HEAT_RUNS_DIR`` is
   not set, the fallback mount is an empty volume, so such paths fail with
   file-not-found.

Example: a machine named myTokamak, whose batchFile rows therefore use
``MachFlag = myTokamak`` (the subdirectory must be named exactly after the
MachFlag, which in a real run is one of the supported machine flags —
``sparc``, ``arc``, ``d3d``, ``nstx``, ...).  With
``HEAT_RUNS_DIR=/home/me/HEAT_runs`` in docker/.env and this layout on the
host::

    /home/me/HEAT_runs/myRun/batchFile.dat
    /home/me/HEAT_runs/myRun/myTokamak/run_input.csv
    /home/me/HEAT_runs/myRun/myTokamak/radsource.csv

``./run.sh tui /home/me/HEAT_runs/myRun/batchFile.dat`` mounts
``/home/me/HEAT_runs/myRun`` at /root/terminal.  The batchFile's Input
column holds the bare name ``run_input.csv`` (resolved to
``/root/terminal/myTokamak/run_input.csv``), and run_input.csv may reference
the radiation file as either ``/root/terminal/myTokamak/radsource.csv`` or
``/root/HEAT_runs/myRun/myTokamak/radsource.csv`` — both name the same host
file.

These path conventions apply to TUI mode only.  In GUI mode every input
(including the radiation source file) is uploaded through the browser, so
container paths written inside an input file are not used; the mounts are
still present in the GUI container but the GUI does not read from them.

Permissions in Docker on Linux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Applications in a docker container run as root, so files HEAT writes to the
mounted data directory would be owned by root on the host.  To avoid this,
HEAT chowns its output to the IDs given by the ``dockerUID`` / ``dockerGID``
environment variables.  ``./docker/setup.sh`` records your IDs in docker/.env
(preferring the ``docker`` group's GID when it exists); run.sh also fills them
in automatically at launch when they are not set anywhere.

If you end up with root-owned files anyway (e.g. from a raw ``docker run``
without these variables), reclaim them with::

    ./run.sh fix_permissions

