TUI Tutorial
############
The tutorials provided below are meant to guide a user through a typical HEAT
workflow in the Terminal User Interface (TUI).  For your own HEAT runs, you will 
need to get GEQDSK and CAD files, and then create PFC files and Input files.  
The tutorials below will provide a reference using the test case provided on github.  
The install page provides information on downloading the test case.

Start docker container in terminal mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To run HEAT in terminal mode, the user needs to launch a bash shell inside the
docker container.  In order to ensure that the docker container maps the correct
user id and group id into the container, the ``runDockerCompose`` script is used.
This script assigns the user id (UID) and group id (GID) from the host machine
into environment variables in the container, then launches the docker container.  

To run HEAT in terminal mode, ensure the last line of ``runDockerCompose`` reads 
``docker-compose run HEAT /bin/bash``.  This command uses the docker-compose.yml
and launches a bash shell inside the container.

Mapping host directories into the container
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
One benefit of docker is that you can bind mount directories from your host
machine directly into the container.  This is achieved by editing the ``volumes``
section in the docker-compose.yml file.  The steps are as follows

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



Running the optical heat flux test case
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
HEAT ships with a series of test cases.  These test cases are used by HEAT 
developers as integration tests, but also can be used for validating a HEAT
install or getting familiar with running the code in terminal mode.  All the 
HEAT integration tests live at this location in the container: 
``/root/source/HEAT/tests/integrationTests``.  After launching a bash terminal
in the container (see section above), navigate to the HEAT source directory, 
and run a provided bash script that launches the optical heat flux test case

  .. code-block:: bash

    cd /root/source/HEAT/source
    ./runTerminalModeTestOptical

You may notice that there are other "runTerminal..." scripts in this directory.
These scripts provide a quick and convenient way to run integration tests or
to run a HEAT case that a user has bind mounted into the container.

If you would like to see what one of these integration tests look like (to use
as a template for making your own HEAT runs), you can navigate to the integration
tests directory inside the container:

  .. code-block:: bash

    cd /root/source/HEAT/tests/integrationTests

This integrationTests directory ships with the HEAT source code, and can be found
here:  <HEATdir>/tests/integrationTests

The following video provides an introduction to HEAT and guides you through running
the test case.  Note that the video is slightly dated so some of the fields in the
batchFile and input file reflect an older version of HEAT.  Users no longer need to
download the test case independently from github, nor do they need to bind mount
the test case into the container (the test case is already inside container in 
HEAT v4.0+). Nevertheless the video still provides a good overview of running HEAT
in terminal mode.

    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/VMBBddutibo" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>



Running an bind mounted HEAT case in terminal mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To run a new HEAT case in terminal mode, the user must bind bound the case into the
container.  Usually, this requires changing the "Batch mode directory" section in the 
docker-compose.yaml file:

  .. code-block:: yaml

    #      # Batch mode directory
          - <path/to/HEATrun/on/host/machine>:/root/terminal

This` will bind the ``<path/to/HEATrun/on/host/machine>`` path on the host machine to
the ``/root/terminal`` path inside the container.  After launching a bash shell in the 
container (see section above for instructions) the user can navigate to the HEAT case 
inside the container:

  .. code-block:: bash

    cd /root/terminal

The HEAT case should be located at that location if all went as expected.  In order to 
run the HEAT case after it has been mounted into the container, the user navigates 
to the HEAT source directory and runs HEAT:

  .. code-block:: bash

    cd /root/source/HEAT/source
    ./runTerminalMode

This will launch HEAT and run the HEAT case bind mounted at ``/root/terminal``.


Running a filament heat flux simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TO BE COMPLETED

Running an Elmer FEM simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Next to the previous cases, HEAT runs can also be connected with FEM simulations utilizing HEAT outputs. The test case presented here connects HEAT optical heat flux calculations with FEM temperature simulations realized with the open-source FEM solver ELMER. To run this test case, the user has to execute in the terminal (in the same directory as in the previous test cases above):

  .. code-block:: bash

    ./ runTerminalModeTestOpticalElmer

For an overview of the requirements for ELMER, the user needs to change the directoty by typing into the terminal: 

  .. code-block:: bash

    cd  /root/source/HEAT/ tests/integrationTests/nstxuTestCase

In this directory, the batch files of any test cases are provided. The batch file to execute the ELMER test case is called batchFile_optical_elmer.dat. To activate the ELMER capabilities in any HEAT batch file, the flag elmer is added under the column output. Any input required for the ELMER run must be stored in the subdirectory elmer of the respective machine directory,  here nstx. To access this folder, the user must enter the following into the terminal:

  .. code-block:: bash

    cd  /root/source/HEAT/ tests/integrationTests/nstxuTestCase/nstx/elmer

The subdirectory contains three key files: 

  1) the custom HEAT library HEATLibrary.so, which interpolates and maps the HEAT-calculated heat flux profile to the FEM mesh surface in order to obtain the temperature gradient which serves as an input to solve Fourier’s heat conduction equation.
  2) the boundary ELMER file caseTemp.sif, which summarizes any information to initiate the ELMER simulation including material properties, the partial differential equations to be solved, and the initial and all surface boundary conditions
  3)  the ELMER input file elmerFile.csv, which summarizes the respective input files and names required for the simulation (PFCname, .SIF, meshFile). PFCname specifies the name of the PFC element to be examined within the .stp file provided for the HEAT run.

A FEM mesh file is not mandatory to provide, but is recommended to incorporate FEM boundary conditions, which must be predefined during mesh generation: the FEM mesh generator incorporated into HEAT provides a generic mesh based on the provided .stp file and PFCname without dedicated boundary conditions. The ELMER output is stored in a data format compatible with HEAT and ParaView for post-processing.

Running a 3D Plasma heat flux simulation using M3DC1 output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TO BE COMPLETED
