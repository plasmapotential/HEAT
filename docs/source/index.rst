**********************
H.E.A.T. Documentation
**********************

Contents
--------
.. toctree::
  install



General Information
-------------------

The Heat flux Engineering Analysis Toolkit (HEAT) is a suite of tools for predicting the heat flux
incident upon PFCs in tokamaks, and the associated temperature distribution throughout the PFCs.
The toolkit connects CAD, FVM, MHD, Plasma Physics, Visualization, HPC, and more, in one streamlined package.
The objective is to enable engineers and physicists to quickly ascertain heat loads given specific magnetic
configurations and geometric configurations.

In its current form, HEAT can predict 3D heat loads from 2D plasmas for limited and diverted discharges.
It can calculate time varying heat loads and temperature profiles.  HEAT can also be used to perform
field line traces.  In the coming year, the following modules are scheduled to be added to HEAT:
  * Gyro orbit effects
  * 3D plasmas
  * Radiated power (detachment) scenarios
  * ELMs

There were three objectives taken into consideration when designing HEAT:
  1) pre-shot estimation of heat loads and temperature
  2) post-shot confirmation of heat loads and associated physics models
  3) design optimization

HEAT has been compiled into an appImage.  For appImage documentation see `<https://appimage.org/>`_ .
The appImage allows HEAT to be deployed as a single file, for use across linux distros.
The user sets the appImage to executable, then runs it as if it were a binary file.
The appImage is available under 'releases' on github and has been tested on Ubuntu18.04, Ubuntu20.04,
Centos7, Centos8.  The appImage release is still in beta version, so there are likely still bugs.
Bugs can be submitted to the github repo.  See Installation section for more info on running
HEAT with the appImage.

Running HEAT requires the following (recommended version in parentheses):
  * Linux machine (Ubuntu18-20,Centos7-8)
  * Web Browser (Google Chrome)

To visualize HEAT results, the user will need an installation of ParaVIEW.
There has been some work to include paraview into the HEAT html interface using
paraviewweb, but this is NOT included in the appImage.  More information on
ParaVIEW can be found at `<https://www.paraview.org/>`_ and ParaVIEW can be
downloaded here `<https://www.paraview.org/download/>`_.


.. toctree::
   :maxdepth: 1
   :caption: Further Information:

   license
   contact





Installation
------------

HEAT installation is relatively simple, and consists of downloading the latest HEAT
appImage from github (`<https://github.com/plasmapotential/HEAT>`_).  For
visualizing HEAT output, the user should have a local copy of ParaVIEW (preferably 5.8+).

To download the appImage
^^^^^^^^^^^^^^^^^^^^^^^^
Commands for installing HEAT are given below.  A video installation tutorial is provided below the commands.

In linux terminal run the following command from a directory of your choosing::

    wget https://github.com/plasmapotential/HEAT/releases/download/v1.2-beta/HEAT_AppImage-v1.2-beta-x86_64.AppImage

Alternatively, you may download the latest release directly from github:
`<https://github.com/plasmapotential/HEAT/releases>`_

After downloading, make file executable::

    chmod +x HEAT_AppImage-v1.2-beta-x86_64.AppImage


HEAT Installation Video:

    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/mDui3_z_2oM" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>


To download test HEAT case
^^^^^^^^^^^^^^^^^^^^^^^^^^
After downloading, you can test your HEAT installation by using a test case we
have prepared.  The test case can be downloaded and extracted by using the following commands::

    wget https://github.com/plasmapotential/HEAT/releases/download/v1.2-beta/testRun.tar.gz
    tar -xvzf testRun.tar.gz


Test Case Installation Video:

    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/SQXY8lI4s-o" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>


To download ParaVIEW
^^^^^^^^^^^^^^^^^^^^
In a web browser, go to: `<https://www.paraview.org/download/>`_.  Download version
for your operating system and follow instructions to run.



Tutorials
---------
Accessing the HEAT graphical user interface
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following video provides an introduction to HEAT and guides you through running
the test case:


    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/ezdfC84ANTk" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>





After downloading appImage and making it executable, run the appImage to start
the HEAT server::

    ./HEAT_AppImage-v1.2-beta-x86_64.AppImage

By default, this will start HEAT and launch a Flask server on the localhost at port 8050.
To access the GUI, open a web browser and navigate to 127.0.0.1:8050::

    google-chrome 127.0.0.1:8050

You should now see the HEAT user interface.  Be sure to install the testRun folder as described in the Installation
section beforehand. For the test case, the files in the testRun folder can be used in the following workflow:

  1) Select a machine
  2) Drag NSTXU_TEST_input.csv into "Select input file" box
  3) Drag gFile g204118.000004 into "Select gFiles" box
  4) Click "Load MHD" button
  5) Click "Load res settings" button to load CAD resolution
  6) Drag STP file IBDH_2tiles.step into "Select STP file" box
  7) Drag PFC file IBDH_2tiles.csv into "PFC Settings > Select Files"
  8) Choose heat flux setting "Gaussian Spreading".  HEAT will default to Eich profile settings
  9) Click "Load HF" button
  10) Click "Load OF Settings" button
  11) Click on "Run HEAT" tab in left ribbon
  12) Select "Heat flux point cloud"
  13) Select "openFOAM thermal analysis"
  14) Click "Run HEAT" button
  15) To watch logFiles, either watch original linux terminal or click "LogFile" tab from left ribbon
  16) All results are saved into directory $HOME/HEAT/data, where $HOME is user linux home directory
  17) For heat flux results, open paraview and import: $HOME/HEAT/data/nstx_204118/000004/paraview/HeatFlux_all.vtk
  18) For temperature results, open paraview and import: $HOME/HEAT/data/nstx_204118/openFoam/heatFoam/SOLID843/SOLID843.foam

A series of videos will be added soon to describe this process and others.



Running a field line trace
^^^^^^^^^^^^^^^^^^^^^^^^^^

The following video provides an example of running a field line trace.


    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/lrSQdDQlXR8" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>


After downloading appImage and making it executable, run the appImage to start
the HEAT server::

    ./HEAT_AppImage-v1.2-beta-x86_64.AppImage

By default, this will start HEAT and launch a Flask server on the localhost at port 8050.
To access the GUI, open a web browser and navigate to 127.0.0.1:8050::

    google-chrome 127.0.0.1:8050

You should now see the HEAT user interface.  You will need two things: a GEQDSK file and the (x,y,z) coordinates of the point you want to launch a trace from.


In the HEAT user interface:

  1) Select a machine
  2) Enter shot number corresponding to GEQDSK file
  3) Enter time bounds (minimum/maximum timestep [ms]) that include the timestep of your GEQDSK file
  4) Enter 0 in Number of Trace Steps (degrees)
  5) Drag your GEQDSK file from a file browser into the Drag and Drop box
  6) Ensure 2D Plasmas is selected
  7) Click Load MHD

Now that you have loaded the MHD into HEAT, you need to set up the field line trace and run it.

  1) Click on the Run HEAT tab
  2) Unselect all checkboxes
  3) Select checkbox for B-field Trace
  4) Enter coordinates to launch field line trace from
  5) Enter ionDirection (1 for positive toroidal direction, -1 for negative toroidal direction, 0 for both)
  6) Enter number of degrees you want to trace for
  7) Click RUN HEAT button
  8) Data will be saved in Data Directory (usually $HOME/HEAT/data/<machine>_<shot>/<timestep>)

One thing to note is that if the field line trace goes very far outside the machine,
or if a GS integration step fails for some reason, then the field line trace will end at that point.

If you get field line traces that look wonky and go in crazy directions, then your GEQDSK file might
need to be adjusted.  HEAT has a suite of tools for dealing with GEQDSK files on the gFile Tools tab.
Email Tom for more information on using HEAT to fix gFiles.




Running HEAT on a Local Area Network (LAN)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The default HEAT settings enable the user to access the GUI in a web browser at
the localhost via port 8050 (127.0.0.1:8050).  When HEAT is being utilized on a
local machine, this is ideal, but when the user desires to run HEAT as a server
(such as on a cluster) and connect to it from a client machine across a network,
then the IP address and port number must be included as command line arguments.

This method can be useful for users who do not have linux operating systems
(ie macOS, windows, android), or for users who want to use HEAT without
downloading and configuring on their local machine.  The limitations to this method
are:
    1) Only one user can access HEAT at a time.  There is work underway to enable
       multiple user sessions, but HEAT can only serve a single session currently.
    2) The IP address assigned to HEAT must be the IP address assigned to that
       machine's network interface card (NIC) by the network DHCP server.  The
       user cannot assign a random IP address to HEAT.
    3) The port must allow TCP/IP traffic.

To do this, include the IP address and desired port number as
switches when executing the HEAT appImage command in a terminal::

    ./HEAT_AppImage-v1.2-beta-x86_64.AppImage -a <address> -p <port>

The -a switch precedes the intended IP address and the -p switch precedes the
port number.  For example, to run HEAT at 192.168.0.100 on port 7500, the following command would
be used::
    ./HEAT_AppImage-v1.2-beta-x86_64.AppImage -a 192.168.0.100 -p 7500

Then, on the client machine, HEAT can be accessed by opening a web browser and
navigating to 192.168.0.100:7500


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
