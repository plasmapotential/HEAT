Tutorials
========

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
