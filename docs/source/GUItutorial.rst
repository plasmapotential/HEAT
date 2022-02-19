GUI Tutorial
============
The tutorials provided below are meant to guide a user through a typical HEAT
workflow in the GUI.  For your own HEAT runs, you will need to get GEQDSK and
CAD files, and then create PFC files and Input files.  The tutorials below
will provide a reference using the test case provided on github.  The install
page provide information on downloading the test case.

Running an optical heat flux simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following video provides an introduction to HEAT and guides you through running
the test case:


    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/7xDlGlWEy8g" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>


Running a field line trace
^^^^^^^^^^^^^^^^^^^^^^^^^^

The following video provides an example of running a field line trace.


    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/lrSQdDQlXR8" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>


One thing to note is that if the field line trace goes very far outside the GEQDSK,
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
