HEAT via the appImage
=====================
This page provides information on downloading and running HEAT from the appImage.

Installing the appImage
=======================

HEAT appImage installation is relatively simple, and consists of downloading the latest HEAT
appImage from github (`<https://github.com/plasmapotential/HEAT>`_).  For
visualizing HEAT output, the user should have a local copy of ParaVIEW (preferably 5.8+).

To download the appImage
------------------------

Commands for installing HEAT are given below.  A video installation tutorial is provided below the commands.

In linux terminal run the following command from a directory of your choosing, with the <tag>
modified to reflect the latest version::

    wget https://github.com/plasmapotential/HEAT/releases/download/<tag>/HEAT_AppImage-<tag>-x86_64.AppImage

For example, to download <tag> = v2.0.0, the command would appear as follows::

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


Running HEAT from the appImage
=================================
In HEAT v2+ here are two ways a user can run HEAT:
 - In an html5 based Graphical User Interface (GUI)
 - From a Terminal User Interface (TUI)

Starting HEAT from both interfaces is covered in the following sections.


Start HEAT in GUI mode
------------------------
After downloading appImage and making it executable, run the appImage to start the HEAT server,
where <tag> is replaced with the version you have (ie v2.0.0):

    ./HEAT_AppImage-<tag>-x86_64.AppImage

By default, this will start HEAT and launch a Flask server on the localhost at port 8050.
To access the GUI, open a web browser and navigate to 127.0.0.1:8050::

    google-chrome 127.0.0.1:8050

You should now see the HEAT user interface.  See the tutorial using the test
case for information on how to run HEAT from the GUI.
Test Case Tutorial:
:doc:`testCase`


Running HEAT in TUI mode
------------------------
After downloading the appImage and making it executable, you can run in TUI mode
by following these steps:

  1) Create a batch mode directory, as described in the test case tutorials :doc:`testCase`
  2) Run the following command from the terminal, replacing the path to the
     batchFile as appropriate:

      ./HEAT_AppImage-v2.0.0-x86_64.AppImage -m t -f /path/to/batchFile.dat
