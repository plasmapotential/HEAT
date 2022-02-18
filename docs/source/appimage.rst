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


Running HEAT from the appImage
=================================
In HEAT v2+ here are two ways a user can run HEAT:
 - In an html5 based Graphical User Interface (GUI)
 - From a Terminal User Interface (TUI)

Starting HEAT from both interfaces is covered in the following sections.


Running HEAT in GUI mode
------------------------

Running HEAT in TUI mode
------------------------
