**********************
H.E.A.T. Documentation
**********************

Contents
--------
.. toctree::
   :maxdepth: 1

   install
   tutorial



General Information
-------------------

The Heat flux Engineering Analysis Toolkit (HEAT) is a suite of tools for predicting the heat flux
incident upon PFCs in tokamaks, and the associated PFC state (ie temperature).
The toolkit connects CAD, FVM, MHD, Plasma Physics, Visualization, HPC, and more, in one streamlined package.
The objective is to enable engineers and physicists to quickly ascertain heat loads given specific magnetic
configurations and geometric configurations.

Some examples of what HEAT can predict:
 - 3D heat loads from 2D plasmas for limited and diverted discharges
 - Heat fluxes from the optical and ion gyro orbit approximations
 - Time varying heat loads and temperature profiles
 - Magnetic field line traces
 - Many other quantities

The following physics modules are scheduled to be added to HEAT soon:
1) Photon tracing, ie radiated power using CHERAB
2) 3D plasmas using M3DC1

A HEAT paper has been published by the journal Fusion Science and Technology under
open access and can be found here: https://doi.org/10.1080/15361055.2021.1951532

For users who want to run HEAT, there are two options:
 - An appImage for running with a single executable on Linux
 - A docker container, which also allows HEAT development


To visualize HEAT results, the user will need an installation of ParaVIEW.
There has been some work to include paraview into the HEAT html interface using
paraviewweb, but this is NOT included in the releases.  More information on
ParaVIEW can be found at `<https://www.paraview.org/>`_ and ParaVIEW can be
downloaded here `<https://www.paraview.org/download/>`_.  Download version
for your operating system and follow instructions to run.


.. toctree::
   :maxdepth: 1
   :caption: Further Information:

   license
   contact


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
