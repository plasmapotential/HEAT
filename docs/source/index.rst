**********************
H.E.A.T. Documentation
**********************

Contents
--------
.. toctree::
   :maxdepth: 1

   install
   GUItutorial
   TUItutorial
   inputFile



General Information
-------------------
The Heat flux Engineering Analysis Toolkit (HEAT) is a suite of tools for predicting the heat flux
incident upon PFCs in tokamaks, and the associated PFC state (ie temperature, stress, etc).
The toolkit connects CAD, FVM, FEM, MHD, ray tracing, plasma physics, and more, in one streamlined package.
The objective is to enable engineers and physicists to quickly calculate heat loads given specific magnetic
configurations and geometric configurations.

Some examples of what HEAT can predict:
 * 3D heat loads from 2D and 3D plasmas for limited and diverted discharges
 * Heat fluxes from the optical approximation, ion gyro orbit approximation, and photon flux
 * Heat and particle fluxes from filaments
 * 3D heat flux profiles from RMPs and Error Fields
 * Temperature, Stress, Recrystallization Kinetics
 * Time varying heat fluxes and thermal state
 * Magnetic field line traces
 * Many other quantities


The latest release of HEAT is v4.2.3.

To cite HEAT, you can use a paper published by the journal Fusion Science and Technology under open access.  The paper can be found here: https://doi.org/10.1080/15361055.2021.1951532

Other recent HEAT related publications:
 * 3D ion gyro-orbit heat load predictions for NSTX-U, Looby et al, https://iopscience.iop.org/article/10.1088/1741-4326/ac8a05
 * 3D PFC Power Exhaust Predictions for the SPARC Tokamak, Looby et al, https://meetings.aps.org/Meeting/DPP22/Session/NO03.11
 * Measurements of multiple heat flux components at the divertor target by using surface eroding thermocouples (invited), Ren et al, https://aip.scitation.org/doi/full/10.1063/5.0101719

For users who want to run HEAT, you will need to download the HEAT docker container from dockerhub.
Dockerhub repo:  https://hub.docker.com/r/plasmapotential/heat
There is a companion repo to this one, which provides some HEAT pre/post processing functions:
https://github.com/plasmapotential/HEATtools.git

The developer is Tom Looby, Senior Scientist at Commonwealth Fusion Systems.

This project is open source under the MIT license.

Tom's email:  tlooby@cfs.energy

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
