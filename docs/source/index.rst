**********************
H.E.A.T. Documentation
**********************

Contents
--------
.. toctree::
   :maxdepth: 1

   install
   docker
   GUItutorial
   TUItutorial
   inputFile


General Information
-------------------
The Heat flux Engineering Analysis Toolkit (HEAT) is a suite of tools for predicting the heat flux
incident upon PFCs in tokamaks, and the associated PFC state (e.g. temperature).
The toolkit connects CAD, FVM, FEM, MHD, ray tracing, plasma physics, and more, in one streamlined package.
The objective is to enable engineers and physicists to quickly ascertain heat loads given specific magnetic
and geometric configurations. HEAT has been used to design the SPARC PFCs and continues to be developed for control room use as SPARC begins operations.

Some examples of what HEAT can predict:
 * 3D heat loads from 2D and 3D plasmas for limited and diverted discharges
 * Heat fluxes from the optical approximation, ion gyro orbit approximation, and photon flux
 * Heat and particle fluxes from filaments and runaway electrons
 * 3D heat flux profiles from RMPs and Error Fields
 * Time varying heat loads and temperature profiles
 * Magnetic field line traces
 * Many other quantities

The latest release of HEAT is **v4.2**, which includes notable additions such as:
 * Mitsuba3 for photon tracing
 * Ability to read STLs directly instead of STEPs (Bring Your Own Mesh)
 * A Runaway Electron module (A. Feyrer, MIT)
 * Ability to read arbitrary R,Z,qPar profiles from CSV (E. Tinacba, ORNL)

Full installation instructions and tutorials (including Docker) are in the **Install** and **Docker** sections below and on Read the Docs: https://heat-flux-engineering-analysis-toolkit-heat.readthedocs.io/en/latest/


Citing HEAT
-----------
To cite HEAT, use the open-access paper in Fusion Science and Technology:
https://doi.org/10.1080/15361055.2021.1951532

Other HEAT-related publications (recent first):
 * `The Experimental Validation of HEAT on the ASDEX Upgrade Tokamak <https://doi.org/10.1080/15361055.2025.2478720>`_, A. Redl et al. 2025
 * `3D modeling of n = 1 RMP driven heat fluxes on the SPARC tokamak PFCs using HEAT <https://iopscience.iop.org/article/10.1088/1741-4326/adf760/meta>`_, M. D'Abusco et al. 2025
 * `Development and validation of non-axisymmetric heat flux simulations with 3D fields using the HEAT code <https://iopscience.iop.org/article/10.1088/1741-4326/adeff1>`_, A. Wingen et al. 2025
 * `Shadow masks predictions in SPARC tokamak plasma-facing components using HEAT code and machine learning methods <https://doi.org/10.1016/j.fusengdes.2025.115010>`_, D. Corona et al. 2025
 * `HEAT simulation and IR data comparison for ST40 plasma-facing components <https://doi.org/10.1016/j.nme.2024.101791>`_, E. Tinacba et al. 2024
 * SPARC power exhaust workflows using open source tools (APS DPP Tutorial), `presentation <https://docs.google.com/presentation/d/1QsxlfUS6zo_vAwRgFoKKztvUsIzq238u2JmYD6xtum4/edit?usp=sharing>`_
 * 3D ion gyro-orbit heat load predictions for NSTX-U, Looby et al., https://iopscience.iop.org/article/10.1088/1741-4326/ac8a05
 * 3D PFC Power Exhaust Predictions for the SPARC Tokamak, Looby et al., https://meetings.aps.org/Meeting/DPP22/Session/NO03.11
 * Measurements of multiple heat flux components at the divertor target by using surface eroding thermocouples (invited), Ren et al., https://aip.scitation.org/doi/full/10.1063/5.0101719
 * SHEFT (HEAT predecessor), A. Wingen 2019, https://info.fusion.ciemat.es/OCS/EPS2019PAP/pdf/P2.1040.pdf
 * The MAFOT code (A. Wingen) is documented at https://github.com/ORNL-Fusion/MAFOT/tree/master/doc


Running HEAT
------------
To run HEAT, download the HEAT Docker container from Docker Hub (tested on Linux, macOS, and Windows). See the **Install** and **Docker** pages for details.

 * Docker Hub: https://hub.docker.com/r/plasmapotential/heat
 * HEATtools (pre/post processing): https://github.com/plasmapotential/HEATtools

This repository uses **Git LFS** for large data files (e.g. the 3D fields test case). After cloning, run::

    git lfs install
    git lfs pull

before running the 3D fields test case (e.g. ``./runTerminalModeTest3Dfields``).


Developer and license
---------------------
The developer is Tom Looby, Scientist at Commonwealth Fusion Systems.
Contact: tlooby@cfs.energy

This project is open source under the **MIT license**. Contributions are welcome (documentation, code, issues); see the GitHub repository for more.


Visualization (ParaVIEW)
------------------------
To visualize HEAT results you need ParaView. HEAT can produce time-varying 3D heat fluxes and visualizations that work with ParaView.

 * ParaView: https://www.paraview.org/
 * Download: https://www.paraview.org/download/

Note: Work has been done to integrate ParaviewWeb into the HEAT HTML interface; this is **not** included in the current releases.


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
