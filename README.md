# HEAT

The Heat flux Engineering Analysis Toolkit (HEAT) is a suite of tools for predicting the heat flux
incident upon PFCs in tokamaks.  The toolkit connects CAD, FEM, MHD, Plasma Physics, Visualization,
HPC, and more, in one streamlined package.  The objective is to enable engineers and physicists to
quickly ascertain heat loads given specific magnetic configurations and geometric configurations.

If a user wants to use HEAT without building it on his/her local machine, then they can directly
access HEAT from a web browser inside the Princeton Plasma Physics Lab (PPPL) VPN.  For instructions
on where to point your web browser email Tom Looby: tlooby@vols.utk.edu

If a user wants to set up HEAT on a local machine there are a variety of configuration steps that need to be completed (ie compiling MAFOT, building OF modules from source, installing FLASK, network configuration and proxy mapping, ParaViewWeb install, etc.).  For more info, again, contact Tom.

The author engineer is Tom Looby, a PhD candidate on assignment at NSTX-U for Oak Ridge National Lab.
This project is an Oak Ridge National Lab tool built by the Fusion Energy Division, but it is also
openSource under the MIT license.

Tom's email:  tlooby@vols.utk.edu

Below is a very simple example HEAT output.  More specifically, this is a real NSTX-U castellated graphite
 tile CAD drawing from the lower outer divertor, coupled to an Eich profile physics module.  Note the complicated 3D structure of the heat flux, due to the castellated tiles.
![Alt text](exampleHF.png?raw=true "NSTXU Heat Flux onto CAD")
