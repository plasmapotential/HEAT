# HEAT

The Heat flux Engineering Analysis Toolkit (HEAT) is a suite of tools for predicting the heat flux
incident upon PFCs in tokamaks.  The toolkit connects CAD, FEM, MHD, Plasma Physics, Visualization,
HPC, and more, in one streamlined package.  The objective is to enable engineers and physicists to
quickly ascertain heat loads given specific magnetic configurations and geometric configurations.

Note that before running HEAT there are a variety of configuration steps that need to be completed
on the local machine (ie compiling MAFOT, building OF modules from source, installing FLASK,
network configuration and proxy mapping, etc.).  

For more info contact Tom Looby, a PhD candidate on assignment at NSTX-U for Oak Ridge National Lab.
This project is an Oak Ridge National Lab tool, but it is openSource under the MIT license.

email:  tlooby@vols.utk.edu

Below is an example HEAT output.  More specifically, this is a real NSTX-U castellated graphite
 tile CAD drawing from the lower outer divertor, coupled to an Eich profile physics module.  
 Note the complicated 3D structure of the heat flux, due to the castellated tiles.
![Alt text](exampleHF.png?raw=true "NSTXU Heat Flux onto CAD")
