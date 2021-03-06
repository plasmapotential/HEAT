# HEAT
## Description:
The Heat flux Engineering Analysis Toolkit (HEAT) is a suite of tools for predicting the heat flux
incident upon PFCs in tokamaks, and the associated temperature distribution throughout the PFCs.  
The toolkit connects CAD, FVM, MHD, Plasma Physics, Visualization, HPC, and more, in one streamlined package.  
The objective is to enable engineers and physicists to quickly ascertain heat loads given specific magnetic
configurations and geometric configurations.

In its current form, HEAT can predict 3D heat loads from 2D plasmas for limited and diverted discharges.  
It can calculate time varying heat loads and temperature profiles.  HEAT can also be used to perform
field line traces.  In the coming year, the following modules are scheduled to be added to HEAT:
-Gyro orbit effects
-3D plasmas
-Radiated power (detachment) scenarios
-ELMs

There were three objectives taken into consideration when designing HEAT:
1) pre-shot estimation of heat loads
2) post-shot confirmation of heat loads and associated physics models
3) design optimization

The first HEAT paper is currently under review in the journal Fusion Science and Technology, and will be published under open access.  
A link will be provided when the paper is published.

An appImage is now available so that users can run HEAT on a linux machine with a single
file.  appImage is available under 'releases' and has been tested on Ubuntu18.04, Ubuntu20.04,
Centos7, Centos8.  The appImage release is still in beta version, so there are likely still bugs.
Bugs can be submitted to this github repo.

The author engineer is Tom Looby, a PhD candidate for Oak Ridge National Lab.
This project is openSource under the MIT license.

Tom's email:  tlooby@vols.utk.edu

## Installation and Tutorials
HEAT installation instructions and tutorials can be found here:
https://heat-flux-engineering-analysis-toolkit-heat.readthedocs.io/en/latest/

## Examples:
Below are a few examples of HEAT output.  HEAT produces time varying 3D heat fluxes, and can easily create visualizations leveraging the power of ParaVIEW.  There is also a HEAT presentation from Aug 2020 available [here](https://docs.google.com/presentation/d/1aqJRaxt97P6R4Kqz7xyaoegtxssHQQPuwvJgVM4cCII/edit?usp=sharing).


Example output for 30 degree section of the NSTX-U divertor with Equilibrium, Heat Flux, Temperature:
![Alt text](HF_T_EQ.gif "Example output of EQ, HF, T, video")

Example output of PFC tile temperature for various strike points sweep frequencies:
![Alt text](sideBySide.gif "Example output of EQ, HF, T, video")

Example output for limited discharges:
![Alt text](limiter.gif "Example output of EQ, HF, T, video")

HEAT Dash / plotly GUI:
![Alt text](gui1.png "HEAT DASH GUI")
![Alt text](gui2.png "HEAT DASH GUI")
