# HEAT
## Description:
The Heat flux Engineering Analysis Toolkit (HEAT) is a suite of tools for predicting the heat flux
incident upon PFCs in tokamaks, and the associated temperature distribution throughout the PFCs.  
The toolkit connects CAD, FVM, MHD, Plasma Physics, Visualization, HPC, and more, in one streamlined package.  
The objective is to enable engineers and physicists to quickly ascertain heat loads given specific magnetic
configurations and geometric configurations.

In its current form, HEAT can predict 3D heat loads from 2D plasmas for limited and diverted discharges.  
It can calculate time varying heat loads and temperature profiles.  HEAT can also be used to perform
field line traces.  The most recent HEAT module to be added to the toolkit is ion gyro orbit
heat load predictions.

The following physics modules are scheduled to be added to HEAT soon:
1) 3D plasmas
2) Radiated power (detachment) scenarios
3) ELMs

There were three objectives taken into consideration when designing HEAT:
1) pre-shot estimation of heat loads
2) post-shot confirmation of heat loads and associated physics models
3) design optimization

The first HEAT paper has been published by the journal Fusion Science and Technology.  The paper is available to the public under open access and can be found here: https://doi.org/10.1080/15361055.2021.1951532

An appImage is available so that users can run HEAT on a linux machine with a single
file.  appImage is available under 'releases' and has been tested on Ubuntu18.04, Ubuntu20.04,
Centos7, Centos8.  The appImage release is still in beta version, so there are likely still bugs.
Bugs can be submitted to this github repo.

The author engineer is Tom Looby,  a Postdoctoral Researcher at Oak Ridge National Laboratory.  The project began during Tom's PhD.  This project is open source under the MIT license.

Tom's email:  loobytp@ornl.gov

## Recent Updates
The latest HEAT appImage is available under the "Releases" page.  Some recent updates, bug fixes, and user feature requests include:
 - Ion gyro orbit heat loads
 - Ion gyro orbit trajectory tracing
 - Ion gyro orbit plots on Outputs tab
 - Multiple material selections for OpenFOAM
 - Updated OpenFOAM to v2106 from v1912

## Installation and Tutorials
HEAT installation instructions and tutorials can be found here:
https://heat-flux-engineering-analysis-toolkit-heat.readthedocs.io/en/latest/

Note that the videos are from an earlier version of HEAT, and do not reflect the most recent appImage perfectly.  They will be updated in early 2022.

## Examples:
Below are a few examples of HEAT output.  HEAT produces time varying 3D heat fluxes, and can easily create visualizations leveraging the power of ParaVIEW.  

There is a HEAT presentation from Aug 2020 available [here](https://docs.google.com/presentation/d/1aqJRaxt97P6R4Kqz7xyaoegtxssHQQPuwvJgVM4cCII/edit?usp=sharing)

And a presentation from Dec 2021 available [here](https://docs.google.com/presentation/d/1BF2DvYyuPM_ATutrNDVy_r3_vKbj0a8H2UtDaoGvVg8/edit?usp=sharing)

Example output for 30 degree section of the NSTX-U divertor with Equilibrium, Heat Flux, Temperature:
![Alt text](HF_T_EQ.gif "Example output of EQ, HF, T, video")

Example output of PFC tile temperature for various strike points sweep frequencies:
![Alt text](sideBySide.gif "Example output of EQ, HF, T, video")

Example trace for ion gyro orbit tracing:
![Alt text](gyroTrace.png "Example ion gyro orbit trajectory")

Example output for ion gyro orbit tracing:
![Alt text](gyroHF.png "Example ion gyro orbit output")

Example output for limited discharges:
![Alt text](limiter.gif "Example output of EQ, HF, T, video")

HEAT Dash / plotly GUI:
![Alt text](gui1.png "HEAT DASH GUI")
![Alt text](gui2.png "HEAT DASH GUI")
