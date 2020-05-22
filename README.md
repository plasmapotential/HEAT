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


The source directory contains code for the project.  The architecture is as follows:

source
   CADclass.py		main CAD module  
   MHDClass.py		main MHD module  
   heatfluxClass.py	main HF module  
   openFOAMClass.py	main openFOAM FEM module  
   toolsClass.py	module for extra tools  
   gfiles.py		module for gfiles  
   GUIscripts/		directory for GUI plotting scripts  
   htmlGUI/		directory containing all HTML/JS/CSS for GUI  
      HEATgui.py	script that launches HEAT with Flask Bindings  
      GUIclass.py       main GUI module  
   inputs/		location of HEAT input files  
   openFoamTemplates	location of OF template case  
