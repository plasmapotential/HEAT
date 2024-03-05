Input File Reference
####################
This page describes the variables used in the HEAT input file, usually called
X_input.csv where X is the machine name.  

It also describes the format of the PFC input file, which is used for 
defining which PFCs to perform calculations on.

Lastly, it describes the radiated power input file, which is used as an emission
source for photon tracing.


X_input.csv File Description
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Below is a description of the input file variables that can be found in the X_input.csv file,
where X is the machine name.  

.. autofunction:: ioClass.IO_HEAT.allowed_class_vars
.. autofunction:: CADClass.CAD.allowed_class_vars
.. autofunction:: MHDClass.MHD.allowed_class_vars
.. autofunction:: heatfluxClass.heatFlux.allowed_class_vars   
.. autofunction:: gyroClass.GYRO.allowed_class_vars    
.. autofunction:: radClass.RAD.allowed_class_vars    
.. autofunction:: openFOAMclass.OpenFOAM.allowed_class_vars  
.. autofunction:: plasma3DClass.plasma3D.allowed_class_vars  


References:
-----------

#. T Eich et al, 2013 Nuclear Fusion, vol53 no9
#. D Brunner et al, 2013 Nuclear Fusion, vol58 no9
#. J Horacek et al, 2016 Plasma Phys. Control. Fusion 58
#. M Makowski et al, 2012  Physics of Plasmas 19, 056122
#. M Greenwald et al, 2002, Plasma Phys. Control. Fusion, vol44 no8
#. T Looby et al 2022 Nuclear Fusion 62 106020


PFC File Description
^^^^^^^^^^^^^^^^^^^^
The HEAT PFC file is a csv that describes which CAD objects we want to include in the
calculation, the resolution of the meshes, and some other variables related to magnetic
shadowing.  This file contains a separate line for each object in the CAD file that 
we want to use in the calculation.

.. autofunction:: engineClass.engineObj.readPFCfile


Photon Emission File Description
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The photon emission file is used for photon radiation flux calculations.  The contents of this
file are described in a comment of the read function:

.. autofunction:: radClass.RAD.read2DSourceFile


Batch File Description
^^^^^^^^^^^^^^^^^^^^^^
The batchFile is used to run HEAT in batch mode (terminal mode).  Running HEAT in batch 
mode is usually much faster, and can be used to perform a large number of HEAT runs
serially with a single command.  The following comment describes the batchFile used by
HEAT.  For more information about running in batch mode, see the TUI Tutorial section.

.. autofunction:: terminalUI.TUI.simulationSchedule