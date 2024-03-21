Installation
############

To run HEAT you will need to install the HEAT docker container from dockerhub.

The link below explains how to install the HEAT docker container and get HEAT up and running.

.. toctree::
   :maxdepth: 0

   docker

To download test HEAT case (deprecated)
***************************************
For HEAT versions equal or greater to v4.0, there are test cases inside the container.
These tests cases serve two purposes:
 * They enable integration testing for HEAT developers in CI
 * They enable HEAT users to test their local setup

These test cases can be found in the <HEATsource>/tests directory, where <HEATsource> is the 
location of your source code.  Inside the container, the <HEATsource> directory can be found at: 
/root/source/HEAT

Alternatively, you can test your HEAT installation by using a test case provided on the releases tab on github.  
The test case can be downloaded and extracted by using the following commands
(again replacing <tag> with latest version ie vX.X.X)::

    wget https://github.com/plasmapotential/HEAT/releases/download/<tag>/testCase_<tag>.tar.gz
    tar -xvzf testCase_<tag>.tar.gz


Test Case Installation Video (from release tab on github):

    .. raw:: html

        <div style="position: relative; padding-bottom: 2%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
            <iframe width="560" height="315" src="https://www.youtube.com/embed/PiU2yCTtfq0" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
        </div>
