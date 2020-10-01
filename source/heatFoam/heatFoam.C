/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    heatFoam

Group
    grpBasicSolvers

Description
    Engineer:  Tom Looby
    Date: Summer / Fall 2020 (covid days...)

    Laplace equation solver for a scalar quantity with temperature dependent
    diffusion constant.

    Usually this is being used for heat flux calculations in tokamaks.
    Built for Heat flux Engineering Analysis Toolkit (HEAT), a code that
    calculates heat flux and temperature in tokamak plasma facing components (PFCs).

    HEAT outputs a heat flux calculated on PFC surface, then calls this package
    to solve for temperature.

    heat flux (HF) is assigned to FVM boundary (STLpatch) using
    timeVaryingMappedFixedValue and a groovy boundary coundition (groovyBC)
    that finds boundary gradT by solving Fourier's law using temperature
    dependent thermal conductivity.  All of that is handled in the OpenFOAM
    case directory.

    for more info on HEAT see github: https://github.com/plasmapotential/HEAT

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( DT \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        DT   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

//Need these for DT(T), thermCond(T)
#include "IFstream.H"
#include "graph.H"
#include "interpolateXY.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    //Need these for DT(T), thermCond(T), Cp(T)
    #include "initContinuityErrs.H"
    #include "interpolateProperties.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    label patchID = mesh.boundaryMesh().findPatchID("STLpatch");

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info<< "Reading Heat Flux HF\n" << endl;
        volScalarField HF
        (
            IOobject
            (
                "HF",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
        qFVM = HF;
        Info << max(HF) << nl << endl;
        //interpolate DT and thermCond
        #include "interpolateProperties.H"

        while (simple.correctNonOrthogonal())
        {

            //minmax handled in controlDict functionObject (fieldMinMax.H)
            //Info << max(HF) << nl << endl;
            //Info << max(T) << nl << endl;
            //Info << min(DT) << nl << endl;
            //Info << min(thermCond) << nl << endl;

            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        #include "write.H"

        qINV=thermCond*mag(fvc::grad(T))*1000;

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
