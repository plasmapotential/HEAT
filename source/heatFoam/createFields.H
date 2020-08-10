Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Cp\n" << endl;

volScalarField Cp
(
    IOobject
    (
        "Cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("Cp",dimensionSet (0,2,-2,-1,0,0,0), 1)
);

Info<< "Reading diffusivity DT\n" << endl;

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
//        IOobject::READ_IF_PRESENT,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("DT",dimensionSet (0,2,-1,0,0,0,0), 1)
//    dimensionedScalar(dimViscosity, Zero)

);


Info<< "Reading conductivity thermCond\n" << endl;

volScalarField thermCond
(
    IOobject
    (
        "thermCond",
        runTime.timeName(),
        mesh,
//        IOobject::READ_IF_PRESENT,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
//    dimensionedScalar ("thermCond",dimensionSet (1,0,-3,1,0,0,0), 1)
    dimensionedScalar ("thermCond",dimensionSet (1,1,-3,-1,0,0,0), 1)
//    dimensionedScalar(dimViscosity, Zero)

);


volScalarField qINV
(
    IOobject
    (
        "qINV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("qINV",dimensionSet (1,0,-3,0,0,0,0), 1)
);

volScalarField qFVM
(
    IOobject
    (
        "qFVM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar ("qFVM",dimensionSet (1,0,-3,0,0,0,0), 1)
);


//Info<< "Reading Heat Flux HF\n" << endl;
//volScalarField HF
//(
//    IOobject
//    (
//        "HF",
//        runTime.timeName(),
//        mesh,
//        IOobject::MUST_READ,
//        IOobject::AUTO_WRITE
//    ),
//    mesh
//);


#include "createFvOptions.H"
