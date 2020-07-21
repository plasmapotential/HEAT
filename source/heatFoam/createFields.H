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

#include "createFvOptions.H"
