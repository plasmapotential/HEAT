    if (runTime.writeTime())
    {
        volVectorField gradT(fvc::grad(T));

        volScalarField gradTx
        (
            IOobject
            (
                "gradTx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradT.component(vector::X)
        );

        volScalarField gradTy
        (
            IOobject
            (
                "gradTy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradT.component(vector::Y)
        );

        volScalarField gradTz
        (
            IOobject
            (
                "gradTz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            gradT.component(vector::Z)
        );

//        volVectorField thermCondgradT
//        (
//            IOobject
//            (
//                "flux",
//                runTime.timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::AUTO_WRITE
//            ),
//            //multiply by 1000 to handle the fact that HEAT does everything in mm
//            thermCond*gradT*1000
//        );

        volScalarField thermCondgradT
        (
            IOobject
            (
                "fluxMag",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            //multiply by 1000 to handle the fact that HEAT does everything in mm
            thermCond*mag(gradT)*1000
        );




//        dimensionedScalar maxT
//        (
//            IOobject
//            (
//                "maxT",
//                runTime.timeName(),
//                mesh,
//                IOobject::NO_READ,
//                IOobject::AUTO_WRITE
//            ),
//            max(T).value()
//        );


        runTime.write();
    }
