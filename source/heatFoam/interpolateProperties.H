Info << "Reading DT (thermal diffusivity) property file and interpolate" << endl;
IFstream file_DT(runTime.path()/runTime.constant()/"DT");
//Read file
graph DTGraph
(
    "DT_data_file",
    "T_data",
    "DT_data",
    file_DT
);

//Interpolation internal field
DT.field() = interpolateXY
(
    T.field(),
    DTGraph.x(),
    DTGraph.y()
);

//Interpolation boundary field
forAll(DT.boundaryFieldRef(), patchi)
{
    DT.boundaryFieldRef()[patchi] = interpolateXY
    (
      T.boundaryFieldRef()[patchi],
      DTGraph.x(),
      DTGraph.y()
    );
}


Info << "Reading thermCond (thermal conductivity) property file and interpolate" << endl;
IFstream file_thermCond(runTime.path()/runTime.constant()/"thermCond");
//Read file
graph thermCondGraph
(
    "thermCond_data_file",
    "T_data",
    "thermCond_data",
    file_thermCond
);

//Interpolation internal field
thermCond.field() = interpolateXY
(
    T.field(),
    thermCondGraph.x(),
    thermCondGraph.y()
);

//Interpolation boundary field
forAll(thermCond.boundaryFieldRef(), patchi)
{
    thermCond.boundaryFieldRef()[patchi] = interpolateXY
    (
      T.boundaryFieldRef()[patchi],
      thermCondGraph.x(),
      thermCondGraph.y()
    );
}


//We don't need Cp here, but I included it for future reference
//You would need to add a Cp file to the /constants directory in the OF case
//
//Info << "Reading Cp (specific heat) property file and interpolate" << endl;
//IFstream file_Cp(runTime.path()/runTime.constant()/"Cp");
////Read file
//graph CpGraph
//(
//    "Cp_data_file",
//    "T_data",
//    "Cp_data",
//    file_Cp
//);
//
////Interpolation internal field
//Cp.field() = interpolateXY
//(
//    T.field(),
//    CpGraph.x(),
//    CpGraph.y()
//);
//
////Interpolation boundary field
//forAll(Cp.boundaryFieldRef(), patchi)
//{
//    Cp.boundaryFieldRef()[patchi] = interpolateXY
//    (
//      T.boundaryFieldRef()[patchi],
//      CpGraph.x(),
//      CpGraph.y()
//    );
//}
