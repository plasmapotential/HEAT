Info << "Reading DT property file and interpolate" << endl;
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
