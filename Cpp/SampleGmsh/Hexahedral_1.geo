//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 10, 10, 0};
//+
Compound Curve {1, 2, 3, 4};
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Extrude {0, 0, 10} {
  Surface{1}; Layers{5}; Recombine;
}
//+
Physical Point("PointSet1") = {8, 4, 5, 1};
//+
Physical Point("PointSet2") = {3, 7, 6, 2};
//+
Physical Curve("CurveSet1") = {11, 12, 7, 9, 10, 5};
//+
Physical Curve("CurveSet2") = {3, 2, 1, 4, 8, 6};
//+
Physical Surface("SurfaceSet1") = {5, 6, 4};
//+
Physical Surface("SurfaceSet2") = {2, 3, 1};
//+
Physical Volume("VolumeSet1") = {1};
