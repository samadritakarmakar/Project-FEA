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
Physical Curve("CurveSet1") = {3, 8, 10, 11};
//+
Physical Curve("CurveSet2") = {1, 5, 6, 7};
//+
Physical Surface("SurfaceSet1") = {4};
//+
Physical Surface("SurfaceSet2") = {2};
//+
Physical Volume("VolumeSet1") = {1};

