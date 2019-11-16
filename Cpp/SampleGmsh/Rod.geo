 
//+
SetFactory("OpenCASCADE");
//+
Circle(1) = {0, 0, 0, 0.05, 0, 2*Pi};

//+
Curve Loop(1) = {1};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 0.1} {
  Surface{1}; 
}
//+
Physical Surface("air_side") = {3, 2};
//+
Physical Surface("hot_side") = {1};

//+
Physical Volume("whole_vol") = {1};
