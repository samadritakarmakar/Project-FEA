 
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {200, 0, 0, 1.0};
//+
Point(3) = {200, 2000, 0, 1.0};
//+
Point(4) = {0, 2000, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Transfinite Curve {3, 1} = 5 Using Progression 1;
//+
Transfinite Curve {4, 2} = 25 Using Progression 1;
//+
Compound Curve {4, 3, 2, 1};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Recombine Surface {1};
//+
Extrude {0, 0, 60} {
  Surface{1}; Layers{2}; Recombine;
}
//+
Transfinite Volume{1} = {1, 2, 10, 6, 5, 14, 3, 4};
//+
Recombine Surface {100017, 100013, 100026, 1, 100025, 100021};
//+
Physical Curve("forceline") = {100009};
//+
Physical Surface("lockedsurface") = {17};
//+
Physical Surface("forcesurface") = {25};
