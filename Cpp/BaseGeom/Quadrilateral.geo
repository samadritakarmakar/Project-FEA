// Gmsh project created on Mon Mar 18 23:34:24 2019
SetFactory("OpenCASCADE");
//+
Point(1) = {-1, -1, 0, 1.0};
//+
Point(2) = {1, -1, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {-1, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {3, 4, 1, 2};
//+
Plane Surface(1) = {2};
//+
Transfinite Surface {1};
//+
Transfinite Surface {1};
//+
Transfinite Surface {1} Right;
//+
Transfinite Surface {1};
//+
Transfinite Surface {1};
//+
Characteristic Length {1} = 1;
//+
Field[1] = Distance;
//+
Delete Field [1];
//+
Transfinite Surface {1};
//+
Transfinite Surface {1};
