//+
SetFactory("Built-in");
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
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Transfinite Curve {1, 2, 3, 4} = 1 Using Progression 1;

//+
Transfinite Surface {1} = {1, 2, 3, 4};
//+
Recombine Surface {1};
//+
Physical Point("pt1", 1) = {1};
//+
Physical Point("pt2", 2) = {2};
//+
Physical Point("pt3", 3) = {3};
//+
Physical Point("pt4", 4) = {4};
//+
Physical Curve("line1", 5) = {1};
//+
Physical Curve("line2", 6) = {2};
//+
Physical Curve("line3", 7) = {3};
//+
Physical Curve("line4", 8) = {4};
