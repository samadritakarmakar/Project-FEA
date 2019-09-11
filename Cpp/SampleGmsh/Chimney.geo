 
//+
Point(1) = {0, -35, 0, 1.0};
//+
Point(2) = {0, -0.3, 0, 1.0};
//+
Point(3) = {0, 0.3, 0, 1.0};
//+
Point(4) = {7, 0.3, 0, 1.0};
//+
Point(5) = {7, -0.3, 0, 1.0};

//+
Point(6) = {0, 35, 0, 1.0};
//+
Point(7) = {20, 35, 0, 1.0};
//+
Point(8) = {20, -35, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 5};
//+
Line(3) = {5, 4};
//+
Line(4) = {3, 4};
//+
Line(5) = {6, 3};
//+
Line(6) = {7, 6};
//+
Line(7) = {8, 7};
//+
Line(8) = {8, 1};
//+
Curve Loop(1) = {8, 1, 2, 3, -4, -5, -6, -7};
//+
Plane Surface(1) = {1};
//+
Physical Curve("FluxLine") = {3};
//+
Physical Curve("ZeroLine") = {8, 1, 2, 4, 5};
//+
Characteristic Length {5, 4} = 0.01;
