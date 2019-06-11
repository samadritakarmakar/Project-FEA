//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 1};
//+
Curve Loop(1) = {1, 2, 3};
//+
Surface(1) = {1};
//+
Characteristic Length {1} = 2;
//+
Characteristic Length {3, 2, 1} = 3;
