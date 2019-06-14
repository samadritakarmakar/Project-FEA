//+
Point(1) = {0, 0, 0};
//+
Point(2) = {1, 0, 0};
//+
Point(3) = {0, 0, 1};
//+
Point(4) = {0, 1, 0};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 2};
//+
Line(3) = {1, 3};
//+
Line(4) = {2, 4};
//+
Line(5) = {4, 1};
//+
Line(6) = {3, 4};
//+
Line Loop(1) = {1, -2, -3};
//+
Plane Surface(1) = {1};
//+
Line Loop(2) = {1, 4, 5};
//+
Plane Surface(2) = {2};
//+
Line Loop(3) = {2, 4, -6};
//+
Plane Surface(3) = {3};
//+
Line Loop(4) = {5, 3, 6};
//+
Plane Surface(4) = {4};
//+
Surface Loop(1) = {1, 2, 3, 4};
//+
Volume(1) = {1};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
Transfinite Surface {2};
//+
Transfinite Surface {4};

//+
Characteristic Length {4, 3, 2, 1} = 3;
//+
Transfinite Curve {3, 5, 6, 4, 2, 1} = 1 Using Progression 1;
