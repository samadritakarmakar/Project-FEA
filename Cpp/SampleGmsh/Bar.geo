 
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {200, 0, 0, 1.0};
//+
Point(3) = {200, 2000, 0, 1.0};
//+
Point(4) = {2, 2000, 0, 1.0};
//+
Point(5) = {0, 2000, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 1};
//+
Physical Curve("ForceLine", 1) = {4};
//+
Curve Loop(1) = {3, 4, 5, 1, 2};
//+
Plane Surface(1) = {1};
//+
Extrude {0, 0, 60} {
  Surface{1}; 
}
//+
Characteristic Length {7} = 0.1;
//+
Physical Surface("LockedSurf", 2) = {27};
//+
Physical Surface("forceSurf") = {15};
//+
Transfinite Curve {7, 10} = 5 Using Progression 1;
//+
Transfinite Curve {9, 11} = 25 Using Progression 1;
//+
Transfinite Curve {13, 18, 22, 26} = 2 Using Progression 1;
//+
Transfinite Curve {1, 3} = 5 Using Progression 1;
//+
Transfinite Curve {13, 14, 26, 22} = 3 Using Progression 1;
//+
Transfinite Curve {11, 2, 5, 9} = 25 Using Progression 1;
