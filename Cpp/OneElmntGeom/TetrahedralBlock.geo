//+
Point(1) = {-1, -1, -1, 1.0};
//+
Point(2) = {1, -1, -1, 1.0};
//+
Point(3) = {1, 1, -1, 1.0};
//+
Point(4) = {-1, 1, -1, 1.0};
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
Extrude {0, 0, 2} {
  Surface{1}; 
}
//+
Transfinite Curve {6, 11, 1, 12, 9, 8, 20, 3, 4, 16, 2, 7} = 1 Using Progression 1;
//+
Characteristic Length {5, 1, 14, 4, 3, 10, 6, 2} = 4;
//+
Transfinite Surface {26} = {14, 5, 6, 10};
//+
Transfinite Surface {1} = {1, 4, 3, 2};
//+
Transfinite Surface {13} = {1, 5, 6, 4};
//+
Transfinite Surface {17} = {6, 2, 3, 10};
//+
Transfinite Surface {13} = {6, 5, 2, 1};
//+
Transfinite Surface {21} = {3, 10, 14, 4};
//+
Transfinite Surface {25} = {4, 14, 5, 1};
//+
Physical Curve("locked") = {4};
//+
Physical Curve("lockedx") = {9};
//+
Physical Surface("forceSurf") = {17};
