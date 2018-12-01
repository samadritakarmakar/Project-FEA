// Gmsh project created on Mon Oct 29 20:54:31 2018
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 4, 0, 1.0};
//+
Point(4) = {4, 4, 0, 1.0};
//+
Point(5) = {4, 10, 0, 1.0};
//+
Point(6) = {0, 10, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {3, 4};
//+
Line(3) = {5, 4};
//+
Line(4) = {5, 6};
//+
Line(5) = {1, 6};
//+
Line(6) = {3, 2};
