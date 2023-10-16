Point(1) = {0,  0.35, 0, 0.005};
Point(2) = {-Sqrt(3)*0.35/2, 0.35/2, 0, 0.005};
Point(3) = {-Sqrt(3)*0.35/2, -0.35/2, 0, 0.005};
Point(4) = {0, -0.35, 0, 0.005};
Point(5) = {Sqrt(3)*0.35/2, -0.35/2, 0, 0.005};

Point(6) = {Sqrt(3)*0.35/2,  -0.35/12, 0, 0.005};
Point(7) = {Sqrt(3)*0.35/2 - 0.35/6, -0.35/12, 0, 0.005};
Point(8) = {Sqrt(3)*0.35/2 - 0.35/6, 0.35/12, 0, 0.005};
Point(9) = {Sqrt(3)*0.35/2,  0.35/12, 0, 0.005};

Point(10) = {Sqrt(3)*0.35/2, 0.35/2, 0, 0.005};



Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 1};
Line(11) = {6, 9};

Curve Loop(12) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
Curve Loop(13) = {-7, -6, 11, -8};


Plane Surface(1) = {12};
Plane Surface(2) = {13};

//Physical Curve("left", 7) = {1, 2, 3};
Physical Curve("left", 7) = {2};
Physical Curve("load", 8) = {5, 9, 11};
Physical Surface("beam", 3) = {1, 2};
Physical Surface("omega", 4) = {2};

