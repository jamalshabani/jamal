Point(1) = {0,  0.35, 0, 0.005};
Point(2) = {-Sqrt(3)*0.35/2, 0.35/2, 0, 0.005};
Point(3) = {-Sqrt(3)*0.35/2, 0, 0, 0.005};
Point(4) = {-Sqrt(3)*0.35/2, -0.35/2, 0, 0.005};
Point(5) = {0, -0.35, 0, 0.005};
Point(6) = {Sqrt(3)*0.35/2, -0.35/2, 0, 0.005};
Point(7) = {Sqrt(3)*0.35/2, 0, 0, 0.005};
Point(8) = {Sqrt(3)*0.35/2, 0.35/2, 0, 0.005};


Point(9) = {0,  0.35/10, 0, 0.005};
Point(10) = {-Sqrt(3)*0.35/20, 0.35/20, 0, 0.005};
Point(11) = {-Sqrt(3)*0.35/20, 0, 0, 0.005};
Point(12) = {-Sqrt(3)*0.35/20, -0.35/20, 0, 0.005};
Point(13) = {0, -0.35/10, 0, 0.005};
Point(14) = {Sqrt(3)*0.35/20, -0.35/20, 0, 0.005};
Point(15) = {Sqrt(3)*0.35/20, 0, 0, 0.005};
Point(16) = {Sqrt(3)*0.35/20, 0.35/20, 0, 0.005};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

Line(9) = {3, 11};
Line(10) = {11, 10};
Line(11) = {10, 9};
Line(12) = {9, 16};
Line(13) = {16, 15};
Line(14) = {15, 14};
Line(15) = {14, 13};
Line(16) = {13, 12};
Line(17) = {12, 11};
Line(18) = {15, 7};

//Upper hexa
Curve Loop(19) = {1, 2, 9, 10, 11, 12, 13, 18, 7, 8};

//Lower hexa
Curve Loop(20) = {3, 4, 5, 6, -18, 14, 15, 16, 17, -9};

//Omega_0
Curve Loop(21) = {-11, -10, -17, -16, -15, -14, -13, -12};


Plane Surface(1) = {19};
Plane Surface(2) = {20};
Plane Surface(3) = {21}; //Omega_0

Physical Curve("clamped", 7) = {2, 3, 5, 8};
Physical Curve("loadxx", 8) = {6, 7};
Physical Curve("loadxy", 9) = {1};
Physical Curve("loadyx", 10) = {4};
Physical Surface("beam", 3) = {1, 2};
Physical Surface("omega", 4) = {3};

