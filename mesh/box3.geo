SetFactory("OpenCASCADE");

//Outer box
Point(1) = {0, 0, 0, 0.01};
Point(2) = {0.5, 0, 0, 0.01};
Point(3) = {0.5, 1/9, 0, 0.01};
Point(4) = {0, 1/9, 0, 0.01};
Point(5) = {0, 0, 1/9, 0.01};
Point(6) = {0.5, 0, 1/9, 0.01};
Point(7) = {0.5, 1/9, 1/9, 0.01};
Point(8) = {0, 1/9, 1/9, 0.01};

//Inner box
Point(111) = {0.5-1/45, 2/45, 2/45, 0.01};
Point(112) = {0.5, 2/45, 2/45, 0.01};
Point(113) = {0.5, 3/45, 2/45, 0.01};
Point(114) = {0.5-1/45, 3/45, 2/45, 0.01};
Point(115) = {0.5-1/45, 2/45, 3/45, 0.01};
Point(116) = {0.5, 2/45, 3/45, 0.01};
Point(117) = {0.5, 3/45, 3/45, 0.01};
Point(118) = {0.5-1/45, 3/45, 3/45, 0.01};

//Outer box
Line(1) = {1, 2};
Line(2) = {2, 6};
Line(3) = {6, 5};
Line(4) = {5, 1};
Line(5) = {5, 8};
Line(6) = {8, 4};
Line(7) = {4, 1};
Line(8) = {8, 7};
Line(9) = {7, 3};
Line(10) = {3, 4};
Line(11) = {6, 7};
Line(12) = {3, 2};

//Inner box
Line(111) = {111, 112};
Line(112) = {112, 116};
Line(113) = {116, 115};
Line(114) = {115, 111};
Line(115) = {115, 118};
Line(116) = {118, 114};
Line(117) = {114, 111};
Line(118) = {118, 117};
Line(119) = {117, 113};
Line(1110) = {113, 114};
Line(1111) = {116, 117};
Line(1112) = {113, 112};

//Outer box
Line Loop(1) = {4, -7, -6, -5}; //left wall
Plane Surface(1) = {1};

Line Loop(2) = {2, 11, 9, 12}; // right wall
Plane Surface(2) = {2};

Line Loop(3) = {3, 5, 8, -11}; // front wall
Plane Surface(3) = {3};

Line Loop(4) = {10, 7, 1, -12}; // back wall
Plane Surface(4) = {4};

Line Loop(5) = {9, 10, -6, 8}; // top wall
Plane Surface(5) = {5};

Line Loop(6) = {3, 4, 1, 2}; // bottom wall
Plane Surface(6) = {6};


//Inner box
Line Loop(111) = {114, -117, -116, -115}; //left wall
Plane Surface(111) = {111};

Line Loop(112) = {112, 1111, 119, 1112}; // right wall
Plane Surface(112) = {112};

Line Loop(113) = {113, 115, 118, -1111}; // front wall
Plane Surface(113) = {113};

Line Loop(114) = {1110, 117, 111, -1112}; // back wall
Plane Surface(114) = {114};

Line Loop(115) = {119, 1110, -116, 118}; // top wall
Plane Surface(115) = {115};

Line Loop(116) = {113, 114, 111, 112}; // bottom wall
Plane Surface(116) = {116};


//Outer box
Surface Loop(1) = {4, 5, 2, 6, 3, 1};
Volume(1) = {1};

//Inner box
Surface Loop(111) = {114, 115, 112, 116, 113, 111};
Volume(111) = {111};
Volume(3) = BooleanFragments{ Volume{1}; Delete; }{ Volume{111}; Delete; };
Coherence;

Physical Surface("leftwall", 7) = {1};
Physical Surface("rightwall", 8) = {2};
Physical Volume("beam", 3) = {2};
Physical Volume("omega", 4) = {111};

//SetFactory("OpenCASCADE");
//Box(14) = {0.5-1/45, 2/45, 2/45, 1/45, 1/45, 1/45};
