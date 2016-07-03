Point(1)={0,0,0};
Point(2)={0,1,0};
Point(3)={1,1,0};
Point(4)={1,0,0};

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

Line loop(5)={1,2,3,4};

Plane Surface(6)={5};

Periodic Line {1}={3};
Periodic Line {2}={4};

Physical Surface("interface") = {6};
Physical Line("BZboundary") = {1,2,3,4};
