cl__1 = 1;
Point(1) = {0, 0, 0, 1};
Point(2) = {-1, 0, 0, 1};
Point(3) = {1, 0, 0, 1};
Point(4) = {0, 1, 0, 1};
Point(5) = {0, -1, 0, 1};
Point(6) = {0, 0, 1, 1};
Point(7) = {0, 0, -1, 1};
Circle(1) = {3, 1, 4};
Circle(2) = {3, 1, 6};
Circle(3) = {2, 1, 6};
Circle(4) = {2, 1, 7};
Circle(5) = {3, 1, 7};
Circle(6) = {4, 1, 6};
Circle(7) = {2, 1, 4};
Circle(8) = {4, 1, 7};
Circle(9) = {6, 1, 5};
Circle(10) = {5, 1, 3};
Circle(11) = {5, 1, 2};
Circle(12) = {7, 1, 5};
Line Loop(14) = {7, 6, -3};
Ruled Surface(14) = {14};
Line Loop(16) = {6, -2, 1};
Ruled Surface(16) = {16};
Line Loop(18) = {1, 8, -5};
Ruled Surface(18) = {18};
Line Loop(20) = {7, 8, -4};
Ruled Surface(20) = {20};
Line Loop(22) = {9, 10, 2};
Ruled Surface(22) = {22};
Line Loop(24) = {3, 9, 11};
Ruled Surface(24) = {24};
Line Loop(27) = {10, 5, 12};
Ruled Surface(27) = {27};
Line Loop(29) = {11, 4, 12};
Ruled Surface(29) = {29};
Surface Loop(31) = {16, 14, 20, 18, 27, 22, 24, 29};
Volume(31) = {31};
