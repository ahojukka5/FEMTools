// Geometria

R = 0.2;
L = 1;
t = 0.100;
lc = 0.050;

Point(1) = {R,0,0,lc};
Point(2) = {R+t,0,0,lc};
Point(3) = {R+t,L/2,0,lc};
Point(4) = {R+t,L,0,lc};
Point(5) = {R,L,0,lc};
Point(6) = {R,L/2,0,lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line Loop(7) = {3, 4, 5, 6, 1, 2};
Plane Surface(8) = {7};
Physical Line("GA1") = {1, 4};
Physical Point("P1") = {3};
Physical Surface("OM1") = {8};
