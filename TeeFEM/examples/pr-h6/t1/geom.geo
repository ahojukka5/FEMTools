// Pintarakenteet, harjoitustehtävät 6, geometria
D = 0.6; // Pilarin halkaisija
L = 20;
H = 10; 
hx = 4;
hy = 2;
lc = 0.2;

// Ulkomitat
Point(1) = {0,0,0,lc};
Point(2) = {L,0,0,lc};
Point(3) = {L,H,0,lc};
Point(4) = {0,H,0,lc};

// Pilari1
Point(5) = {hx,hy,0,lc};
Point(6) = {hx-D/2,hy,0,lc};
Point(7) = {hx+D/2,hy,0,lc};
// Pilari2
Point(8) = {L-hx,hy,0,lc};
Point(9) = {L-hx-D/2,hy,0,lc};
Point(10) = {L-hx+D/2,hy,0,lc};
// Pilari3
Point(11) = {hx,H-hy,0,lc};
Point(12) = {hx-D/2,H-hy,0,lc};
Point(13) = {hx+D/2,H-hy,0,lc};
// Pilari4
Point(14) = {L-hx,H-hy,0,lc};
Point(15) = {L-hx-D/2,H-hy,0,lc};
Point(16) = {L-hx+D/2,H-hy,0,lc};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {7, 5, 6};
Circle(6) = {6, 5, 7};
Circle(7) = {13, 11, 12};
Circle(8) = {12, 11, 13};
Circle(9) = {10, 8, 9};
Circle(10) = {9, 8, 10};
Circle(11) = {16, 14, 15};
Circle(12) = {15, 14, 16};
Line Loop(13) = {8, 7};
Plane Surface(14) = {13};
Line Loop(15) = {12, 11};
Plane Surface(16) = {15};
Line Loop(17) = {9, 10};
Plane Surface(18) = {17};
Line Loop(19) = {5, 6};
Plane Surface(20) = {19};
Line Loop(21) = {4, 1, 2, 3};
Plane Surface(22) = {13, 15, 17, 19, 21};

// Pinnat ja reunat
Physical Surface("OM1") = {14};
Physical Surface("OM2") = {20};
Physical Surface("OM3") = {16};
Physical Surface("OM4") = {18};
Physical Surface("OM5") = {22};
Physical Line("GA1") = {7, 8};
Physical Line("GA2") = {5, 6};
Physical Line("GA3") = {11, 12};
Physical Line("GA4") = {9, 10};
Physical Line("GA5") = {3};
Physical Line("GA6") = {2};
Physical Line("GA7") = {1};
Physical Line("GA8") = {4};
