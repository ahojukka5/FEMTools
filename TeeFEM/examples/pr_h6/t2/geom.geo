// Pintarakenteet, harjoitustehtävät 6, tehtävä 2, geometria
L = 9;
lc = 1;
n = 50;

Point(1) = {0,0,0,lc};
Extrude {L, 0, 0} {
  Point{1};
  Layers{n};
}
Extrude {0, L, 0} {
  Line{1};
  Layers{n};
}
Physical Line("GA1") = {2};
Physical Line("GA2") = {4};
Physical Line("GA3") = {1};
Physical Line("GA4") = {3};
Physical Surface("OM1") = {5};
