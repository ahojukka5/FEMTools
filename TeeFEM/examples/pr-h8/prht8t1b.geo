// Pintarakenteet, harjoitus 8, tehtävä 1 a)

a = 1;
b = 2;
lc = 1;

n = 20;
nx = n;
ny = b/a*nx;

Point(1) = {0,0,0,lc};
Extrude {a, 0, 0} {
  Point{1};
  Layers{nx};
}
Extrude {0, b, 0} {
  Line{1};
  Layers{ny};
}
Physical Line("GA1") = {3};
Physical Line("GA2") = {1};
Physical Line("GA3") = {4};
Physical Line("GA4") = {2};
Physical Surface("OM1") = {5};
