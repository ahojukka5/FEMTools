// TeeFEM yksikkösuorakaide elementtien testaukseen. 
a = 1;
b = 1;
lc = 1;
n = 20;
na = a/b*n;
nb = n;

Point(1) = {0,0,0,lc};
Extrude {a, 0, 0} {
  Point{1};
  Layers{na};
}
Extrude {0, b, 0} {
  Line{1};
  Layers{nb};
}
Physical Surface("OM1") = {5};

Physical Line("GA1") = {3};
Physical Line("GA2") = {1};
Physical Line("GA3") = {4};
Physical Line("GA4") = {2};
