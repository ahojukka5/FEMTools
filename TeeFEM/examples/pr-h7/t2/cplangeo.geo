// TeeFEM yksikkösuorakaide elementtien testaukseen. 
a = 1;
b = 1;
lc = 1;
n = 50;
na = n;
nb = 2*b/a*n;

Point(1) = {0,0,0,lc};
Extrude {a/2, 0, 0} {
  Point{1};
  Layers{na};
}

Extrude {a/2, 0, 0} {
  Point{2};
  Layers{na};
}

Extrude {0, b, 0} {
 Line{1,2};
 Layers{nb};
}

Physical Line("GA1") = {4};
Physical Line("GA2") = {1, 2};
Physical Line("GA3") = {9};
Physical Line("GA4") = {7, 3};
Physical Line("GA5") = {5};
Physical Surface("OM1") = {6, 10};
