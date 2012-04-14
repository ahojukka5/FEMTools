// Sauva
lc = 1;
r = 1;
L = 10;
Point(1) = {0, 0, 0, lc};
Extrude {r, 0, 0} {
  Point{1};
}
Extrude {0, L, 0} {
  Line{1};
}
Extrude {{0, 1, 0}, {0, 0, 0}, 2*Pi/3} {
  Surface{5};
}
Extrude {{0, 1, 0}, {0, 0, 0}, 2*Pi/3} {
  Surface{22};
}
Extrude {{0, 1, 0}, {0, 0, 0}, 2*Pi/3} {
  Surface{39};
}
Physical Volume("OM1") = {1, 2, 3};
