// add start and end point (x,y,z,h-factor)
// note: only the x-coordinates matter
Point(1) = {-3, 0, 0};
Point(2) = {6, 0, 0};

// connect the two points with a line
Line(1) = {1, 2};

// assign physical contour 'outer_bc' to the two points
Physical Point("outer_bc") = {1, 2};

// assign physical region 'the_only_region' to the line segment
Physical Curve("the_only_region") = {1};

Mesh.CharacteristicLengthMin = 3;
