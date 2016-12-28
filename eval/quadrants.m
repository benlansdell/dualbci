function quads = quadrants(angles)
	quads = zeros(size(angles));
	fe1 = angles < (pi/4) | angles > (7*pi/4);
	ru2 = angles > (pi/4) & angles < (3*pi/4);
	fe3 = angles > (3*pi/4) & angles < (5*pi/4);
	ru4 = angles > (5*pi/4) & angles < (7*pi/4);
	quads(fe1) = 1;
	quads(ru2) = 2;
	quads(fe3) = 3;
	quads(ru4) = 4;
end