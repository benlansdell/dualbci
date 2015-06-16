function theta = taskTheta(dirstr)
	switch dirstr
	case 'east'
		theta = 0;
	case 'west'
		theta = pi;
	case 'north'
		theta = pi/2;
	case 'south'
		theta = 3*pi/2;
	end
end