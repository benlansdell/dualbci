function trans_angles2 = translate_angles(angles1, angles2)
	%Shift the second set of angles 
	trans_angles2 = angles2;
	abovepi1 = angles1 > pi;
	belowangles1 = angles2 < angles1-pi;
	trans_angles2(abovepi1&belowangles1) = angles2(abovepi1&belowangles1)+2*pi;

	belowpi1 = angles1 < pi;
	belowangles1 = angles2 > angles1+pi;
	trans_angles2(belowpi1&belowangles1) = angles2(belowpi1&belowangles1)-2*pi;
end