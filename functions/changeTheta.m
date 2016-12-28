function [deltaH, theta, torquesize] = changeTheta(torqueFE, torqueRU, BCItheta)
	%Compute change in angle
	torqueH = atan(torqueRU./torqueFE); %theta
	torquesize = sqrt(torqueRU.^2 + torqueFE.^2); %speed
	piplus = torqueFE < 0 & torqueRU > 0;
	piminus = torqueFE < 0 & torqueRU < 0;
	%Between plus/minus pi
	torqueH(piplus) = torqueH(piplus) + pi;
	torqueH(piminus) = torqueH(piminus) - pi;
	%Between 0 and 2pi
	torqueH(torqueH<0) = torqueH(torqueH<0)+2*pi;
	theta = torqueH;
	%Compute difference, take modulus to be between 0 and 2pi
	deltaH = mod(torqueH - BCItheta, 2*pi);
	%Change back to plus/minus 2pi, take abs value
	deltaH(deltaH>pi) = abs(deltaH(deltaH>pi)-2*pi);
end