function [x, y, array] = map_to_pos2(electrode)
	%Conver to number, if not already
	if isstr(electrode)
		electrode = str2num(electrode);
	end
	electrode = floor(electrode);
	%Determine array
	if electrode > 96
		%Ipsilateral
		array = 2;
		electrode = electrode - 96;
	else
		%Contralateral
		array = 1;
	end
	%Determine position
	if electrode <= 8
		x = electrode + 1;
		y = 1;
	elseif electrode > 8 & electrode <= 88
		x = mod(electrode+1, 10)+1;
		y = floor((electrode+1)/10)+1;
		if mod(y,2) == 0;
			x = 11-x;
		end
	else
		x = mod(electrode+2, 10)+1;
		y = 10;
	end
end