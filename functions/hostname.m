function name = hostname()
	[ret, name] = system('hostname');   
	if ret ~= 0,
	   if ispc
	      name = getenv('COMPUTERNAME');
	   else      
	      name = getenv('HOSTNAME');      
	   end
	end
	name = strtrim(lower(name));
	nB = min(50, length(name));
	name = name(1:nB);
end
