function label = findLabel(idx, k)
	label = -1;
	for i = 1:size(k, 1)
		if any(k{i,2} == idx)
			label = k{i,1};
		end
	end
	assert(isstr(label), ['Cannot find label for ' num2str(idx) 'th parameter. Check data structure.'])
end