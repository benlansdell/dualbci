function coeffID = getCoeffID(conn)
	coeff = exec(conn, ['SELECT MAX(`fitID`) FROM `Fits`']);
	coeff = fetch(coeff);
	if isnan(coeff.Data{1}) | strcmp(coeff.Data{1}, 'No data')
		coeffID = 0;
	else
		coeffID = coeff.Data{1}+1;
	end
end