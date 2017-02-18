function fitID = getFitID(conn)
	fit = exec(conn, ['SELECT MAX(`fitID`) FROM `Fits`']);
	fit = fetch(fit);
	if isnan(fit.Data{1}) | strcmp(fit.Data{1}, 'No data')
		fitID = 0;
	else
		fitID = fit.Data{1}+1;
	end
end