%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
data = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`, fl.dir, bci.angle, MOD(fl.dir-bci.angle, 2*PI()) FROM GrangerEstimates ge INNER JOIN Fits f ON f.fitID = ge.fitID INNER JOIN Recordings rec ON rec.`nev file` = f.`nev file` INNER JOIN BCIUnits bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit INNER JOIN AnalysisLinear al ON al.`1DBCrecording` = f.`nev file` INNER JOIN Fits f1 ON f1.`nev file` = al.`manualrecording` AND f1.`modelID` = 1 and f1.unit = ge.fromunit INNER JOIN FitsLinear fl ON f1.`fitID` = fl.fitID WHERE f.modelID = 3']);
data = fetch(data);
data = data.Data;

hold on 
for idx = 1:(size(data,1)-1)
	if strcmp(data{idx, 4}, data{idx+1, 4})
		x = data{idx,7}
		y = data{idx+1,7};
		c = abs(data{idx,2}-data{idx+1,2});
		scatter(x,y,[],c);
	end
end