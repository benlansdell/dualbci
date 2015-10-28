SELECT mc1.causaldensity mc1cd, bc.causaldensity bccd, mc2.causaldensity mc2cd, f1.`nev file` mc1nev
FROM
FitsMVGC mc1 
INNER JOIN Fits f1
ON f1.id = mc1.id 
INNER JOIN AnalysisLinear al 
ON al.`manualrecording` = f1.`nev file`
INNER JOIN Fits f2
ON f2.`nev file` = al.`1DBCrecording` AND f2.modelID = 13
INNER JOIN FitsMVGC bc 
ON bc.id = f2.id 
INNER JOIN Fits f3
ON f3.`nev file` = al.`manualrecordingafter` AND f3.modelID = 13
INNER JOIN FitsMVGC mc2
ON mc2.id = f3.id
WHERE f1.modelID = 13






SELECT mc1.score mc1score, bc.score bcscore, mc2.score mc2score, f1.`nev file` mc1nev, mc1.fromunit
FROM
GrangerEstimates mc1 
INNER JOIN Fits f1
ON f1.id = mc1.id 
INNER JOIN AnalysisLinear al 
ON al.`manualrecording` = f1.`nev file`
INNER JOIN Fits f2
ON f2.`nev file` = al.`1DBCrecording` AND f2.modelID = 13
INNER JOIN GrangerEstimates bc 
ON bc.id = f2.id AND mc1.fromunit = bc.fromunit 
INNER JOIN Fits f3
ON f3.`nev file` = al.`manualrecordingafter` AND f3.modelID = 13
INNER JOIN GrangerEstimates mc2
ON mc2.id = f3.id AND mc2.fromunit = mc1.fromunit
WHERE f1.modelID = 13




SELECT mc1.score mc1score, bc.score bcscore, mc2.score mc2score, f1.`nev file` mc1nev, mc1.fromunit
FROM
GrangerEstimates mc1 
INNER JOIN Fits f1
ON f1.id = mc1.id 
INNER JOIN AnalysisLinear al 

ON al.`manualrecording` = f1.`nev file`
INNER JOIN Fits f2
ON f2.`nev file` = al.`1DBCrecording` AND f2.modelID = 13
INNER JOIN GrangerEstimates bc 
ON bc.id = f2.id AND mc1.fromunit = bc.fromunit 
INNER JOIN Fits f3
ON f3.`nev file` = al.`manualrecordingafter` AND f3.modelID = 13
INNER JOIN GrangerEstimates mc2
ON mc2.id = f3.id AND mc2.fromunit = mc1.fromunit
INNER JOIN BCIUnits bci
ON bci.ID = al.`1DBCrecording` AND bci.unit = mc1.fromunit
WHERE f1.modelID = 13



SELECT mc1.score mc1score, bc.score bcscore, mc2.score mc2score, bci.unit IS NOT NULL bciunit, f1.`nev file` mc1nev, mc1.fromunit
FROM
GrangerEstimates mc1 
INNER JOIN Fits f1
ON f1.id = mc1.id 
INNER JOIN AnalysisLinear al 
ON al.`manualrecording` = f1.`nev file`
INNER JOIN Fits f2
ON f2.`nev file` = al.`1DBCrecording` AND f2.modelID = 13
INNER JOIN GrangerEstimates bc 
ON bc.id = f2.id AND mc1.fromunit = bc.fromunit 
INNER JOIN Fits f3
ON f3.`nev file` = al.`manualrecordingafter` AND f3.modelID = 13
INNER JOIN GrangerEstimates mc2
ON mc2.id = f3.id AND mc2.fromunit = mc1.fromunit
LEFT JOIN BCIUnits bci
ON bci.ID = al.`1DBCrecording` AND bci.unit = mc1.fromunit
WHERE f1.modelID = 13