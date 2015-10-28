
%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
results = exec(conn, ['SELECT * FROM fits_glm_vs_lm_decode']);
results = fetch(results);
results = results.Data;
results = cell2mat(results(:,2:3));
plot(results(:,1), results(:,2), '.', [0 0], [1 1], [-0.5 1], [-0.5 1])
xlabel('Correlation GLM')
ylabel('Correlation LM')
saveplot(gcf, './worksheets/2015_07_28-LinearVsGLMDecoding/results.eps')
