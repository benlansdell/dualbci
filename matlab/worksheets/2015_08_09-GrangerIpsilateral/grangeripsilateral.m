%Compute sum of granger scores for out degree for given unit and nev file
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%Connections to ipsi units
data = exec(conn, ['SELECT f.unit, f.`nev file`, ge.score, ge.fromunit FROM `fits` f '...
'INNER JOIN estimates_granger ge '...
'ON f.id = ge.id '...
'WHERE f.modelID = 2 AND f.unit > 96']);
data = fetch(data);
ipsidata = data.Data;
ipsidata = cell2mat(ipsidata(:,3));

%Connections to contra units
data = exec(conn, ['SELECT f.unit, f.`nev file`, ge.score, ge.fromunit FROM `fits` f '...
'INNER JOIN estimates_granger ge '...
'ON f.id = ge.id '...
'WHERE f.modelID = 2 AND f.unit < 96']);
data = fetch(data);
contradata = data.Data;
contradata = cell2mat(contradata(:,3));

%Histograms of each
figure
hist((contradata), 400)
xlabel('Contra granger in')
%xlim([0 150])
figure 
hist((ipsidata), 400)
xlabel('Ipsi granger in')
%xlim([0 150])

%Histogram  
figure
b1 = bar(hist(log(ipsidata), 400) ./ sum(hist(log(ipsidata), 400)))
set(get(gca,'child'),'FaceColor',[0 0.5 0.5],'EdgeColor',[0 .5 .5]);
hold on
b2 = bar(hist(log(contradata), 400) ./ sum(hist(log(contradata), 400)))
set(get(b2,'child'),'FaceColor',[0 0.5 0.5],'EdgeColor',[0 0 .7]);
legend('Ipsi-', 'Contra-')
xlabel('log_{10}(Granger score)')
ylabel('density')
%xlim([0 150])
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/MVGCanalysis1.eps')
