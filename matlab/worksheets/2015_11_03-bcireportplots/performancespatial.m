%Compute sum of granger scores for out degree for given unit and nev file
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%data = exec(conn, ['SELECT f.unit, f.`nev file`, ge.score, ge.fromunit FROM `fits` f '...
%'INNER JOIN estimates_granger ge '...
%'ON f.id = ge.id '...
%'INNER JOIN experiment_tuning et '...
%'ON et.manualrecording = f.`nev file` '...
%'WHERE f.modelID = 2']);

data = exec(conn, ['SELECT ge.score, fl1.dir, fl2.dir, MOD(fl1.dir- fl2.dir, 2*PI()), flin1.`nev file`, ge.`fromunit`, flin1.`unit` FROM `fits` fgranger '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = fgranger.`nev file` AND flin1.`unit` = fgranger.`unit` '...
'INNER JOIN `estimates_granger` ge '...
'ON fgranger.id = ge.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = fgranger.`nev file` AND flin2.`unit` = ge.`fromunit` '...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
' INNER JOIN experiment_tuning et '...
' ON et.`manualrecording` = flin1.`nev file` '...
'WHERE flin1.modelID = 1 AND flin2.modelID = 1 AND fgranger.modelID = 2']);

data = fetch(data);
data = data.Data;

scores = cell2mat(data(:,1));
fromunits = cellfun(@str2num, data(:,6));
tounits = cellfun(@str2num, data(:,7));
angles = cell2mat(data(:,4));
tuningsim = cos(angles);
nG = length(fromunits);
fromposx = zeros(nG,1);
fromposy = zeros(nG,1);
fromarrays = zeros(nG,1);
toposx = zeros(nG,1);
toposy = zeros(nG,1);
toarrays = zeros(nG,1);
for idx = 1:nG
	[x,y,array] = map_to_pos(fromunits(idx));
	fromposx(idx) = x;
	fromposy(idx) = y;
	fromarrays(idx) = array;
	[x,y,array] = map_to_pos(tounits(idx));
	toposx(idx) = x;
	toposy(idx) = y;
	toarrays(idx) = array;
end
%Find the units on the same array:
contraidx = (toarrays == 1) & (fromarrays == 1);
fromposx = fromposx(contraidx);
fromposy = fromposy(contraidx);
toposx = toposx(contraidx);
toposy = toposy(contraidx);
scores = scores(contraidx);
tuningsim = tuningsim(contraidx);
distances = sqrt((fromposx-toposx).^2+(fromposy-toposy).^2);
angles = atan((fromposy-toposy)./(fromposx-toposx));



figure 
indices = randsample(length(distances), 2000);
scatter(distances(indices), scores(indices))

%figure
%scatter(angles(indices), scores(indices));

%lscores = log(scores);
smoothed = smooth(distances, scores, .1, 'rloess');
%figure
[xx,ind] = sort(distances);
hold on 
plot(xx,smoothed(ind),'r-'); %, xx,scores(ind),'b.')

%Bin into bins...
maxdist = 13;
nbins = 100;
bins = cell(nbins,1);
%Set scores to zero
for idx = 1:nbins
	bins{idx} = [];
end
for idx = 1:length(distances)
	%Figure out which bin to put score in...
	distance = distances(idx);
	tsim = tuningsim(idx);
	score = scores(idx);
	bin = ceil(nbins*(distance)/(maxdist))+1;
	bins{bin} = [bins{bin},score];

end
meanscores = [];
stdscores = [];
for idx = 1:nbins
	meanscores(idx) = mean(bins{idx});
	stdscores(idx) = std(bins{idx});
end
figure 
plot(meanscores, '.')

%2D heatmap
%Bin into bins...
maxdist = 13;
maxtuning = 2;
nbins = 10;
bins = cell(nbins,nbins);
%Set scores to zero
for idx = 1:nbins
	for j = 1:nbins
		bins{idx,j} = [];
	end
end
for idx = 1:length(distances)
	%Figure out which bin to put score in...
	distance = distances(idx);
	tsim = tuningsim(idx);
	score = scores(idx);
	bini = ceil(nbins*(distance)/(maxdist))+1;
	binj = ceil(nbins*(tsim+1)/(maxtuning));
	bins{bini, binj} = [bins{bini,binj},score];
end
meanscores = [];
stdscores = [];
for idx = 1:nbins
	for j = 1:nbins
		meanscores(idx,j) = mean(bins{idx,j});
		stdscores(idx,j) = std(bins{idx,j});
	end
end
figure 
imagesc(meanscores)
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/grangerspatialtuning.eps')


%Add noise to angles:
anglesn = angles + randn(size(angles))*.1;
%Bin into bins...
maxangle = pi/2;
nbins = 50;
bins = cell(nbins,1);
%Set scores to zero
for idx = 1:nbins
	bins{idx} = [];
end
for idx = 1:length(anglesn)
	%Figure out which bin to put score in...
	theta = anglesn(idx);
	score = scores(idx);
	bin = ceil(nbins*(theta+maxangle)/(2*maxangle));
	bin = min(max(bin, 1), nbins);
	bins{bin} = [bins{bin},score];
end
meanscores = [];
stdscores = [];
for idx = 1:nbins
	meanscores(idx) = mean(bins{idx});
	stdscores(idx) = std(bins{idx});
end
figure 
plot(meanscores, '.')
ylim([0 20])