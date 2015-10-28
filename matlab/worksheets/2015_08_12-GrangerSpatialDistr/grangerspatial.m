%Compute sum of granger scores for out degree for given unit and nev file
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

data = exec(conn, ['SELECT f.unit, f.`nev file`, ge.score, ge.fromunit FROM `fits` f '...
'INNER JOIN estimates_granger ge '...
'ON f.id = ge.id '...
'INNER JOIN experiment_tuning et '...
'ON et.manualrecording = f.`nev file` '...
'WHERE f.modelID = 2']);
data = fetch(data);
data = data.Data;

scores = cell2mat(data(:,3));
fromunits = cellfun(@str2num, data(:,4));
tounits = cellfun(@str2num, data(:,1));
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
distances = sqrt((fromposx-toposx).^2+(fromposy-toposy).^2);
angles = atan((fromposy-toposy)./(fromposx-toposx));

indices = randsample(length(distances), 2000);
scatter(distances(indices), scores(indices))

figure
scatter(angles(indices), scores(indices));

lscores = log(scores);
smoothed = smooth(distances, scores, .1, 'rloess');
%figure
[xx,ind] = sort(distances);
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