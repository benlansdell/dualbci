matfilein = './expts/2dmanualpos_1day.mat';

load(matfilein);

%Load all condition info
conds = {'1D Horiz Brain  Position',    '1D Horiz Brain  Velocity',    '1D Horiz Manual Position',    '1D Horiz Manual Velocity',    '1D Vert  Brain  Position',    '1D Vert  Brain  Velocity',    '1D Vert  Manual Position',    '1D Vert  Manual Velocity',    '2D Brain  Velocity',    '2D Manual Position',    '2D Manual Velocity', 'Dual Control'};
nC = length(conds);
nE = length(expts);

%Extract total time for each condition
totalcond_secs = zeros(1,nC);
eachcond_secs = zeros(nE, nC);
for idx = 1:nE
	totalcond_secs = totalcond_secs + expts(idx).cond_secs;
	eachcond_secs(idx,:) = expts(idx).cond_secs;
end

%Plot bar plot of all recordings worth
clf
bar(totalcond_secs/3600);
xlabel('Condition');
ylabel('Hours of recording');
set(gca, 'XTickLabel',conds, 'XTick',1:numel(conds))
%Rotate xlabels
rotateXLabels(gca, 90)
%scale = 0.9;
%pos = get(gca, 'Position');
%pos(2) = pos(2)+scale*pos(4);
%pos(4) = (1-scale)*pos(4);
%set(gca, 'Position', pos)
saveplot(gcf, './worksheets/02_18_2015/plots/2dmanualpos_1day_totalconds.eps', 'eps', [6 12])

%Plot histogram for each condition
for idx = 1:nC
	clf
	nz = eachcond_secs(:,idx) ~= 0;
	hist(eachcond_secs(nz,idx)/60);
	xlabel('Minutes of recording')
	ylabel('Count')
	title(['Condition: ' conds{idx} ' nonzero expts=' num2str(sum(nz))])
	saveplot(gcf, ['./worksheets/02_18_2015/plots/2dmanualpos_1day_conds_' num2str(idx) '.eps'])
end