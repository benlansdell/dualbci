%Find all units common between start and end recording

%%%%%%%
%Setup%
%%%%%%%

conds = {'1D Horiz Brain  Position',    '1D Horiz Brain  Velocity',    '1D Horiz Manual Position',    '1D Horiz Manual Velocity',    '1D Vert  Brain  Position',    '1D Vert  Brain  Velocity',    '1D Vert  Manual Position',    '1D Vert  Manual Velocity',    '2D Brain  Velocity',    '2D Manual Position',    '2D Manual Velocity', 'Dual Control'};
mat_sprc_pos = load('./expts/2dmanualpos_sprc_pos_lv_def.mat');
mat_sprc = load('./expts/2dmanualpos_sprc_def.mat');
expts_sprc_pos = mat_sprc_pos.expts;
expts_sprc = mat_sprc.expts;

nE = length(expts_sprc);
nC = length(conds);
nEl = 128;

DP = zeros(nEl, nE, 4+nC+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Note their change in deviance between sprc and sprc_pos_lv%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For each unit, find all experiments containing that unit
for i = 1:nE
	expt_sprc = expts_sprc(i);
	expt_sprc_pos = expts_sprc_pos(i);
	unitsstart = cellfun(@str2num, expt_sprc_pos.modelstart.unitnames);
	unitsend = cellfun(@str2num, expt_sprc_pos.modelend.unitnames);
	daten = expt_sprc.startdate;
	%Find units in common in both start and end of expt
	%Append this experiment's info (deviances for each case, and the seconds of recording in between)
	for j = 1:nEl
		startidx = find(unitsstart == j);
		endidx = find(unitsend == j);
		if any(unitsstart == j) & any(unitsend == j)
			%Pull out events
			evts = expt_sprc.events;
			bcisecs = 0;
			for k = 1:length(evts)
				evt = evts{k};
				nsettings = length(evt.dur);
				for idx = 1:nsettings
					%See if current electrode was used in this settings, if it was, add to BCI duration
					if any(evt.bciunits{idx} == j)
						bcisecs = bcisecs + evt.dur(idx);
					end
				end
			end
			%Find number of hours this unit was used as a BCI
			DP(j, i, 1) = expt_sprc.modelstart.dev{startidx};
			DP(j, i, 2) = expt_sprc.modelend.dev{endidx};
			DP(j, i, 3) = expt_sprc_pos.modelstart.dev{startidx};
			DP(j, i, 4) = expt_sprc_pos.modelend.dev{endidx};
			DP(j, i, 5:16) = expt_sprc.cond_secs/3600;
			DP(j, i, 17) = bcisecs/3600;
			DP(j, i, 18) = daten;
		end
	end
end

%%%%%%
%Plot%
%%%%%%
totidx = 5:12;
manidx = 4+[3 4 7 8 10 11];
braidx = 4+[1 2 5 6 9 12];
bciidx = 4+[13];
datidx = 4+[14];

for idx = 1:nEl
	nz = squeeze(DP(idx,:,1))>0;
	%Only plot electrodes that were used during at least one experiment
	if sum(nz) > 0
		display(['Plotting ' num2str(idx) ' n=' num2str(sum(nz))])
		nzExpts = squeeze(DP(idx,nz,:));
		if size(nzExpts, 2)==1
			nzExpts = nzExpts';
		end
		%Scatter plot showing DeltaD_p before (xaxis) and after (yaxis) for all expts
		%Color this plot in different ways to highlight (hopefully) interesting features
		x = nzExpts(:,2)-nzExpts(:,4);
		y = nzExpts(:,1)-nzExpts(:,3);
		ct = sum(nzExpts(:,totidx),2);
		cm = sum(nzExpts(:,manidx),2);
		cb = sum(nzExpts(:,braidx),2);
		cc = sum(nzExpts(:,bciidx),2);
		da = nzExpts(:,datidx);
		xmax = min(max(x), 2000);
		ymax = min(max(y), 2000);
		tt = linspace(0,max(xmax,ymax), 100);
		n = length(x);
		%Color according to total expt time
		clf
		subplot(4,1,1);
		scatter(x,y,[],ct);	
		colorbar
		hold on
		plot(tt,tt,'--')
		xlim([0 xmax])
		ylim([0 ymax])
		ylabel('Day 0')
		title(['Unit: ' num2str(idx) ' n=' num2str(n) ' Total recorded time.'])
		%Color according to manual control time
		subplot(4,1,2);
		scatter(x,y,[],cm);	
		colorbar
		hold on
		plot(tt,tt,'--')
		xlim([0 xmax])
		ylim([0 ymax])
		ylabel('Day 0')
		title(['Manual control time'])
		%Color according to brain control time
		subplot(4,1,3);
		scatter(x,y,[],cb);	
		colorbar
		hold on
		plot(tt,tt,'--')
		xlim([0 xmax])
		ylim([0 ymax])
		ylabel('Day 0')
		title(['Brain control time'])
		%Color according to BCI unit time specifically
		subplot(4,1,4);
		scatter(x,y,[],cc);	
		colorbar
		hold on
		plot(tt,tt,'--')
		xlim([0 xmax])
		ylim([0 ymax])
		xlabel('Day 1')
		ylabel('Day 0')
		title(['Time driving BCI'])		
		saveplot(gcf, ['./worksheets/02_19_2015/plots/2dmanualpos_1day_sprc-sprcpos_unit_' num2str(idx) '_n_' num2str(n) '.eps'], 'eps', [5 20])
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make the same plots but with all units' information aggregated%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = [];
Y = [];
CT = [];
CM = [];
CB = [];
CC = [];
N = 0;

for idx = 1:nEl
	nz = squeeze(DP(idx,:,1))>0;
	%Only plot electrodes that were used during at least one experiment
	if sum(nz) > 0
		display(['Plotting ' num2str(idx) ' n=' num2str(sum(nz))])
		nzExpts = squeeze(DP(idx,nz,:));
		if size(nzExpts, 2)==1
			nzExpts = nzExpts';
		end
		%Scatter plot showing DeltaD_p before (xaxis) and after (yaxis) for all expts
		%Color this plot in different ways to highlight (hopefully) interesting features
		x = nzExpts(:,2)-nzExpts(:,4);
		y = nzExpts(:,1)-nzExpts(:,3);
		ct = sum(nzExpts(:,totidx),2);
		cm = sum(nzExpts(:,manidx),2);
		cb = sum(nzExpts(:,braidx),2);
		cc = sum(nzExpts(:,bciidx),2);
		da = nzExpts(:,datidx);
		n = length(x);
		X = [X(:); x];
		Y = [Y(:); y];
		CT = [CT(:); ct];
		CM = [CM(:); cm];
		CB = [CB(:); cb];
		CC = [CC(:); cc];
		N = N + n;
	end
end

xmax = min(max(X), 2000);
ymax = min(max(Y), 2000);
TT = linspace(0,max(xmax,ymax), 100);
%Color according to total expt time
clf
subplot(4,1,1);
scatter(X,Y,[],CT);	
xlim([0 xmax])
ylim([0 ymax])
colorbar
hold on
plot(TT,TT,'--')
ylabel('Day 0')
title(['All units. N=' num2str(N) ' Total recorded time.'])
%Color according to manual control time
subplot(4,1,2);
scatter(X,Y,[],CM);	
xlim([0 xmax])
ylim([0 ymax])
colorbar
hold on
plot(TT,TT,'--')
ylabel('Day 0')
title(['Manual control time'])
%Color according to brain control time
subplot(4,1,3);
scatter(X,Y,[],CB);	
xlim([0 xmax])
ylim([0 ymax])
colorbar
hold on
plot(TT,TT,'--')
ylabel('Day 0')
title(['Brain control time'])
%Color according to BCI unit time specifically
subplot(4,1,4);
scatter(X,Y,[],CC);	
xlim([0 xmax])
ylim([0 ymax])
colorbar
hold on
plot(TT,TT,'--')
xlabel('Day 1')
ylabel('Day 0')
title(['Time driving BCI'])		
saveplot(gcf, ['./worksheets/02_19_2015/plots/2dmanualpos_1day_sprc-sprcpos_all_N_' num2str(N) '.eps'], 'eps', [5 20])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot as a function of date%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx = 1:nEl
	nz = squeeze(DP(idx,:,1))>0;
	%Only plot electrodes that were used during at least one experiment
	if sum(nz) > 0
		display(['Plotting ' num2str(idx) ' n=' num2str(sum(nz))])
		nzExpts = squeeze(DP(idx,nz,:));
		if size(nzExpts, 2)==1
			nzExpts = nzExpts';
		end
		%Scatter plot showing DeltaD_p before (xaxis) and after (yaxis) for all expts
		%Color this plot in different ways to highlight (hopefully) interesting features
		x = nzExpts(:,2)-nzExpts(:,4);
		y = nzExpts(:,1)-nzExpts(:,3);
		z = x-y;

		ct = sum(nzExpts(:,totidx),2);
		cm = sum(nzExpts(:,manidx),2);
		cb = sum(nzExpts(:,braidx),2);
		cc = sum(nzExpts(:,bciidx),2);
		da = nzExpts(:,datidx);
		%ymin = max(min(z), -2000);
		%ymax = min(max(z), 2000);
		ymin = -1000;
		ymax = 1000;
		n = length(x);
		%Color according to total expt time
		clf
		subplot(4,1,1);
		scatter(da,z,[],ct);	
		colorbar
		ylim([ymin ymax])
		datetick('x', 'mm-dd', 'keeplimits')
		ylabel('DP')
		title(['Unit: ' num2str(idx) ' n=' num2str(n) ' Total recorded time.'])
		%Color according to manual control time
		subplot(4,1,2);
		scatter(da,z,[],cm);	
		datetick('x', 'mm-dd', 'keeplimits')
		colorbar
		ylim([ymin ymax])
		ylabel('DP')
		title(['Manual control time'])
		%Color according to brain control time
		subplot(4,1,3);
		scatter(da,z,[],cb);	
		datetick('x', 'mm-dd', 'keeplimits')
		colorbar
		ylim([ymin ymax])
		ylabel('DP')
		title(['Brain control time'])
		%Color according to BCI unit time specifically
		subplot(4,1,4);
		scatter(da,z,[],cc);	
		datetick('x', 'mm-dd', 'keeplimits')
		colorbar
		ylim([ymin ymax])
		xlabel('Date (2013)')
		ylabel('DP')
		title(['Time driving BCI'])		
		saveplot(gcf, ['./worksheets/02_19_2015/plots/2dmanualpos_1day_sprc-sprcpos_date_unit_' num2str(idx) '_n_' num2str(n) '.eps'], 'eps', [5 20])
	end
end

X = [];
Y = [];
CT = [];
CM = [];
CB = [];
CC = [];
DA = [];
N = 0;

for idx = 1:nEl
	nz = squeeze(DP(idx,:,1))>0;
	%Only plot electrodes that were used during at least one experiment
	if sum(nz) > 0
		display(['Plotting ' num2str(idx) ' n=' num2str(sum(nz))])
		nzExpts = squeeze(DP(idx,nz,:));
		if size(nzExpts, 2)==1
			nzExpts = nzExpts';
		end
		%Scatter plot showing DeltaD_p before (xaxis) and after (yaxis) for all expts
		%Color this plot in different ways to highlight (hopefully) interesting features
		x = nzExpts(:,2)-nzExpts(:,4);
		y = nzExpts(:,1)-nzExpts(:,3);
		ct = sum(nzExpts(:,totidx),2);
		cm = sum(nzExpts(:,manidx),2);
		cb = sum(nzExpts(:,braidx),2);
		cc = sum(nzExpts(:,bciidx),2);
		da = nzExpts(:,datidx);
		n = length(x);
		X = [X(:); x];
		Y = [Y(:); y];
		CT = [CT(:); ct];
		CM = [CM(:); cm];
		CB = [CB(:); cb];
		CC = [CC(:); cc];
		DA = [DA(:); da];
		N = N + n;
	end
end

ymin = -1000;
ymax = 1000;
Z = X-Y;
%Color according to total expt time
clf
subplot(4,1,1);
scatter(DA,Z,[],CT);	
ylim([ymin ymax])
datetick('x', 'mm-dd', 'keeplimits')
colorbar
ylabel('DP')
title(['All units. N=' num2str(N) ' Total recorded time.'])
%Color according to manual control time
subplot(4,1,2);
scatter(DA,Z,[],CM);	
ylim([ymin ymax])
datetick('x', 'mm-dd', 'keeplimits')
colorbar
ylabel('DP')
title(['Manual control time'])
%Color according to brain control time
subplot(4,1,3);
scatter(DA,Z,[],CB);	
ylim([ymin ymax])
datetick('x', 'mm-dd', 'keeplimits')
colorbar
ylabel('DP')
title(['Brain control time'])
%Color according to BCI unit time specifically
subplot(4,1,4);
scatter(DA,Z,[],CC);	
ylim([ymin ymax])
datetick('x', 'mm-dd', 'keeplimits')
colorbar
xlabel('Date (2013)')
ylabel('DP')
title(['Time driving BCI'])		
saveplot(gcf, ['./worksheets/02_19_2015/plots/2dmanualpos_1day_sprc-sprcpos_all_date_N_' num2str(N) '.eps'], 'eps', [5 20])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot change in DP between Day 0 and Day 1%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for idx = 1:nEl
	nz = squeeze(DP(idx,:,1))>0;
	%Only plot electrodes that were used during at least one experiment
	if sum(nz) > 0
		display(['Plotting ' num2str(idx) ' n=' num2str(sum(nz))])
		nzExpts = squeeze(DP(idx,nz,:));
		if size(nzExpts, 2)==1
			nzExpts = nzExpts';
		end
		%Scatter plot showing DeltaD_p before (xaxis) and after (yaxis) for all expts
		%Color this plot in different ways to highlight (hopefully) interesting features
		x = nzExpts(:,2)-nzExpts(:,4);
		y = nzExpts(:,1)-nzExpts(:,3);
		z = x-y;

		ct = sum(nzExpts(:,totidx),2);
		cm = sum(nzExpts(:,manidx),2);
		cb = sum(nzExpts(:,braidx),2);
		cc = sum(nzExpts(:,bciidx),2);
		da = nzExpts(:,datidx);
		%ymin = max(min(z), -2000);
		%ymax = min(max(z), 2000);
		ymin = -1000;
		ymax = 1000;
		n = length(x);
		%Color according to total expt time
		clf
		subplot(4,1,1);
		scatter(da,z,[],ct);	
		colorbar
		ylim([ymin ymax])
		datetick('x', 'mm-dd', 'keeplimits')
		ylabel('DP')
		title(['Unit: ' num2str(idx) ' n=' num2str(n) ' Total recorded time.'])
		%Color according to manual control time
		subplot(4,1,2);
		scatter(da,z,[],cm);	
		datetick('x', 'mm-dd', 'keeplimits')
		colorbar
		ylim([ymin ymax])
		ylabel('DP')
		title(['Manual control time'])
		%Color according to brain control time
		subplot(4,1,3);
		scatter(da,z,[],cb);	
		datetick('x', 'mm-dd', 'keeplimits')
		colorbar
		ylim([ymin ymax])
		ylabel('DP')
		title(['Brain control time'])
		%Color according to BCI unit time specifically
		subplot(4,1,4);
		scatter(da,z,[],cc);	
		datetick('x', 'mm-dd', 'keeplimits')
		colorbar
		ylim([ymin ymax])
		xlabel('Date (2013)')
		ylabel('DP')
		title(['Time driving BCI'])		
		saveplot(gcf, ['./worksheets/02_19_2015/plots/2dmanualpos_1day_sprc-sprcpos_date_unit_' num2str(idx) '_n_' num2str(n) '.eps'], 'eps', [5 20])
	end
end

%%%%%%%%%
%Fit LMs%
%%%%%%%%%

%Do some stats on this...