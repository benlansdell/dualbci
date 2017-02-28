nevfile = './testdata/20130117SpankyUtah001.nev';
binsize = 0.05;
offset = -0.075;
threshold = 5;
verbose = 0;

nK_sp = 1;
nK_pos = 1;
processed = preprocess(nevfile, binsize, threshold, offset);
data = filters_sp_vel(processed, nK_sp, nK_pos);

%Plot polar plot of firing rate vs angle

nB = 10000;
nU = size(data.y,1);
nbins = 10;

firing = zeros(nU,nbins);

for nrn = 1:nU
	a = {};
	for n = 1:nbins
		a{n} = [];
	end 
	for idx = 1:nB
		%Compute angle
		an = atan(data.dtorque(idx,2)/data.dtorque(idx,1));
		if data.dtorque(idx,2) < 0
			an = an + pi;
		end
		if data.dtorque(idx,1) < 0 & data.dtorque(idx,2)>0
			an = an + 2*pi;
		end
		bi = ceil(an/(2*pi)*nbins);
	
		%Compute rate
		ra = data.y(nrn, idx)/binsize;
		a{bi}(end+1) = ra;
	end 
	
	%plot this baby
	for n = 1:nbins
		firing(nrn, n) = mean(a{n});
	end
end

nx = 4;
ny = 6;

firing(:,nbins+1) = firing(:,1);

lnsp = linspace(0,2*pi, nbins+1);
for idx = 1:nU
	subplot(nx,ny,idx);
	polarplot(lnsp, firing(idx,:));
	title(num2str(idx))
end

%Just plot 4 of them... 
toplot = [5, 4, 2, 15];
figure
for idx = 1:length(toplot)
	i = toplot(idx);
	subplot(2,2,idx);
	polarplot(lnsp, firing(i,:));
	title(num2str(processed.unitnames{i}))
end	

saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/tuning_example.eps')