%Compute sum of granger scores for out degree for given unit and nev file
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
modelID = 3;
blackrock = './blackrock/';
%Fetch paramcode to load
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};
eval(paramcode);

%Connections to ipsi units
%Pick a bunch of pairs of BCI neurons from curated dataset
data = exec(conn, ['SELECT bc.`ID`, bc.`unit`, bc2.`unit` FROM `bci_units` bc '...
'INNER JOIN experiment_tuning et '...
'ON et.`1DBCrecording` = bc.`ID` '...
'INNER JOIN `bci_units` bc2 '... 
'ON bc2.`ID` = bc.`ID` AND bc2.`unit` != bc.`unit` '...
'GROUP BY bc.`ID`']);
data = fetch(data);
data = data.Data;
nR = size(data,1);
wgts = 0.05:0.05:.95;
nW = length(wgts);

samplerate = 1/binsize;
sigma_fr = 0.25;
sigma_fr = sigma_fr*samplerate;
sz = sigma_fr*3*2;
x = linspace(-sz/2, sz/2, sz);
gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);

%Estimate based on approximately last 5s of data (assuming binsize is 1/25 s)...
decayrate = 1/125;
sz = 200;
x = linspace(-sz, 0, sz);
expFilter = exp(x*decayrate);
expFilter = expFilter/sum(expFilter);
expFilter = fliplr(expFilter);
results = zeros(nR, nW, 2);
results_method2 = zeros(nR, nW, 2);

momax = 40;
acmaxlags = 2000;

for idx = 1:nR 
	nevfile = data{idx,1};
	bcunit1 = data{idx,2};
	bcunit2 = data{idx,3};
	units = {bcunit1, bcunit2};
	nevpath = [blackrock nevfile];
	%Load nev file
	display(['Loading ' nevfile])
	processed = preprocess_spline(nevpath, binsize, threshold, offset, [], [], units);
	for i = 1:length(wgts)
		weights = [wgts(i), 1-wgts(i)];
		%For each unit:
		del = zeros(size(processed.binnedspikes,1),1);
		cursor = zeros(size(del));
		for j = 1:2
			%Decode cursor with different relative weights according to exponential filter BCI
			%Smooth firing with Gaussian filter
			gf = conv(processed.binnedspikes(:,j), gaussFilter_fr, 'same');
			%Estimate baseline firing with exponential filter
			baseline = conv(processed.binnedspikes(:,j), expFilter);
			baseline = baseline(1:end-sz+1);
			%Test plot
			%ii = 2000:3000;       
			%tt = ii*binsize;
			%plot(tt, processed.binnedspikes(ii,j), tt, baseline(ii), tt, gf(ii));
			del = del+weights(j)*(gf-baseline);
		end
		cursor = cumsum(del);
		%Run MVGC on these cursor and neuron pairs
		X = [cursor, processed.rates]';
		mvgc = mvgcBCI(X, [], [], pval, momax, acmaxlags);
		results(idx,i,:) = mvgc.pwcgc(1,2:3);
	end
end

%Take average of results and plot
averesults = squeeze(mean(results,1));
stdresults = squeeze(std(results,1));
sresults = sum(averesults,2);
presults = averesults(:,1)./sresults;
stdresults = stdresults(:,1)./sresults;
%Note the symmetry!
%plot(wgts, presults, [0 1], [0 1], 'k:', wgts, presults+stdresults, 'r--', wgts, presults-stdresults, 'r--');
plot(wgts, presults, [0 1], [0 1], 'k:');
xlabel('BCI weight of unit A')
ylabel('Proportion of total GC score of unit A')
legend('Simulation')
saveplot(gcf, './worksheets/2015_08_10-cursorsim/cursorsimresults.eps')
%Looks good enough...

%However, the other way to do this is to take, at each time step, a random choice as to whether to
%'listen' to neuron 1 or neuron 2. I'll try this method too:

for idx = 1:nR 
	nevfile = data{idx,1};
	bcunit1 = data{idx,2};
	bcunit2 = data{idx,3};
	units = {bcunit1, bcunit2};
	nevpath = [blackrock nevfile];
	%Load nev file
	display(['Loading ' nevfile])
	processed = preprocess_spline(nevpath, binsize, threshold, offset, [], [], units);
	for i = 1:length(wgts)
		weight = wgts(i);
		del = zeros(size(processed.binnedspikes,1),1);
		cursor = zeros(size(del));
		draws = rand(size(del));
		unita = draws<weight;
		unitb = draws>=weight;
		%Unit a's contributions
		gf = conv(processed.binnedspikes(:,1), gaussFilter_fr, 'same');
		baseline = conv(processed.binnedspikes(:,1), expFilter);
		baseline = baseline(1:end-sz+1);
		del = del+unita.*(gf-baseline);
		%Unit b's contributions
		gf = conv(processed.binnedspikes(:,2), gaussFilter_fr, 'same');
		baseline = conv(processed.binnedspikes(:,2), expFilter);
		baseline = baseline(1:end-sz+1);
		del = del+unitb.*(gf-baseline);
		%Sum them
		cursor = cumsum(del);
		%Run MVGC on these cursor and neuron pairs
		X = [cursor, processed.rates]';
		mvgc = mvgcBCI(X, [], [], pval, momax, acmaxlags);
		results_method2(idx,i,:) = mvgc.pwcgc(1,2:3);
	end
end