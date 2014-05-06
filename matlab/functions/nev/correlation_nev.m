function [scoreFE, scoreRU] = correlation_nev(nevfile, fn_out, threshold, binsize, sigma_fr, sigma_trq)
	%correlation_nev	Function to fit the following model to spike and torque data:
	%			lambda_i(t) = \lambda_0 + \sum_j^N k_j^1 . x_i^1(t+\tau+jh) + \sum_j^N k_j^2 . x_i^2(t+\tau+jh)
	%		That is, it fits a linear filter to the torque data. In the above formula, time-step size h
	%		is given by parameter binsize; max value of N is given by kernellength; \tau is given by offset.
	%		Program fits all models with kernels of length between 1 and N, meaning that for each single-unit N models
	%		will be fit.
	%
	% 		Input:
	%			nevfile = file to process
	%			fn_out = base name for output plots
	%			threshold = (optional, default = 5) threshold firing rate below which unit is ignored
	%			binsize = (optional, default = 0.002) time bin size for spikes (torque resmapled to this rate0)
	%			sigma_fr = (optional, default = 0) width of gaussian filter to apply to spikes for firing rate. If 0 then no filter applied
	%			sigma_trq = (optional, default = 0.25) width of gaussian filter to apply to torque. If 0 then no filter applied
	%				Note: both sigmas are in units of seconds, and then are scaled according to binsize
	%		
	%		Output:
	%			(none) produces plots of the filter for each single-unit channel, along with summary of quality of fits for 
	%			each channel
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
	%			threshold = 5;
	%			binsize = 0.002;
	%			sigma_fr = 0; 
	%			sigma_trq = 0.25;
	%			fn_out = './worksheets/tuning/crosscorr/20130117SpankyUtah001';
	%			correlation_nev(nevfile, fn_out, threshold, binsize, sigma_fr, sigma_trq);
	
	%Optional arguments
	if (nargin < 3)	threshold = 5; end
	if (nargin < 4) binsize = 0.002; end
	if (nargin < 5) sigma_fr = 0; end
	if (nargin < 6) sigma_trq = 0.25; end
	%Set this to above 0 if want debug info/plots
	verbosity = 1;
	%Total number of possible units recorded from
	nE = 128;
	nunits = 5; 
	nU = nE*nunits;
	%Offset to apply
	%offset = 0;
	%Put firing rate ahead of torque
	offset = -0.150;
	%Size of gaussian filter to apply
	samplerate = 1/binsize;
	%Max lag for correlations
	maxlag = 90;
	maxlag = round(maxlag*samplerate)/samplerate;
	maxpeak = 3;
	%Make sure we can perform the sample rate conversion easily
	assert(rem(samplerate,1) == 0, 'Select a binsize corresponding to an integer sample rate.');
	ns3file = [nevfile(1:end-3) 'ns3'];

	%Filters
    if sigma_fr > 0
    	sigma_fr = sigma_fr*samplerate;
   		sz = sigma_fr*3;
	    x = linspace(-sz/2, sz/2, sz);
	    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    	gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    end
    if sigma_trq > 0
    	sigma_trq = sigma_trq*samplerate;
   		sz = sigma_trq*3;
	    x = linspace(-sz/2, sz/2, sz);
	    gaussFilter_trq = exp(-x.^2/(2*sigma_trq^2));
    	gaussFilter_trq = gaussFilter_trq/sum(gaussFilter_trq);
    end

	%%%%%%%%%%%%%%%%%%%%%%
	%Process spiking data%
	%%%%%%%%%%%%%%%%%%%%%%
	NEV = openNEV(nevfile);
	%Find the duration and sample rate of the nev file recording
	nevsamplerate = NEV.MetaTags.TimeRes;
	dur = NEV.MetaTags.DataDuration/nevsamplerate;
	%Convert spike times into array of binned spikes, one for each spike sorted channel
	spiketimes = double(NEV.Data.Spikes.TimeStamp)/nevsamplerate;
    elecs = cell(1,nU);
    spikemuas = struct('times', elecs);
    unitnames = cell(1,nU);
    averate = zeros(1,nU);
    for idx=1:nU
        spikemuas(idx).times = [0];    
    end
    for i=1:length(spiketimes)
       	E = NEV.Data.Spikes.Electrode(i);
       	unit = NEV.Data.Spikes.Unit(i);
       	U = single((E-1)*nunits)+single(unit)+1;
       	spikemuas(U).times = [spikemuas(U).times; spiketimes(i)];
       	unitnames{U} = [num2str(E) '.' num2str(unit)];
    end
    %Check which channels are doing stuff
    for idx=1:nU
    	averate(idx) = (length(spikemuas(idx).times)-1)/dur;
    	if (verbosity) display(['Electrode.Unit: ' unitnames{idx} ' Spike count: ' num2str(length(spikemuas(idx).times)-1) ' Mean firing rate (Hz): ' num2str(averate(idx))]); end
    end
    %Set a threshold firing rate, below which we ignore that unit
    abovethresh = averate > threshold;
    %Update nU
    nU = sum(abovethresh);
    unitnames = unitnames(abovethresh);
    spikemuas = spikemuas(abovethresh);
    averate = averate(abovethresh);
   	%Bin spikes (chronux function)
    binnedspikes = binspikes(spikemuas, samplerate);
    %From this apply gaussian filter to spike train for each electrode
    for idx=1:nU
    	if sigma_fr > 0
    	    gf = conv(binnedspikes(:,idx), gaussFilter_fr, 'same');
			rates(:,idx)=gf*samplerate;
	    else
    	    rates(:,idx)=binnedspikes(:,idx)*samplerate;
	    end
    end
    %rates

	%%%%%%%%%%%%%%%%%%%%%
	%Process torque data%
	%%%%%%%%%%%%%%%%%%%%%
	NS3 = openNSx(ns3file, 'read', 'c:138:139');
	nsxtorque = double(NS3.Data);
	nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
	%Switch sign of FE axis for coordinate consistency
	nsxtorque(2,:)=-nsxtorque(2,:);
    for j=1:2
        %Scale from uint16 value to proportion
        nsxtorque(j,:) = nsxtorque(j,:)/(2^15);
        %Subtract mean
        nsxtorque(j,:) = nsxtorque(j,:)-mean(nsxtorque(j,:));
        %Smooth
        if sigma_trq > 0
	        nsxtorque(j,:) = conv(nsxtorque(j,:),gaussFilter_trq,'same');
	    end
    end
	%Resample at rate of binsize
	torque = resample(nsxtorque', samplerate, nsxsamplerate);	
    plot(torque(1:(5*samplerate),1), torque(1:(5*samplerate),2));	
	%Check they're the same length, and trim
	nsamp = min(size(torque,1), size(rates,1));
	torque=torque(1:nsamp,:);
	rates = rates(1:nsamp,:);
	%Apply offset to data
	delaysamples = round(offset*samplerate);
	%If offset is positive, then firing rate proceeds torque
	if (delaysamples > 0)
    	binnedspikes = binnedspikes(1+delaysamples:end,:);
    	rates = rates(1+delaysamples:end,:);
    	torque = torque(1:end-delaysamples,:);
    %If offset is negative then firing rate precedes torque
	elseif (delaysamples < 0)
	    binnedspikes = binnedspikes(1:end+delaysamples,:);
	    rates = rates(1:end+delaysamples,:);
    	torque = torque(1-delaysamples:end,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute p(x,y|spike)/p(x,y) ~ p(spike|x,y)%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dtorquex = [diff(torque(:,1)); 0]; dtorquey = [diff(torque(:,2)); 0];
	ddtorquex = [diff(dtorquex); 0]; ddtorquey = [diff(dtorquey); 0];
	dtorquex(1:5) = 0; dtorquey(1:5) = 0;
	ddtorquex(1:5) = 0; ddtorquey(1:5) = 0;
	figure 
	subplot(2,2,1)
    [ctrs1,ctrs2,nc, priorF] = smoothhist2D([torque(:,1), torque(:,2)], 5, [100, 100], 0.05);
    title('Torque position density')
	subplot(2,2,2)
    [ctrs1d,ctrs2d,ncd, priorFd] = smoothhist2D([dtorquex, dtorquey], 5, [100, 100], 0.05);
    title('Torque velocity density')
	subplot(2,2,3)
    [ctrs1dd,ctrs2dd,ncdd, priorFdd] = smoothhist2D([ddtorquex, ddtorquey], 5, [100, 100], 0.05);
    title('Torque accel density')
    saveplot(gcf, [fn_out '_torque_density.eps'], 'eps', [2 6]);
    for i=1:nU
    	spikeidx = find(binnedspikes(:,i)>0);
    	figure
    	%Position
    	subplot(2,2,1)
     	[ctrs1c,ctrs2c,nc, relF] = smoothhist2D([torque(spikeidx,1), torque(spikeidx,2)], 5, [100, 100], 0.05);
     	normF = relF./priorF;
     	normFscaled = normF/max(max(normF(:)));
     	image(ctrs1,ctrs2,floor(nc.*normFscaled) + 1);
%     	image(ctrs1,ctrs2,floor(nc.*normF) + 1);
    	title(unitnames{i})
%    	xlim([-0.5 0.5])
%    	ylim([-0.5 0.5])
    	%Velocity
    	subplot(2,2,2)
     	[ctrs1cd,ctrs2cd,nc, relF] = smoothhist2D([dtorquex(spikeidx), dtorquey(spikeidx)], 5, [100, 100], 0.05);
     	normF = relF./priorFd;
     	normFscaled = normF/max(max(normF(:)));
     	image(ctrs1d,ctrs2d,floor(nc.*normFscaled) + 1);
%     	image(ctrs1d,ctrs2d,floor(nc.*normF) + 1);
	   	%xlim([-0.02 0.02])
    	%ylim([-0.02 0.02])
    	%Accel
    	subplot(2,2,3)
     	[ctrs1cdd,ctrs2cdd,nc, relF] = smoothhist2D([ddtorquex(spikeidx), ddtorquey(spikeidx)], 5, [100, 100], 0.05);
     	normF = relF./priorFdd;
     	normFscaled = normF/max(max(normF(:)));
		image(ctrs1dd,ctrs2dd,floor(nc.*normFscaled) + 1);
%		image(ctrs1dd,ctrs2dd,floor(nc.*normF) + 1);
    	title(unitnames{i})
    	%xlim([-0.02 0.02])
    	%ylim([-0.02 0.02])
   		saveplot(gcf, [fn_out '_spikedensity_unit_' unitnames{i} '.eps'], 'eps', [5 5]);
%    	pause
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Cross- and auto-correlation%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tt = -maxlag:binsize:maxlag;
	peakrange = find(tt > -maxpeak & tt < maxpeak);    % area to search for peak
    peakRU = zeros(1,nU);
    peakFE = zeros(1,nU);
    lagRU = zeros(1,nU);
    lagFE = zeros(1,nU);
    stdRU = zeros(1,nU);
    stdFE = zeros(1,nU);
    scoreRU = zeros(1,nU);
    scoreFE = zeros(1,nU);

   	autotorqueFE = xcov(torque(:,1),samplerate*maxlag);%, 'coeff');
   	autotorqueRU = xcov(torque(:,2),samplerate*maxlag);%, 'coeff');

	figure 
	subplot(1,2,1);
   	plot(tt,autotorqueFE);
   	xlim([-maxpeak*2 maxpeak*2])
   	title('Auto-corr torque FE');
	subplot(1,2,2);
   	plot(tt,autotorqueRU);
   	xlim([-maxpeak*2 maxpeak*2])
   	title('Auto-corr torque RU');
	saveplot(gcf, [fn_out '_auto-torque.eps'], 'eps', [4 2]);

    for i=1:nU
    	%Compute cross correlation
	    covFE = xcov(rates(:,i), torque(:,1),samplerate*maxlag,'unbiased');
    	% normalize against spikes auto-covariance
    	autorate = xcov(rates(:,i),samplerate*maxlag);%, 'coeff');
    	covFE = covFE / sqrt(xcov(rates(:,i),0));
    	covFE = covFE / sqrt(xcov(torque(:,1),0));
    	peakFE(i) = covFE(abs(covFE) == max(abs(covFE(peakrange))));
	    stdFE(i) = std(covFE);
    	lagFE(i) = tt((covFE == peakFE(i)) | (covFE == -peakFE(i)));
    	avg = mean(covFE);
    	peakFE(i) = peakFE(i) - avg;
    	scoreFE(i) = peakFE(i)/stdFE(i);
    	%Plot cross correlation
    	figure
    	subplot(3,2,1);
    	plot(tt,covFE,[-maxlag maxlag],[avg avg],...
               [lagFE(i) lagFE(i)],[avg avg+peakFE(i)],...
                  [-maxlag maxlag],avg+[stdFE(i) stdFE(i)]*sign(peakFE(i)));
		title(['FE Unit: ' unitnames{i} ' Score: ' num2str(scoreFE(i))])
		xlim([-maxlag maxlag])
    	subplot(3,2,3);
    	plot(tt,covFE,[-maxlag maxlag],[avg avg],...
               [lagFE(i) lagFE(i)],[avg avg+peakFE(i)],...
                  [-maxlag maxlag],avg+[stdFE(i) stdFE(i)]*sign(peakFE(i)));
		xlim([-maxpeak maxpeak])

		subplot(3,2,5);
    	plot(tt,autorate);
    	xlim([-maxpeak*2 maxpeak*2])
    	title('Auto-corr rate');

	    covRU = xcov(rates(:,i), torque(:,2),samplerate*maxlag,'unbiased');
    	% normalize against spikes auto-covariance
    	covRU = covRU / sqrt(xcov(rates(:,i),0));
    	covRU = covRU / sqrt(xcov(torque(:,2),0));
    	peakRU(i) = covRU(abs(covRU) == max(abs(covRU(peakrange))));
    	stdRU(i) = std(covRU);
    	lagRU(i) = tt((covRU == peakRU(i)) | (covRU == -peakRU(i)));
    	avg = mean(covRU);
    	peakRU(i) = peakRU(i) - avg;
    	scoreRU(i) = peakRU(i)/stdRU(i);
    	%Plot cross- and auto-correlations
    	subplot(3,2,2);
    	plot(tt,covRU,[-maxlag maxlag],[avg avg],...
               [lagRU(i) lagRU(i)],[avg avg+peakRU(i)],...
                  [-maxlag maxlag],avg+[stdRU(i) stdRU(i)]*sign(peakRU(i)));
		title(['RU Score: ' num2str(scoreRU(i))])
		xlim([-maxlag maxlag])
    	subplot(3,2,4);
    	plot(tt,covRU,[-maxlag maxlag],[avg avg],...
               [lagRU(i) lagRU(i)],[avg avg+peakRU(i)],...
                  [-maxlag maxlag],avg+[stdRU(i) stdRU(i)]*sign(peakRU(i)));
		xlim([-maxpeak maxpeak])

		saveplot(gcf, [fn_out '_unit_' unitnames{i} '_cross_maxscore_' num2str(max(scoreFE(i), scoreRU(i))) '.eps'], 'eps', [6 6]);
%		pause
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Cross-correlation with torque velocity and firing rate%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    peakDRU = zeros(1,nU);
    peakDFE = zeros(1,nU);
    lagDRU = zeros(1,nU);
    lagDFE = zeros(1,nU);
    stdDRU = zeros(1,nU);
    stdDFE = zeros(1,nU);
    scoreDRU = zeros(1,nU);
    scoreDFE = zeros(1,nU);

    for i=1:nU
    	%Compute cross correlation
	    covDFE = xcov(rates(:,i), dtorquex,samplerate*maxlag,'unbiased');
    	% normalize against spikes auto-covariance
    	covDFE = covDFE / sqrt(xcov(rates(:,i),0));
    	covDFE = covDFE / sqrt(xcov(dtorquex,0));
    	peakDFE(i) = covDFE(abs(covDFE) == max(abs(covDFE(peakrange))));
	    stdDFE(i) = std(covDFE);
    	lagDFE(i) = tt((covDFE == peakDFE(i)) | (covDFE == -peakDFE(i)));
    	avg = mean(covDFE);
    	peakDFE(i) = peakDFE(i) - avg;
    	scoreDFE(i) = peakDFE(i)/stdDFE(i);
    	%Plot cross correlation
    	figure
    	subplot(3,2,1);
    	plot(tt,covDFE,[-maxlag maxlag],[avg avg],...
               [lagDFE(i) lagDFE(i)],[avg avg+peakDFE(i)],...
                  [-maxlag maxlag],avg+[stdDFE(i) stdDFE(i)]*sign(peakDFE(i)));
		title(['FE Unit: ' unitnames{i} ' Score: ' num2str(scoreDFE(i))])
		xlim([-maxlag maxlag])
    	subplot(3,2,3);
    	plot(tt,covDFE,[-maxlag maxlag],[avg avg],...
               [lagDFE(i) lagDFE(i)],[avg avg+peakDFE(i)],...
                  [-maxlag maxlag],avg+[stdDFE(i) stdDFE(i)]*sign(peakDFE(i)));
		xlim([-maxpeak maxpeak])

	    covDRU = xcov(rates(:,i), dtorquey,samplerate*maxlag,'unbiased');
    	% normalize against spikes auto-covariance
    	covDRU = covDRU / sqrt(xcov(rates(:,i),0));
    	covDRU = covDRU / sqrt(xcov(dtorquey,0));
    	peakDRU(i) = covDRU(abs(covDRU) == max(abs(covDRU(peakrange))));
    	stdDRU(i) = std(covDRU);
    	lagDRU(i) = tt((covDRU == peakDRU(i)) | (covDRU == -peakDRU(i)));
    	avg = mean(covDRU);
    	peakDRU(i) = peakDRU(i) - avg;
    	scoreDRU(i) = peakDRU(i)/stdRU(i);
    	%Plot cross- and auto-correlations
    	subplot(3,2,2);
    	plot(tt,covDRU,[-maxlag maxlag],[avg avg],...
               [lagDRU(i) lagDRU(i)],[avg avg+peakDRU(i)],...
                  [-maxlag maxlag],avg+[stdDRU(i) stdDRU(i)]*sign(peakDRU(i)));
		title(['RU Score: ' num2str(scoreDRU(i))])
		xlim([-maxlag maxlag])
    	subplot(3,2,4);
    	plot(tt,covDRU,[-maxlag maxlag],[avg avg],...
               [lagDRU(i) lagDRU(i)],[avg avg+peakDRU(i)],...
                  [-maxlag maxlag],avg+[stdDRU(i) stdDRU(i)]*sign(peakDRU(i)));
		xlim([-maxpeak maxpeak])

		saveplot(gcf, [fn_out '_unit_' unitnames{i} '_cross_vel_maxscore_' num2str(max(scoreDFE(i), scoreDRU(i))) '.eps'], 'eps', [6 6]);
%		pause
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Heat map of tau and xcov for all units%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	figure
	gridx = floor(sqrt(nU));
	if (gridx < sqrt(nU))
		gridx = gridx+1;
		gridy = gridx; 
		peakFE = [peakFE, zeros(1, gridx*gridy-nU)];
		peakRU = [peakRU, zeros(1, gridx*gridy-nU)];
		lagFE = [lagFE, zeros(1, gridx*gridy-nU)];
		lagRU = [lagRU, zeros(1, gridx*gridy-nU)];
		scoreFE = [scoreFE, zeros(1, gridx*gridy-nU)];
		scoreRU = [scoreRU, zeros(1, gridx*gridy-nU)];
	end
	subplot(3,2,1)
	image(reshape(abs(peakFE),gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [0 max(abs(peakFE))];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Peak FE cross-correlation'])
	colorbar
	subplot(3,2,2)
	image(reshape(abs(peakRU),gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [0 max(abs(peakRU))];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Peak RU cross-correlation'])
	colorbar
	subplot(3,2,3)
	image(reshape(lagFE,gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [min(lagFE) max(lagFE)];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Lag FE'])
	colorbar
	subplot(3,2,4)
	image(reshape(lagRU,gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [min(lagRU) max(lagRU)];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Lag RU'])
	colorbar
	subplot(3,2,5)
	image(reshape(abs(scoreFE),gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [0 max(abs(scoreFE))];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Score FE cross-correlation' num2str(i)])
	colorbar
	subplot(3,2,6)
	image(reshape(abs(scoreRU),gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [0 max(abs(scoreRU))];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Score RU cross-correlation' num2str(i)])
	colorbar
	saveplot(gcf, [fn_out '_summary.eps'], 'eps', [6 4]);	

	%Plot lag vs scores
	figure
	plot(abs(scoreFE), lagFE, '.b', abs(scoreRU), lagRU, '.r')
	xlabel('|score|')
	ylabel('lag (s)')
	legend('FE', 'RU')
	saveplot(gcf, [fn_out '_lag_vs_score.eps'])

	%Plot lag vs scores
	figure
	plot(abs(scoreDFE), lagDFE, '.b', abs(scoreDRU), lagDRU, '.r')
	xlabel('|score|')
	ylabel('lag (s)')
	legend('FE', 'RU')
	saveplot(gcf, [fn_out '_vel_lag_vs_score.eps'])

	%figure
	%plot(zeros(nU,2), [scoreFE'; scoreRU'], '.b')
	%xlabel('score FE')
	%ylabel('score RU')
	%saveplot(gcf, [fn_out '_popvec.eps'])

end
