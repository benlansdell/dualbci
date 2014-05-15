function [scoreFE, scoreRU] = correlation_shoham_nev(nevfile, fn_out, threshold)
	%correlation_shoham_nev	Function to fit the following model to spike and torque data:
	%			lambda_i(t) = \lambda_0 + \sum_j^N k_j^1 . x_i^1(t+\tau+jh) + \sum_j^N k_j^2 . x_i^2(t+\tau+jh)
	%		That is, it fits a linear filter to the torque data. In the above formula, time-step size h
	%		is given by parameter binsize; max value of N is given by kernellength; \tau is given by offset.
	%		Program fits all models with kernels of length between 1 and N, meaning that for each single-unit N models
	%		will be fit.
	%
	%		Usage:
	%			[scoreFE, scoreRU] = correlation_shoham_nev(nevfile, fn_out, threshold)
	%
	% 	Input:
	%			nevfile = file to process
	%			fn_out = base name for output plots
	%			threshold = (optional, default = 5) threshold firing rate below which unit is ignored
	%		
	%		Output:
	%			(none) produces plots of the filter for each single-unit channel, along with summary of quality of fits for 
	%			each channel
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
	%			threshold = 1;
	%			fn_out = './worksheets/tuning/crosscorr/20130117SpankyUtah001';
	%			correlation_shoham_nev(nevfile, fn_out, threshold);
	
	%Optional arguments
	if (nargin < 3)	threshold = 1; end
	%Set this to above 0 if want debug info/plots
	verbosity = 1;
	%Offset to apply
	%offset = 0;
	%Put firing rate ahead of torque
	offset = -0.125;
  sigma_fr = 0;
  binsize = 0.001;
	%Size of gaussian filter to apply
	samplerate = 1/binsize;
	%Max lag for correlations
	maxlag = 90;
	maxlag = round(maxlag*samplerate)/samplerate;
	maxpeak = 3;
  %Preprocess spike and torque data
  [binnedspikes rates torque dtorque ddtorque unitnames] = preprocess_shoham_nev(nevfile, fn_out, binsize, sigma_fr, threshold, offset, verbosity);
  nU = length(unitnames);

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
	xlabel('time (s)')
   	title('Auto-corr torque FE');
	subplot(1,2,2);
   	plot(tt,autotorqueRU);
   	xlim([-maxpeak*2 maxpeak*2])
	xlabel('time (s)')
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
    	scoreDRU(i) = peakDRU(i)/stdDRU(i);
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

        %Plot lag vel vs lag
        figure
        plot(lagFE, lagDFE, '.b', lagRU, lagDRU, '.r')
        xlabel('lag')
        ylabel('lag vel')
        legend('FE', 'RU')
        saveplot(gcf, [fn_out '_lags.eps'])

        %Plot score vel vs score
        figure
        plot(abs(scoreFE), abs(scoreDFE), '.b', abs(scoreRU), abs(scoreDRU), '.r')
        xlabel('|score|')
        ylabel('|score| (vel)')
        legend('FE', 'RU')
        saveplot(gcf, [fn_out '_scores.eps'])

	%Threshold score above which to look at mean lags
	thr = 4;
	meanlFE = mean(lagFE);
	meanlRU = mean(lagRU);
	meanlDFE = mean(lagDFE);
	meanlDRU = mean(lagDRU);
	meanlFEthr = mean(lagFE(scoreFE>thr));
	meanlRUthr = mean(lagRU(scoreRU>thr));
	meanlDFEthr = mean(lagDFE(scoreDFE>thr));
	meanlDRUthr = mean(lagDRU(scoreDRU>thr));
	display(['mean lag FE: ' num2str(meanlFE)])
	display(['mean lag RU: ' num2str(meanlRU)])
	display(['mean lag vel FE: ' num2str(meanlDFE)])
	display(['mean lag vel RU: ' num2str(meanlDRU)])
	display(['mean lag FE (score > 4): ' num2str(meanlFEthr)])
	display(['mean lag RU (score > 4): ' num2str(meanlRUthr)])
	display(['mean lag vel FE (score > 4): ' num2str(meanlDFEthr)])
	display(['mean lag vel RU (score > 4): ' num2str(meanlDRUthr)])

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

end
