function [b, dev, stats] = fitGLM_LS(nevfile, matfile, fn_out, binsize, nkt, threshold, units)
  %   [b, dev, stats] = fitGLM_LS(nevfile, matfile, fn_out, binsize, threshold)
  %
  %     Fit GLM to spike data from blackrock recording file for each unit above a specified threshold
  %     
  %
  %     Input:
  %     nevfile = Blackrock recording file to process
  %     matfile = Labview trial data file containing corresponding trial information (for intentionality estimation)
  %     fn_out = file to output diagnostic plots to
  %     binsize = (optional, default = 0.002) size of timebins over which to compute regression
  %     nkt = (optional, default = 100) size of filters to fit
  %     threshold = (optional, default = 5) threshold firing rate below which unit is ignored
  %     units = (optional, default = None) cell array of strings. if specified the only compute for units matching the names provided. ie: '55.1'
  %   
  %   Output:
  %     b = [nB x nU] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
  %     dev = deviance
  %     stats = fitting statistics output from glmfit, or from lassoglm
  %  
  %   Test code:
  %     nevfile = './testdata/20130117SpankyUtah001.nev';
  %     matfile = './testdata/Spanky_2013-01-17-1325.mat';
  %     binsize = 0.002;
  %     nkt = 100;
  %     threshold = 5;
  %     fn_out = './worksheets/glm_ls/plots/20130117SpankyUtah001';
  %     units = {'55.1'};
  %     [b, dev, stats] = fitGLM_lasso(nevfile, matfile, fn_out, binsize, nkt, threshold, units);

  assert(nargin >= 3, 'Need more than 2 input arguments.');
  if (nargin < 4) binsize = 0.002; end
  if (nargin < 5) nkt = 100; end
  if (nargin < 6) threshold = 5; end
  if (nargin == 7) unitsubset = 1; end 
  global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
  samplerate = 1/binsize;
  offset = 0;
  coarsebs = 0.05;
  %Make a Gaussian filter to smooth
  coarseSR = 1/coarsebs; %Sample rate for which to smooth values (1/50ms = 20Hz)
  sigma_fr = 0.25;
  sz = 30;
  sigma_fr = sigma_fr*coarseSR;
  sz = sigma_fr*3*2;
  x = linspace(-sz/2, sz/2, sz);
  gaussFilter = exp(-x.^2/(2*sigma_fr^2));
  gaussFilter = gaussFilter/sum(gaussFilter);

  %Load spike times, and stimulus data for each unit
  [binnedspikes rates torque dt ddt unitnames tspks] = preprocess_shoham_nev(nevfile, fn_out, binsize, threshold, offset);
  %Bin spikes more coarsely for comparison later...
  %[binnedspikes_coarse, rates_coarse] = preprocess_shoham_nev(nevfile, fn_out, coarsebs, threshold, offset);
  nU = length(unitnames);
  T = size(binnedspikes, 1)*binsize;

  verbosity = 0;
  withintrial = zeros(size(rates,1),1);
  targetpos = zeros(size(rates,1),2);
  startpos = zeros(size(rates,1),2);
  rtorque = zeros(size(torque));
  %Load trial data
  trials = import_trials(matfile);
  for idx=1:length(trials)
    	%For each trial, use only those corresponding to the given nev file
    	[path1, name1, ext1] = fileparts(nevfile);
    	[path2, name2, ext2] = fileparts(trials(idx).nevfile);
    	if strcmp(name1, name2)
    	  	%mark the within trial times in the above vector
    	  	%Time relative to nev file recording
    	  	trialstart = round(samplerate*(trials(idx).starttime-trials(idx).offset));
    	  	trialend = round(samplerate*(trials(idx).endtime-trials(idx).offset));
    	  	l = length(trialstart:trialend);
    	  	withintrial(trialstart:trialend)=1;
    	  	targetpos(trialstart:trialend,:)=repmat(trials(idx).target, l, 1);
    	  	startpos(trialstart:trialend,:)=repmat(trials(idx).cursorstart, l, 1);
    	  	%Compute torque to subtract from torque to make target position at origin
    	  	trqend = torque(trialend,:);
    	  	%Center torque
    	  	rtorque(trialstart:trialend,:) = torque(trialstart:trialend,:)-repmat(trqend, l, 1);
    	  	trqstr = rtorque(trialstart,:);
    	  	%Compute angle to rotate torque, dtorque and ddtorque by so that cursor start is at directly below target
    	  	theta = atan(trqstr(2)/trqstr(1));
    	  	if trqstr(1)<0
    	  		if trqstr(2)>0
    	  			theta = theta + pi;
    	  		else
    	  			theta = theta - pi;
    	  		end
    	  	end
    	  	%Plot to make sure things make sense
    	  	%Before rotation
    	  	if verbosity > 1
  	  	  	figure
    	  		subplot(1,2,1)
    	  		plot(rtorque(trialstart:trialend,1), rtorque(trialstart:trialend,2), trqstr(1), trqstr(2), 'or');
    	  		title('Cursor position before rotation. red = start')
    	  		xlim([-0.5 0.5])
    	  		ylim([-0.5 0.5])
    	  	end
    	  	%Place directly below torque end position (theta = pi/2)
    	  	%theta = -theta - pi/2;
    	  	%OR don't do any rotation
    	  	theta = 0;
    	  	%Rotate torque
    	  	rotation = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    	  	rtorque(trialstart:trialend,:) = (rotation*(rtorque(trialstart:trialend,:)'))';
    	  	trqstr = (rotation*trqstr')';
    	  	%Plot after rotation
    	  	if verbosity > 1
    	  		subplot(1,2,2)
    	  		plot(torque(trialstart:trialend,1), rtorque(trialstart:trialend,2), trqstr(1), trqstr(2), 'or');
    	  		title('Cursor position after rotation. red = start')
    	  		xlim([-0.5 0.5])
    	  		ylim([-0.5 0.5])
    	    end
    	  	display(['Trial ' num2str(idx) ' within nev file. t_start: ' num2str(trials(idx).starttime) ' t_end: ' num2str(trials(idx).endtime)]);
    	  	%pause
    	end
  end
  
  %Truncate everything so that only data within trial is included
  %only if using rtorque is this needed...
  
  %For each unit, fit a GLM to the torque data
  for idx=1:nU 
    if unitsubset & ~ismember(units, unitnames{idx})
      continue
    end
    display(['Fitting GLM for unit ', unitnames{idx}])
    %Binned spikes for given unit
    tsp = binnedspikes(:,idx);
    sps = find(tsp);
    lsp = length(tsp);
    display(['using ' num2str(sum(tsp)) ' spikes.'])
  
    %Make stimulus vector
    shist = zeros(lsp, nkt);
    ex_torque = zeros(lsp, nkt*2);
    for j = (nkt+1):(lsp-nkt)
      shist(j,:) = tsp(j-nkt:j-1);
      %Use torque
      ex_torque(j,1:nkt) = torque(j:(j+nkt-1),1);
      ex_torque(j,nkt+1:end) = torque(j:(j+nkt-1),2);
      %Or use velocity
      %ex_torque(j,1:nkt) = dt(j:(j+nkt-1),1);
      %ex_torque(j,nkt+1:end) = dt(j:(j+nkt-1),2);
      %Or use acceleration
      %ex_torque(j,1:nkt) = ddt(j:(j+nkt-1),1);
      %ex_torque(j,nkt+1:end) = ddt(j:(j+nkt-1),2);
    end
  
    %Form stim vector
    Stim = [shist ex_torque];%, ex_rtorque];
    %Truncate to exclude start and end of recording
    Stim = Stim(nkt+1:end-nkt,:);
    tsp = tsp(nkt+1:end-nkt);
    %Only focus on stimulus at spike times...
    %Stim = Stim(sps);

    %If only running tests, then take a subsection of data
    tStim = double(Stim);
    ttsp = double(tsp);
    %Add a column of ones for the constant term
    tStim = [ones(size(tStim, 1),1), tStim];
    tStim = tStim(1:min(50000, size(tStim, 1)),:);
    ttsp = tsp(1:min(50000, size(Stim, 1)));
    %[b, dev, stats] = glmfit(tStim,ttsp,'binomial');
    %[b, dev, stats] = glmfit(tStim,ttsp,'poisson');
    [b, fitinfo] = lassoglm(tStim, ttsp, 'poisson', 'Alpha',.25,'CV',3,'NumLambda',20,'LambdaRatio',10^-4)

    selected_lambda=fitinfo.IndexMinDeviance % index of best model based on min deviance
    % selected_lambda=1 % alternatively just pick a model by hand
    bb=b(:,selected_lambda); % bb is [k1; k2; h]
    bbb=[fitinfo.Intercept(selected_lambda); bb]; % bbb is [b; k1, k2; h] This is the vector with all of fit coefficients

    %Extract filters fitted...
    const = bbb(1);
    fil_shist = bbb(2:nkt+1);
    fil_torque1 = bbb(nkt+2:2*nkt+1);
    fil_torque2 = bbb(2*nkt+2:3*nkt+1);
    
  	%Plot results
  	figure
    tt = (1:nkt)*binsize*1000;
  	subplot(241);  % sta % ------------------------
  	plot(tt, fil_shist);
  	title(['unit ' unitnames{idx} '. spike history filter']);
  	xlabel('time (ms)');
  	
    subplot(242);  % sta % ------------------------
    plot(tt, fil_torque1);
    title(['unit ' unitnames{idx} '. RU filter']);
    xlabel('time (ms)');
  
    subplot(243);  % sta % ------------------------
    plot(tt, fil_torque2);
    title(['unit ' unitnames{idx} '. FE filter']);
    xlabel('time (ms)');

    %subplot(244); %Prediction using model
    %Smooth actual spikes and predicted model
    %smthLNfittedrates(i,:) = conv(LNfittedrates(i,:), gaussFilter, 'same');
    %smthfittedrates(i,:) = conv(fittedrates(i,:), gaussFilter, 'same');
    %smthrates(i,:) = conv(squeeze(rates(:,i)), gaussFilter, 'same');

    t_i = 60;
    t_f = 66; %three seconds worth of data
    tt = t_i:binsize:t_f;
    ii = round(tt/binsize)+1;
    rho_hat = glmval(bbb, tStim, 'log'); % Uses the fit coefficients and the original input data to generate the ouput rho
    %Resample for 50ms time bins
    %rho_hat_RS = resample(rho_hat, coarseSR, samplerate);
    %rho_hat_RS = rebin(rho_hat, binsize, coarsebs);
    %hold on
    %plot(tt, rho_hat(ii)*20);
    %%Plot spikes per time bin, also
    %%plot(sps, zeros(size(sps)), 'r.');
    %%Plot smoothed firing rate
    %%smoothedrates = conv(rates_coarse(:,idx), gaussFilter, 'same');
    %%plot(ttt, smoothedrates(iii));
    %plot(ttt, binnedspikes_coarse(iii,idx), 'r');
    %xlim([min(ttt) max(ttt)])
    %xlabel('time (s)')
    %ylabel('estimated firing rate (Hz)')

    subplot(2,4,[5 6 7 8]); %Try the same spike plot but another way...
    hold on
    %Place a bar every time there's a spike
    for j=1:length(tspks(idx).times)
    %for j=1:length(sps)
      %Convert to seconds
      %sp = sps(j)*binsize;
      %Or just use the spike times computed above
      sp = tspks(idx).times(j);
      if (sp < max(tt) & sp > min(tt))
        plot([sp sp], [0 1.2*max(rho_hat(ii))], 'Color', [0.66 0.66 0.66]);
      end
    end
    plotyy(tt, rho_hat(ii), tt, [torque(ii,1), torque(ii,2)]);
    xlim([min(tt) max(tt)])
    ylim([0 1.2*max(rho_hat(ii))])
    xlabel('time (s)')
    %ylabel('estimated firing rate (Hz)')

    saveplot(gcf, [fn_out '_unit_' unitnames{idx} '_filters.eps'], 'eps', [10,5]);
  
  end
