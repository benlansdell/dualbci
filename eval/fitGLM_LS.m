function [b, dev, stats] = fitGLM_LS(nevfile, matfile, fn_out, binsize, nkt, threshold, units)
  %   [b, dev, stats] = fitGLM_LS(nevfile, matfile, fn_out, binsize, threshold)
  %
  %     Fit GLM to spike data from blackrock recording file for each unit above a specified threshold
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
  %     nkt = 200;
  %     threshold = 5;
  %     fn_out = './worksheets/glm_ls/plots/20130117SpankyUtah001';
  %     units = {'55.1'};
  %     [b, dev, stats] = fitGLM(nevfile, matfile, fn_out, binsize, nkt, threshold, units);

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
  sigma_fr = 0.5;
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
    %Truncate other data for comparison, too
    binnedspikes = binnedspikes(nkt+1:end-nkt,:);
    rates = rates(floor(nkt/2)+1:end-nkt,:);
    torque = torque(nkt+1:end-nkt,:); 
    dt = dt(nkt+1:end-nkt,:);
    ddt = dt(nkt+1:end-nkt,:);

    %If only running tests, then take a subsection of data
    tStim = double(Stim);
    ttsp = double(tsp);
    tStim = Stim(1:min(50000, size(Stim, 1)),:);
    ttsp = tsp(1:min(50000, size(Stim, 1)));
    [b, dev, stats] = glmfit(tStim,ttsp,'poisson');
   
    %Extract filters fitted...
    const = b(1);
    fil_shist = b(2:nkt+1);
    fil_torque1 = b(nkt+2:2*nkt+1);
    fil_torque2 = b(2*nkt+2:3*nkt+1);
    
    %Or use lassoglm and compare...
    %[b, fitinfo] = lassoglm(tStim, ttsp, 'binomial', 'Alpha',.25,'CV',3,'NumLambda',20,'LambdaRatio',10^-4)
  	%fil_shist = b(1:nkt);
    %fil_torque1 = b(nkt+1:2*nkt);
    %fil_torque2 = b(2*nkt+1:3*nkt);
    
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

    t_i = 30;
    t_f = 40; %three seconds worth of data
    tt = t_i:binsize:t_f;
    ii = round(tt/binsize)+1;
    rho_hat = glmval(b, tStim, 'log'); % Uses the fit coefficients and the original input data to generate the ouput rho

    subplot(2,4,[5 6 7 8]); %Try the same spike plot but another way...
    hold on

    %smthLNfittedrates(i,:) = conv(LNfittedrates(i,:), gaussFilter, 'same');
    smthfittedrates = conv(rho_hat, gaussFilter, 'same')/binsize;
    smthrates = conv(squeeze(rates(:,idx)), gaussFilter, 'same');
    smthrates2 = conv(ttsp, gaussFilter, 'same')/binsize;
    plot(tt, smthrates2(ii), tt, smthfittedrates(ii))
    %plot(tt, smthrates(ii), tt, smthfittedrates(ii))
    legend('Actual', 'GLM')
    xlim([t_i, t_f])
    xlabel('time (s)')
    ylabel('estimated firing rate (Hz)')

    saveplot(gcf, [fn_out '_unit_' unitnames{idx} '_filters.eps'], 'eps', [10,5]);
  
  end
