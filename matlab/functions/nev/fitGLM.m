%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load some training data %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load spike times, and stimulus data for each unit
nevfile = './testdata/20130117SpankyUtah001.nev';
matfile = './testdata/Spanky_2013-01-17-1325.mat';
binsize = 0.002;
nkt = 100;  % Number of time bins in filter;
global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
samplerate = 1/binsize;
RefreshRate = samplerate; 
offset = 0;
threshold = 5;
verbosity = 1;
flip = true; %Flip causality
fn_out = './worksheets/glm_matlab/20130117SpankyUtah001';
[binnedspikes rates torque unitnames tspks] = preprocess_shoham_nev(nevfile, fn_out, binsize, threshold, offset);
nU = length(unitnames);
T = size(binnedspikes, 1)*binsize;

%Load 'stimulus' data
withintrial = zeros(size(rates,1),1);
targetpos = zeros(size(rates,1),2);
startpos = zeros(size(rates,1),2);
rtorque = zeros(size(torque));
%Load trial info
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

  tsp = binnedspikes(:,idx);
  nsp = length(tsp);

  %Make spike history vector
  shist = zeros(nsp, nkt);
  ex_torque = zeros(nsp, nkt*2);
  for j = nkt:nsp
    shist(j,:) = tsp(j-nkt+1:j);
    ex_torque(j,1:nkt) = torque(j-nkt+1:j,1);
    ex_torque(j,nkt+1:end) = torque(j-nkt+1:j,2);
  end

  %Form stim vector
  Stim = [shist, ex_torque];%, ex_rtorque];
  b = glmfit(Stim,tsp,'poisson')
	%Compute STA and use as initial guess for k

	%Fit GLM (traditional version) via max likelihood
	
  %Extract filters fitted...
	
	%Plot results
	figure
  tt = (1:size(sta,1))*binsize*1000;
  if ~flip
  tt = tt - size(sta,1)*binsize*1000;
  end
	subplot(241);  % sta % ------------------------
	plot(tt, sta(:,1));
	title(['unit ' unitnames{idx} '. filter FE']);
	xlabel('time (ms)');
  ylabel('STA');
	
  subplot(242);  % sta % ------------------------
  plot(tt, sta(:,2));
  title(['unit ' unitnames{idx} '. filter RU']);
  xlabel('time (ms)');
  ylabel('STA');

  subplot(243);  % sta % ------------------------
  plot(tt, sta(:,3));
  title(['unit ' unitnames{idx} '. target centric filter FE']);
  xlabel('time (ms)');
  ylabel('STA');

  subplot(244);  % sta % ------------------------
  plot(tt, sta(:,4));
  title(['unit ' unitnames{idx} '. target centric filter RU']);
  xlabel('time (ms)');
  ylabel('STA');



  %subplot(245)
  %%Plot stc
  %title('STC')
  %imagesc(stc0)
  %subplot(246)
  %%Plot singular values of stc
  %for j = 1:size(stc0,1)
  %  svals(j) = s(j,j);
  %end
  %plot(log(svals(1:100)))
  %ylabel('log|singular values|')

  saveplot(gcf, [fn_out '_unit_' unitnames{idx} '_filters.eps'], 'eps', [10,6]);

	%subplot(233); % sta-projection % ---------------
	%imagesc(gg0.k)
	%title('projected STA');
	%
	%subplot(234); % estimated filter % ---------------
	%imagesc(gg1.k) 
	%title('ML estimate: full filter'); xlabel('space'); ylabel('time');
	%
	%subplot(236); % ----------------------------------
	%plot(ggsim.iht,exp(ggsim.ih),'k', gg1.iht,exp(gg1.ihbas*gg1.ih),'b',...
	%    gg2.iht, exp(gg2.ihbas*gg2.ih), 'r');
	%title('post-spike kernel');
	%axis tight;

end

%Once fit model, test with some test data

%How do we compare performance?