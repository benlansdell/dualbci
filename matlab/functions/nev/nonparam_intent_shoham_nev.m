function nonparam_intent_shoham_nev(nevfile, matfile, fn_out, threshold)
	%nonparam_intent_shoham_nev	Function to estimate the spike density
  %
  %         p(spike| x,y) ~ p (x,y|spike) / p(x,y)
  %
  %   with cursor coordinates always centered around target position. Note that translating position 
  %   won't affect velocity or accel encoding.
	%
	%		Usage:
	%			nonparam_intent_shoham_nev(nevfile, fn_out, threshold)
	%
	% 	Input:
	%			nevfile = file to process
  %     matfile = .mat file containing trial information from labview
	%			fn_out = base name for output plots
	%			threshold = (optional, default = 5) threshold firing rate below which unit is ignored
	%		
	%		Output:
	%			(none) produces plots of the filter for each single-unit channel, along with summary of quality of fits for 
	%			each channel
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
  %     matfile = './testdata/Spanky_2013-01-17-1325.mat';
	%			threshold = 5;
	%			fn_out = './worksheets/glm/july';
	%			nonparam_intent_shoham_nev(nevfile, matfile, fn_out, threshold);
	
	%Optional arguments
	if (nargin < 3)	threshold = 5; end
	%Offset to apply
	%Put firing rate ahead of torque
	offset = -0.125;
  binsize = 0.001;
	%Size of gaussian filter to apply
	samplerate = 1/binsize;
	%Max lag for correlations
	maxlag = 90;
	maxlag = round(maxlag*samplerate)/samplerate;
	maxpeak = 3;
  %Preprocess spike and torque data
  [binnedspikes rates torque dtorque ddtorque unitnames] = preprocess_shoham_nev(nevfile, fn_out, binsize, threshold, offset);
  nU = length(unitnames);

    %Trim to data occurring within a trial
    %Create a vector the length of the number of bins used here
    withintrial = zeros(size(rates,1),1);
    targetpos = zeros(size(rates,1),2);
    startpos = zeros(size(rates,1),2);
    %Load trial info
    trials = import_trials(matfile);
    for idx=1:length(trials)
        %Use only those within the corresponding nev file
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
            torque(trialstart:trialend,:) = torque(trialstart:trialend,:)-repmat(trqend, l, 1);
            trqstr = torque(trialstart,:);
            %Compute angle to rotate torque, dtorque and ddtorque by so that cursor start is at directly below target
            theta = atan(trqstr(2)/trqstr(1));
            if trqstr(1)<0
              if trqstr(2)>0
                theta = theta + pi;
              else
                theta = theta - pi;
              end
            end
            %Before rotation
            %Place directly below torque end position (theta = pi/2)
            theta = -theta - pi/2;
            %OR don't do any rotation
            %theta = 0;
            %Rotate torque, dtorque and ddtorque
            rotation = [cos(theta), -sin(theta); sin(theta), cos(theta)];
            torque(trialstart:trialend,:) = (rotation*(torque(trialstart:trialend,:)'))';
            dtorque(trialstart:trialend,:) = (rotation*(dtorque(trialstart:trialend,:)'))';
            ddtorque(trialstart:trialend,:) = (rotation*(ddtorque(trialstart:trialend,:)'))';
            trqstr = (rotation*trqstr')';
            display(['Trial ' num2str(idx) ' within nev file. t_start: ' num2str(trials(idx).starttime) ' t_end: ' num2str(trials(idx).endtime)]);
            %pause
        end
    end
    %Now have a vector of times that occur within a trial, trim rates and torque accordingly
    rates = rates(withintrial==1,:);
    torque = torque(withintrial==1,:);
    dtorque = dtorque(withintrial==1,:);
    ddtorque = ddtorque(withintrial==1,:);
    binnedspikes = binnedspikes(withintrial==1,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute p(x,y|spike)/p(x,y) ~ p(spike|x,y)%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	dtorque(1:5,:) = 0;
	ddtorque(1:5,:) = 0;
	figure 

  nx = 64;
  %Compute p(x,y) and p(x,y|spike) histograms for (x,y), (v_x,v_y) and (a_x,a_y)
  histxy = zeros(nx);
  histvxy = zeros(nx);
  histaxy = zeros(nx);

  meansxy = zeros(nU, 2);
  meansvxy = zeros(nU, 2);
  meansaxy = zeros(nU, 2);

  maxxy = (std(torque(:,1))+std(torque(:,2)));
  maxvxy = (std(dtorque(:,1))+std(dtorque(:,2)))/2;
  maxaxy = (std(ddtorque(:,1))+std(ddtorque(:,2)))/2;

  win = fspecial('gaussian',21,2);
  win = win ./ sum(sum(win));

  for idx=1:size(torque,1)

    x = torque(idx,1); y = torque(idx,2);
    vx = dtorque(idx,1); vy = dtorque(idx,2);
    ax = ddtorque(idx,1); ay = ddtorque(idx,2);

    x = floor(nx*(x/maxxy+1)/2)+1;
    y = floor(nx*(y/maxxy+1)/2)+1;
    vx = floor(nx*(vx/maxvxy+1)/2)+1;
    vy = floor(nx*(vy/maxvxy+1)/2)+1;
    ax = floor(nx*(ax/maxaxy+1)/2)+1;
    ay = floor(nx*(ay/maxaxy+1)/2)+1;

    if x >= 1 & x <= nx & y >= 1 & y <= nx
      histxy(x,y) = histxy(x,y) + 1;
    end
    if vx >= 1 & vx <= nx & vy >= 1 & vy <= nx
      histvxy(vx,vy) = histvxy(vx,vy) + 1;
    end 
    if ax >= 1 & ax <= nx & ay >= 1 & ay <= nx
      histaxy(ax,ay) = histaxy(ax,ay) + 1;
    end

  end

  histxy = (histxy+1)/sum(sum(histxy));
  histvxy = (histvxy+1)/sum(sum(histvxy));
  histaxy = (histaxy+1)/sum(sum(histaxy));
  zaxis = [0 1];

	subplot(2,2,1)
  %Heat map of position density
  image(histxy, 'CDataMapping', 'scaled');
  caxis(zaxis);
  set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
  colorbar;
  title('Torque position density')
	subplot(2,2,2)
  %Heat map of velocity density
  image(histvxy, 'CDataMapping', 'scaled');
  caxis(zaxis);
  set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
  colorbar;
  title('Torque velocity density')
	subplot(2,2,3)
  %Heat map of accel density
  image(histaxy, 'CDataMapping', 'scaled');
  caxis(zaxis);
  set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
  colorbar;
  title('Torque accel density')
  saveplot(gcf, [fn_out '_torque_density.eps'], 'eps', [6 6]);
  for i=1:nU
    histxy_sp = zeros(nx);
    histvxy_sp = zeros(nx);
    histaxy_sp = zeros(nx);
  	spikeidx = find(binnedspikes(:,i)>0);

    for idx=1:length(spikeidx)

      x = torque(spikeidx(idx),1); y = torque(spikeidx(idx),2);
      vx = dtorque(spikeidx(idx),1); vy = dtorque(spikeidx(idx),2);
      ax = ddtorque(spikeidx(idx),1); ay = ddtorque(spikeidx(idx),2);

      x = floor(nx*(x/maxxy+1)/2)+1;
      y = floor(nx*(y/maxxy+1)/2)+1;
      vx = floor(nx*(vx/maxvxy+1)/2)+1;
      vy = floor(nx*(vy/maxvxy+1)/2)+1;
      ax = floor(nx*(ax/maxaxy+1)/2)+1;
      ay = floor(nx*(ay/maxaxy+1)/2)+1;

      if x >= 1 & x <= nx & y >= 1 & y <= nx
        histxy_sp(x,y) = histxy_sp(x,y) + 1;
      end
      if vx >= 1 & vx <= nx & vy >= 1 & vy <= nx
        histvxy_sp(vx,vy) = histvxy_sp(vx,vy) + 1;
      end 
      if ax >= 1 & ax <= nx & ay >= 1 & ay <= nx
        histaxy_sp(ax,ay) = histaxy_sp(ax,ay) + 1;
      end
    end

    histxy_sp = histxy_sp/sum(sum(histxy_sp));
    histvxy_sp = histvxy_sp/sum(sum(histvxy_sp));
    histaxy_sp = histaxy_sp/sum(sum(histaxy_sp));

  	figure
  	%Position
  	subplot(3,2,1)
    image(histxy_sp./(histxy), 'CDataMapping', 'scaled');
    colorbar;
  	title(unitnames{i})
    subplot(3,2,2)
    image(imfilter(histxy_sp./(histxy),win), 'CDataMapping', 'scaled');
    colorbar;
    title(unitnames{i})

  	%Velocity
  	subplot(3,2,3)
    image(histvxy_sp./(histvxy), 'CDataMapping', 'scaled');
    colorbar;
    title(unitnames{i})
    subplot(3,2,4)
    image(imfilter(histvxy_sp./(histvxy),win), 'CDataMapping', 'scaled');
    colorbar;
    title(unitnames{i})

  	%Accel
  	subplot(3,2,5)
    image(histaxy_sp./(histaxy), 'CDataMapping', 'scaled');
    colorbar;
    title(unitnames{i})
    subplot(3,2,6)
    image(imfilter(histaxy_sp./(histaxy),win), 'CDataMapping', 'scaled');
    colorbar;
    title(unitnames{i})

  	saveplot(gcf, [fn_out '_spikedensity_unit_' unitnames{i} '.eps'], 'eps', [5 5]);
  	%pause

    %For all of them plot a set of arrows for each one's mean 
    a = histxy_sp./(histxy);
    a = a / sum(sum(a));
    b = histvxy_sp./(histvxy);
    b = b / sum(sum(b));
    c = histaxy_sp./(histaxy);
    c = c / sum(sum(c));
    for j = 1:nx
      for k = 1:nx
        meansxy(i, 1) = meansxy(i, 1) + j*a(j,k);
        meansxy(i, 2) = meansxy(i, 2) + k*a(j,k);
        meansvxy(i, 1) = meansvxy(i, 1) + j*b(j,k);
        meansvxy(i, 2) = meansvxy(i, 2) + k*b(j,k);
        meansaxy(i, 1) = meansaxy(i, 1) + j*c(j,k);
        meansaxy(i, 2) = meansaxy(i, 2) + k*c(j,k);
      end
    end
  end

  %Compute mean
  meansxy = meansxy/nx - 0.5;
  meansvxy = meansvxy/nx - 0.5;
  meansaxy = meansaxy/nx - 0.5;

  %Plot arrows
  figure
  for idx = 1:nU
    subplot(131)
    hold on
    plot([0, meansxy(idx, 1)], [0, meansxy(idx, 2)])
    subplot(132)
    hold on
    plot([0, meansvxy(idx, 1)], [0, meansvxy(idx, 2)])
    subplot(133)
    hold on
    plot([0, meansaxy(idx, 1)], [0, meansaxy(idx, 2)])
  end

  saveplot(gcf, [fn_out '_mean_arrows.eps'], 'eps', [6 2])

end
