function [scoreFE, scoreRU] = nonparam_shoahm_nev(nevfile, fn_out, threshold)
	%nonparam_shoahm_nev	Function to fit the following model to spike and torque data:
	%			lambda_i(t) = \lambda_0 + \sum_j^N k_j^1 . x_i^1(t+\tau+jh) + \sum_j^N k_j^2 . x_i^2(t+\tau+jh)
	%		That is, it fits a linear filter to the torque data. In the above formula, time-step size h
	%		is given by parameter binsize; max value of N is given by kernellength; \tau is given by offset.
	%		Program fits all models with kernels of length between 1 and N, meaning that for each single-unit N models
	%		will be fit.
	%
	%		Usage:
	%			[scoreFE, scoreRU] = nonparam_shoahm_nev(nevfile, fn_out, threshold)
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
	%			threshold = 5;
	%			fn_out = './worksheets/tuning/crosscorr/20130117SpankyUtah001';
	%			nonparam_shoahm_nev(nevfile, fn_out, threshold);
	
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

  maxxy = (std(torque(:,1))+std(torque(:,2)));
  maxvxy = (std(dtorque(:,1))+std(dtorque(:,2)));
  maxaxy = (std(ddtorque(:,1))+std(ddtorque(:,2)));

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

  histxy = (histxy+1)/max(max(histxy));
  histvxy = (histvxy+1)/max(max(histvxy));
  histaxy = (histaxy+1)/max(max(histaxy));
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
  saveplot(gcf, [fn_out '_torque_density.eps'], 'eps', [2 6]);
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

    histxy_sp = histxy_sp/max(max(histxy_sp));
    histvxy_sp = histvxy_sp/max(max(histvxy_sp));
    histaxy_sp = histaxy_sp/max(max(histaxy_sp));

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
  end
end
