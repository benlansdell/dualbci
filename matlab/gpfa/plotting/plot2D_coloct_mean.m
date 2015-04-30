function plot2D_coloct_mean(seq, xspec, octs, varargin)
%
% plot2D(seq, xspec, ...)
%
% Plot neural trajectories in a three-dimensional space.
%
% INPUTS:
%
% seq        - data structure containing extracted trajectories
% xspec      - field name of trajectories in 'seq' to be plotted 
%              (e.g., 'xorth' or 'xsm')
%
% OPTIONAL ARGUMENTS:
%
% dimsToPlot - selects three dimensions in seq.(xspec) to plot 
%              (default: 1:3)
% nPlotMax   - maximum number of trials to plot (default: 20)
% redTrials  - vector of trialIds whose trajectories are plotted in red
%              (default: [])
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

  dimsToPlot = 1:2;
  nPlotMax   = 20;
  redTrials  = [];
  assignopts(who, varargin);

  f = figure;
  pos = get(gcf, 'position');
  set(f, 'position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
  
  nPlots = min(length(seq), nPlotMax);
  nMaxBins = 500;

  nC = size(seq(1).xorth,1);
  %nC = 2;
  nT = min(length(seq), nPlotMax);

  %Find min trial length
  for n = 1:min(length(seq), nPlotMax)
    tlen(n) = seq(n).T;
  end
  Tmin = min(tlen);
  for idx = 1:nC
    dat_trunc{idx} = zeros(nT, nC, Tmin);
  end

  %Plot each trial
  for n = 1:nT
    %Projected activity
    dat = seq(n).(xspec);
    target = octs(n);
    dat_trunc{target}(n,:,:) = dat(:,1:Tmin);
  end

  for idx = 1:nC
    dat_mean{idx} = squeeze(mean(dat_trunc{idx},1));
    dat_mean{idx}(1,:) = dat_mean{idx}(1,:)-dat_mean{idx}(1,1);
    dat_mean{idx}(2,:) = dat_mean{idx}(2,:)-dat_mean{idx}(2,1);
    dat_std{idx} = squeeze(std(dat_trunc{idx},1));
  end

  %cm = lines(nPlots);
  cm = jet(8);
  %for n = 1:nC
  for n = 1:2
    lw = 2;
    %dat_std_norm = sqrt(dat_mean{n}(1,:).^2+ dat_mean{n}(2,:).^2);
    %plot3(dat(1,:), dat(2,:), dat(3,:), '.-', 'linewidth', lw, 'color', col);
    %plot(dat(1,:), dat(2,:), '.-', 'linewidth', lw, 'color', cm(o,:));
    %h.Color(4) = 0.5;
    z = zeros(size(dat_mean{n}(1,:)));
    S = surface([dat_mean{n}(1,:);dat_mean{n}(1,:)],[dat_mean{n}(2,:);dat_mean{n}(2,:)],[z;z],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',lw,...
            'edgealpha',.8,...
            'edgecolor',cm(n,:));
    %alpha(alph)
    hold on;
  end
  legend('1', '2', '3', '4', '5', '6', '7', '8');

  axis equal;
  if isequal(xspec, 'xorth')
    str1 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(2));
  else
    str1 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(2));
  end