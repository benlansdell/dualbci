function plot2D_coloct(seq, xspec, octs, varargin)
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

  %cm = lines(nPlots);
  cm = jet(8);
  for n = 1:min(length(seq), nPlotMax)
    dat = seq(n).(xspec)(dimsToPlot,:);
    %Place at same starting point
    %dat(1,:) = dat(1,:)-dat(1,1);
    %dat(2,:) = dat(2,:)-dat(2,1);    
    nB = size(dat,2);
    nB = min(nB, nMaxBins);
    dat = dat(:,1:nB);
    T   = seq(n).T;
        
    %if ismember(seq(n).trialId, redTrials)
    %  col = [1 0 0]; % red
    %  lw  = 3;
    %else
    %  col = 0.2 * [1 1 1]; % gray
    %  lw = 0.5;
    %end
    lw = 1;
    o = octs(n);
    %plot3(dat(1,:), dat(2,:), dat(3,:), '.-', 'linewidth', lw, 'color', col);
    %plot(dat(1,:), dat(2,:), '.-', 'linewidth', lw, 'color', cm(o,:));
    %h.Color(4) = 0.5;
    z = zeros(size(dat(1,:)));
    S = surface([dat(1,:);dat(1,:)],[dat(2,:);dat(2,:)],[z;z],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',lw,...
            'edgealpha',.8,...
            'edgecolor',cm(o,:));
    hold on;
  end

  axis equal;
  if isequal(xspec, 'xorth')
    str1 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(2));
    str3 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(3));
  else
    str1 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(2));
    str3 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(3));
  end
  %xlabel(str1, 'interpreter', 'latex', 'fontsize', 24);
  %ylabel(str2, 'interpreter', 'latex', 'fontsize', 24);
  %zlabel(str3, 'interpreter', 'latex', 'fontsize', 24);
