function plot3D_coloct_intrial(seq, xspec, octs, intrial, varargin)
%
% plot3D(seq, xspec, ...)
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

  dimsToPlot = 1:3;
  nPlotMax   = 20;
  redTrials  = [];
  assignopts(who, varargin);
  nMaxBins = 50000;
  ncol = 8; %number of colors

  if size(seq(1).(xspec), 1) < 3
    fprintf('ERROR: Trajectories have less than 3 dimensions.\n');
    return
  end

  f = figure;
  pos = get(gcf, 'position');
  set(f, 'position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
  
  nPlots = min(length(seq), nPlotMax);

  %cm = lines(nPlots);
  cm = jet(8);
  cm = [cm; [0.8 0.8 0.8]];
  colormap(cm);

  for n = 1:min(length(seq), nPlotMax)
    dat = seq(n).(xspec)(dimsToPlot,:);
    T   = seq(n).T;
        
    %Place at same starting point
    dat(1,:) = dat(1,:)-dat(1,intrial(n,1));
    dat(2,:) = dat(2,:)-dat(2,intrial(n,1));    
    dat(3,:) = dat(3,:)-dat(3,intrial(n,1));    
    nB = size(dat,2);
    nB = min(nB, nMaxBins);
    dat = dat(:,1:nB);
    T   = seq(n).T;        
    lw = 1.5;
    o = octs(n);
    trialidx = intrial(n,1):(min(nB, intrial(n,2)));
    trialcol = o;
    col = (ncol+1)*ones(T, 1);
    col(trialidx) = trialcol;

    size(col);
    S = mesh([dat(1,:);dat(1,:)],[dat(2,:);dat(2,:)],[dat(3,:);dat(3,:)],[col'; col'],...%);
            'facecol','no',...
            'edgecol','interp',...
            'linew',lw,...
            'edgealpha',.6);
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
