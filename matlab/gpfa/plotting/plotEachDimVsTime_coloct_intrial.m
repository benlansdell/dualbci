function plotEachDimVsTime_coloct(seq, xspec, binWidth, octs, intrial, varargin)
%
% plotEachDimVsTime(seq, xspec, binWidth, ...)
%
% Plot each state dimension versus time in a separate panel.
%
% Color traces by octant which trial (target - start) pos lies in 
%
% INPUTS:
%
% seq       - data structure containing extracted trajectories
% xspec     - field name of trajectories in 'seq' to be plotted 
%             (e.g., 'xorth' or 'xsm')
% binWidth  - spike bin width used when fitting model
% octs      - vector indicating which octant trial lies in
%
% OPTIONAL ARGUMENTS:
%
% nPlotMax  - maximum number of trials to plot (default: 20)
% redTrials - vector of trialIds whose trajectories are plotted in red
%             (default: [])
% nCols     - number of subplot columns (default: 4)
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

  nPlotMax  = 20;
  redTrials = [];
  nCols     = 4;
  assignopts(who, varargin);
  nMaxBins = 50000;
  ncol = 8; %number of colors

  f = figure;
  pos = get(gcf, 'position');
  set(f, 'position', [pos(1) pos(2) 2*pos(3) pos(4)]);

  Xall = [seq.(xspec)];
  xMax = ceil(10 * max(abs(Xall(:)))) / 10; % round max value to next highest 1e-1
  
  Tmax    = max([seq.T]);  
  xtkStep = ceil(Tmax/25)*5;
  xtk     = 1:xtkStep:Tmax;
  xtkl    = 0:(xtkStep*binWidth):(Tmax-1)*binWidth;
  ytk     = [-xMax 0 xMax];

  nRows   = ceil(size(Xall, 1) / nCols);
  cm = jet(8);
  cm = [cm; [0.8 0.8 0.8]];
  colormap(cm);


  %Plot each trial
  for n = 1:min(length(seq), nPlotMax)
    %Projected activity
    dat = seq(n).(xspec);
    T   = seq(n).T;
        
    %On each set of axes
    for k = 1:size(dat,1)
      subplot(nRows, nCols, k);
      hold on;
      %Place at same starting point
      tt = ((1:T)-intrial(k,1))*binWidth/1000;
      %dat(k,:) = dat(k,:)-dat(k,intrial(n,1));
      nB = size(dat,2);
      nB = min(nB, nMaxBins);
      dat = dat(:,1:nB);
      lw = 1;
      o = octs(n);
      trialidx = intrial(n,1):(min(nB, intrial(n,2)));
      trialcol = o;
      col = (ncol+1)*ones(T, 1);
      col(trialidx) = trialcol;
  
      size(col);
      z = zeros(size(dat(1,:)));
      S = surface([1:T;1:T],[dat(k,:);dat(k,:)],[z;z],[col'; col'],...%);
              'facecol','no',...
              'edgecol','interp',...
              'linew',lw,...
              'edgealpha',.8);
      %splot(1:T, dat(k,:), 'linewidth', lw, 'color', col);
    end
  end

