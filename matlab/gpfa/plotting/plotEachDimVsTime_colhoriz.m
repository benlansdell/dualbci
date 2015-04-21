function plotEachDimVsTime_colhoriz(seq, xspec, binWidth, quads, varargin)
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
% dwelltime - number of bins before trial starts in each trial.spikes matrix. 
%             (default: 0)
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

  nPlotMax  = 124;
  redTrials = [];
  nCols     = 4;
  dwelltime = 0;
  assignopts(who, varargin);

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
  cm = colormap(lines(2));

  %Plot each trial
  for n = 1:min(length(seq), nPlotMax)
    %Projected activity
    dat = seq(n).(xspec);
    T   = seq(n).T;
        
    %On each set of axes
    for k = 1:size(dat,1)
      subplot((nRows+1), nCols, k);
      hold on;
      lw = 0.05;
      if quads(n) == 1 | quads(n) == 4
        col = cm(1,:); 
      else
        col = cm(2,:); 
      end
      nB = min(dwelltime/20, T);
      plot(1:nB, dat(k,1:nB), 'linewidth', lw, 'color', [0.2 0.2 0.2]);
      %Once trial starts color by the side the target is on relative to the start position of the cursor
      if T > nB
        plot((nB+1):T, dat(k,(nB+1):T), 'linewidth', lw, 'color', col);
      end
    end
  end

%Plot axes
for k = 1:size(dat,1)
  h = subplot((nRows+1), nCols, k);
  axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);

  if isequal(xspec, 'xorth')
    str = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$',k);
  else
    str = sprintf('$${\\mathbf x}_{%d,:}$$',k);
  end
  %title(str, 'interpreter', 'latex', 'fontsize', 16);
      
  set(h, 'xtick', xtk, 'xticklabel', xtkl);
  set(h, 'ytick', ytk, 'yticklabel', ytk);
  xlabel('Time (ms)');
end

h = subplot((nRows+1), nCols, size(dat,1)+1)
plot(zeros(2))
axis off 
legend('Quad 1/4','Quad 2/3')
