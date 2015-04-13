function plotEachDimVsTime_colquad(seq, xspec, binWidth, quads, varargin)
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

  nPlotMax  = 124;
  redTrials = [];
  nCols     = 4;
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
  cm = colormap(lines(4));

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
      col = cm(quads(n),:); 
      plot(1:T, dat(k,:), 'linewidth', lw, 'color', col);
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
plot(zeros(4))
axis off 
legend('Quad 1','Quad 2','Quad 3','Quad 4')
