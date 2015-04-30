function plot3D_coloct_mean(seq, xspec, octs, varargin)
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
% alph       - transparency of tube plots
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

  dimsToPlot = 1:3;
  nPlotMax   = 20;
  redTrials  = [];
  alph       = 0.4;
  assignopts(who, varargin);

  if size(seq(1).(xspec), 1) < 3
    fprintf('ERROR: Trajectories have less than 3 dimensions.\n');
    return
  end

  nC = size(seq(1).xorth,1);
  nT = min(length(seq), nPlotMax);

  f = figure;
  pos = get(gcf, 'position');
  set(f, 'position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
  
  nPlots = min(length(seq), nPlotMax);

  cm = jet(8);

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

    %dat(1,:) = dat(1,:)-dat(1,1);
    %dat(2,:) = dat(2,:)-dat(2,1); 

  for idx = 1:nC
    dat_mean{idx} = squeeze(mean(dat_trunc{idx},1));
    dat_mean{idx}(1,:) = dat_mean{idx}(1,:)-dat_mean{idx}(1,1);
    dat_mean{idx}(2,:) = dat_mean{idx}(2,:)-dat_mean{idx}(2,1);
    dat_mean{idx}(3,:) = dat_mean{idx}(3,:)-dat_mean{idx}(3,1);
    dat_std{idx} = squeeze(std(dat_trunc{idx},1));
  end

  %Plot
  for n = 1:nC
    lw = 1;
    col = repmat(cm(o,:), Tmin, 1);
    %tubeplot(dat_mean{n}(1,:), dat_mean{n}(2,:), dat_mean{n}(3,:), 0.001, col');
    dat_std_norm = sqrt(dat_mean{n}(1,:).^2+ dat_mean{n}(2,:).^2+ dat_mean{n}(3,:).^2);
    a = cell(4, 1);
    [a{:}] = tubeplot(dat_mean{n}(1,:), dat_mean{n}(2,:), dat_mean{n}(3,:), 0.1*dat_std_norm);
    h = a{4};
    alpha(alph)
    %svgSpecularLightingDistant([h], 'res', 1, 1, 1, 225, 45, 'res');
    %svgSpecularLightingDistant(h, 'blur', 1, 16, 2, 225, 45, 'lighting')
    %svgComposite(s, 'lighting', 'SourceGraphic', 'atop', 'obj');
    %   s : Array of plot object handles
%   source : Any previous defined filter result string, 'SourceGraphic',
%            or 'SourceAlpha'.
%   specularConstant : Specular constant
%   specularExponent : Specular exponent
%   surfaceScale : Surface scaling factor
%   azimuth : Light azimuth angle [deg], typical 225.
%   elevation : Light elevation angle [deg], typical 45.
%   result : String that identifies the filter result for following filter
%            stages.
    drop_shadow_lighting(h, cm(n,:));
    hold on;
  end

  shading interp;
  camlight; lighting gouraud

  %light('Position',[1 0 0],'Style','local')

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
end

function drop_shadow_lighting(s, col)
  set(s, 'EdgeColor', 'none', 'FaceLighting', 'gouraud')
  set(s, 'FaceColor', col)
  % Draws a shadow and adds light effects.
  % PRELIMINARY IMPLEMENTATION (Parameters may change)
  % Parameters:
  %   s : Array of plot object handles
  svgBoundingBox(s, 'axes', 12, 'off');
  %svgGaussianBlur(s, 'SourceAlpha', 2, 'blur');
  svgSpecularLightingDistant(s, 'SourceAlpha', 1, 16, 2, 225, 45, 'lighting')
  svgComposite(s, 'lighting', 'SourceGraphic', 'atop', 'obj');
  %svgGaussianBlur(s, 'SourceAlpha', 5, 'blur2');
  %svgOffset(s, 'blur2', [8 7], 'shade');
  %svgComposite(s, 'obj', 'shade', 'over', 'final');
end