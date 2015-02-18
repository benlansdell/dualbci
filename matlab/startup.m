%Charlie's scripts uses these global variables...
%global metaData matpath nevpath
%matpath = '/home/lansdell/projects/bci/matlab/labview';
%nevpath = '/home/lansdell/projects/bci/matlab/blackrock';
%metaData = '/home/lansdell/projects/bci/matlab/labview';

if ismac
	homedir = '/Users/';
else
	homedir = '/home/';
end

addpath([homedir 'lansdell/projects/bci/matlab/functions']);
addpath_recurse([homedir 'lansdell/projects/bci/matlab/functions']);
addpath_recurse([homedir 'lansdell/projects/bci/matlab/preprocess']);
addpath_recurse([homedir 'lansdell/projects/bci/matlab/models']);
addpath_recurse([homedir 'lansdell/projects/bci/matlab/fitting']);
addpath_recurse([homedir 'lansdell/projects/bci/matlab/eval']);
addpath_recurse([homedir 'lansdell/projects/bci/matlab/worksheets']);
addpath_recurse([homedir 'lansdell/projects/bci/matlab/other']);

%nev files, labview files, respectively
addpath([homedir '/lansdell/projects/bci/matlab/blackrock'], [homedir '/lansdell/projects/bci/matlab/labview']);
addpath_recurse([homedir '/matlab/chronux']);
addpath_recurse([homedir '/matlab/arfit']);

%Add extra color
%my_ColorOrder = [   0.00000   0.00000   1.00000;
%   0.00000   0.50000   0.00000;
%   1.00000   0.00000   0.00000;
%   0.00000   0.75000   0.75000;
%   0.75000   0.00000   0.75000;
%   0.75000   0.75000   0.00000;
%   0.25000   0.25000   0.25000;
%   0.00000   0.90000   0.00000;
%   1.00000   0.50000   0.00000;
%   0.95000   0.95000   0.00000;
%   0.60000   0.60000   0.60000;
%   0.00000   0.00000   0.00000;
%   1.00000   0.00000   1.00000];
%set(0,'DefaultAxesColorOrder',my_ColorOrder);

% Change default axes fonts.
%set(0,'DefaultAxesFontName', 'Arial')
%set(0,'DefaultAxesFontSize', 10)

% Change default text fonts.
%set(0,'DefaultTextFontname', 'Arial')
%set(0,'DefaultTextFontSize', 10)

% Change default line width
%set(0,'DefaultLineLineWidth',1.5)
