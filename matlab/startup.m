global metaData matpath nevpath
matpath = '/home/lansdell/projects/bci/matlab/labview';
nevpath = '/home/lansdell/projects/bci/matlab/blackrock';
metaData = '/home/lansdell/projects/bci/matlab/labview';

addpath('./functions', './LogProcess', './NPMK', './testdata/');
addpath_recurse('./worksheets');
addpath_recurse('~/matlab/chronux');
%Add extra code

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Arial')
set(0,'DefaultAxesFontSize', 8)

% Change default text fonts.
set(0,'DefaultTextFontname', 'Arial')
set(0,'DefaultTextFontSize', 8)

% Change default line width
set(0,'DefaultLineLineWidth',1.5)
