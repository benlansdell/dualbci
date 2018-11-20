if ismac
	homedir = '/Users/';
else
	homedir = '/home/';
end

%Change this to local database
databaseuser = 'root';
databasepwd = 'Eskimo1!';
databaseurl = 'jdbc:mysql://localhost:3306/spanky_db';

%Add database jar file...
javaaddpath([homedir 'lansdell/matlab/mysql-connector-java-5.1.35-bin.jar']);
%javaaddpath([homedir 'lansdell/matlab/mysql-connector-java-8.0.11.jar']);

localdir = [homedir 'lansdell/projects/dualbci'];
matlabdir = [homedir 'lansdell/matlab'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add paths to other packages%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath([localdir '/functions']);
addpath_recurse([localdir '/functions']);
addpath_recurse([localdir '/preprocess']);
addpath_recurse([localdir '/models']);
addpath_recurse([localdir '/fitting']);
addpath_recurse([localdir '/eval']);
addpath_recurse([localdir '/scripts']);
addpath_recurse([localdir '/sql']);

%nev files, labview files, respectively
addpath([localdir '/blackrock'], [localdir '/labview']);

addpath_recurse([matlabdir '/schmidt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add paths to other packages%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath([matlabdir '/mvgc_v1.0']);
addpath_recurse([matlabdir '/gpfa']);
%addpath([matlabdir '/te_matlab_0.5/']);

mvgcstartup

%%%%%%%%%%%%%%%%%%%%%%%%
%GPFA startup....
%%%%%%%%%%%%%%%%%%

%path(path,[matlabdir 'gpfa/util/invToeplitz']);
% Create the mex file if necessary.
%if ~exist(sprintf('gpfa/util/invToeplitz/invToeplitzFastZohar.%s',mexext),'file')
%  try
%    eval(sprintf('mex -outdir util/invToeplitz util/invToeplitz/invToeplitzFastZohar.c'));
%    fprintf('NOTE: the relevant invToeplitz mex files were not found.  They have been created.\n');
%  catch
%    fprintf('NOTE: the relevant invToeplitz mex files were not found, and your machine failed to create them.\n');
%    fprintf('      This usually means that you do not have the proper C/MEX compiler setup.\n');
%    fprintf('      The code will still run identically, albeit slower (perhaps considerably).\n');
%    fprintf('      Please read the README file, section Notes on the Use of C/MEX.\n');
%  end
%end

% Posterior Covariance Precomputation  
%path(path,[matlabdir 'gpfa/util/precomp']);
% Create the mex file if necessary.
%if ~exist(sprintf('gpfa/util/precomp/makePautoSumFast.%s',mexext),'file')
%  try
%    eval(sprintf('mex -outdir util/precomp util/precomp/makePautoSumFast.c'));
%    fprintf('NOTE: the relevant precomp mex files were not found.  They have been created.\n');
%  catch
%    fprintf('NOTE: the relevant precomp mex files were not found, and your machine failed to create them.\n');
%    fprintf('      This usually means that you do not have the proper C/MEX compiler setup.\n');
%    fprintf('      The code will still run identically, albeit slower (perhaps considerably).\n');
%    fprintf('      Please read the README file, section Notes on the Use of C/MEX.\n');
%  end
%end
%%%%
