%Partial least square expt
%Loading parameters and file info
modelID = 13;
verbose = 1;
blackrock = './blackrock/';
labview = './labview/';
scriptname = './2015_07_28-PLSexpt/plstest.m';
scriptdesc = ['Testing PLS.'];
%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
paramcode = exec(conn, ['SELECT `description` FROM Models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};
%Fetch files to analyze
tablename = '`AnalysisLinear`';
toprocess = exec(conn, ['SELECT `manualrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);
idx = 1;
nevfile = toprocess{idx};

%Prepare data
%Load paramcode
eval(paramcode);
N_sub = 5000;
nevpath = [blackrock nevfile];
%Load parameters
eval(paramcode);
matfile = fetch(exec(conn, ['SELECT `labview file`,`duration` FROM `Recordings` rec WHERE rec.`nev file` = "' nevfile '"']));
duration = matfile.Data{2};
matfile = matfile.Data{1};
matpath = [labview matfile];

%Load units, preprocess, and prepare matrices
searchstr = ['SELECT unit FROM `Units` WHERE `nev file` = "' nevfile '"'...
' ORDER BY `Units`.`firingrate` DESC LIMIT ' num2str(nU)];
units = fetch(exec(conn, searchstr));
units = units.Data;
%Preprocess
processed = preprocess_spline_lv(nevpath, matfile, binsize, threshold, offset, [], [], units);
%Truncate to specified duration
[processed, processed_novel] = split_recording(processed, dur, dur+testdur);
nUtotal = length(processed.unitnames);
%Prepare data
data = filters_sprc_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);

%Run PLS    