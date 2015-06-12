%Script to run simple linear regression on manual control and 1D brain control
%Quantify the change in rotation of BC units between brain and manual control
%Determine relationship between performance and angle of rotation...

script = './worksheets/2015_06_10/lineartuning.m';
binsize = 0.1;
threshold = 5;
blackrock = './blackrock/';
%100ms offset (torque data is moved 100ms ahead of spiking data)
offset = -0.1;
nK_pos = 1;
nK_sp = 0;
const = 'on';
%Fetch each pair of nev files to run
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
tablename = 'AnalysisLinear';
colnames = {'1DBCrecording', 'manualrecording'};
toprocess = exec(conn, ['SELECT `1DBCrecording`,`manualrecording` FROM AnalysisLinear WHERE regrFE1 IS NULL']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);
for idx = 1:nR
	bcifile = toprocess{idx, 1};
	manfile = toprocess{idx, 2};
	bcipath = [blackrock bcifile];
	manpath = [blackrock manfile];
	%Get the mat file to use
	matfile = fetch(exec(conn, ['SELECT `labview file` FROM Recordings WHERE `nev file` = "' bcifile '"']));
	matfilebci = matfile.Data;
	matfile = fetch(exec(conn, ['SELECT `labview file` FROM Recordings WHERE `nev file` = "' manfile '"']));
	matfileman = matfile.Data;
	%Check they're from the same recording session
	if ~strcmp(matfileman, matfilebci)
		display(['Recordings are from a different trial session, check your data for: ' bcifile])
		continue
	end
	matfile = ['./labview/' matfileman{1}];
	%Get the BC units from the BCI recording
	BCunits = fetch(exec(conn, ['SELECT `unit`, `direction` FROM BCIUnits WHERE `ID` = "' bcifile '"']));
	BCunits = BCunits.Data;
	if size(BCunits,1) ~= 2
		display([bcifile ' is not configured as a 1D brain control task. Continuing'])
		continue
	end
	BCunit1 = num2str(BCunits{1,1});
	BCunit2 = num2str(BCunits{2,1});
	BCdir1 = taskTheta(BCunits{1,2});
	BCdir2 = taskTheta(BCunits{2,2});
	%Preprocess data
	processedbci = preprocess_spline_lv(bcipath, matfile, binsize, threshold, offset);
	processedman = preprocess_spline_lv(manpath, matfile, binsize, threshold, offset);
	%Truncate to just the units of interest
	ind1 = find(ismember(processedbci.unitnames,BCunit1));
	ind2 = find(ismember(processedbci.unitnames,BCunit2));
	processedbci.unitnames = processedbci.unitnames([ind1, ind2]);
	processedbci.binnedspikes = processedbci.binnedspikes(:,[ind1, ind2]);
	processedbci.rates = processedbci.rates(:,[ind1, ind2]);
	ind1 = find(ismember(processedman.unitnames,BCunit1));
	ind2 = find(ismember(processedman.unitnames,BCunit2));
	processedman.unitnames = processedman.unitnames([ind1, ind2]);
	processedman.binnedspikes = processedman.binnedspikes(:,[ind1, ind2]);
	processedman.rates = processedman.rates(:,[ind1, ind2]);
	%Fit a linear model to the pair
	data = filters_sp_pos_lv(processedman, nK_sp, nK_pos);
	modelman = MLE_lmfit(data, const);
	data = filters_sp_pos_lv(processedbci, nK_sp, nK_pos);
	modelbci = MLE_lmfit(data, const);
	%Extract and save regression coefficients
	regrFE1 = modelman.b_hat(1,2);
	regrRU1 = modelman.b_hat(1,3);
	regrFE2 = modelman.b_hat(2,2);
	regrRU2 = modelman.b_hat(2,3);
	regrBCI1 = modelbci.b_hat(1,2);
	regrBCI2 = modelbci.b_hat(2,2);
	%Compute preferred tuning angle and size (magnitude of regression coeffs)
	%Compute rotation angle between units
	[deltaH1, mantheta1, size1] = changeTheta(regrFE1, regrRU1, BCdir1);
	[deltaH2, mantheta2, size2] = changeTheta(regrFE2, regrRU2, BCdir2);
	%Tag with computer run on, date, last git commit, and script name
	host = hostname();
	comm = currCommit();
	stamp = datestr(now);
	%Save results in database
	colnames = {'regrFE1', 'regrRU1', 'regrFE2', 'regrRU2', 'regrBCI1' , 'regrBCI2', 'deltaangle1', 'deltaangle2', 'angle1', 'tuningsize1', 'angle2',...
	'tuningsize2', 'computer', '`analysis date`', 'commit', 'scriptname'};
	sqldata = {regrFE1, regrRU1, regrFE2, regrRU2, regrBCI1, regrBCI2, deltaH1, deltaH2, mantheta1, size1, mantheta2, size2, host, stamp, comm, script};
	whereclause = ['WHERE `1DBCrecording` = "' bcifile '"'];
	update(conn,tablename,colnames,sqldata,whereclause)
end