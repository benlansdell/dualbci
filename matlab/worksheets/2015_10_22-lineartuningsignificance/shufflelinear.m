%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 1;
blackrock = './blackrock/';
labviewpath = './labview/';

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
tablename = '`experiment_tuning`';
colnames = '`1DBCrecording`, `manualrecording`, `manualrecordingafter`';
toprocess = exec(conn, ['SELECT ' colnames ' FROM ' tablename]);
%toprocess = exec(conn, ['SELECT `1DBCrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);
nR = min(nR, 10);

sizes_MC = [];
sizes_MC2 = [];
sizes_BC = [];

nS = 100;

for idx = 1:nR
	BCnevfile = toprocess{idx,1};
	nevfile1 = toprocess{idx,2};
	nevfile2 = toprocess{idx,3};
	display([num2str(idx) '/' num2str(nR) ': Processing ' nevfile1])
	if exist([blackrock nevfile1], 'file') & exist([blackrock BCnevfile], 'file') & exist([blackrock nevfile2], 'file')
		%Manual control 1
		nevpath = [blackrock nevfile1];
		%Load parameters
		eval(paramcode);
		%Preprocess data
		processed = preprocess_smooth(nevpath, binsize, sigma_fr, sigma_trq, threshold, offset);
		for j = 1:nS
			display(num2str(j))
			%Shuffle output
			shuffle = randperm(size(processed.binnedspikes,1));
			processed.binnedspikes = processed.binnedspikes(shuffle, :);
			%Fit a linear model to training data
			data = filters_sp_pos(processed, nK_sp, nK_pos);
			model = MLE_lmfit(data, const);
			%Compute size
			sizes = sqrt(model.b_hat(:,2).^2+model.b_hat(:,3).^2);
			sizes_MC = [sizes_MC; sizes];

			%Extract deviance
			dev = model.dev{idx, 1};
			%Extract SE
			stats = model.stats{idx};
			%Extract MSE. Need to add cross-validation code to do this... add later
			b_hat = model.b_hat(idx,:);
			datanovel = filters_sp_pos(processed_novel, nK_sp, nK_pos);
			rho_hat = glmval(b_hat', squeeze(data.X(idx,:,:)), 'identity');
			rho_hat_novel = glmval(b_hat', squeeze(datanovel.X(idx,:,:)), 'identity');
		end

		%Brain control
		nevpath = [blackrock BCnevfile];
		%Load parameters
		eval(paramcode);
		%Preprocess data
		processed = preprocess_smooth(nevpath, binsize, sigma_fr, sigma_trq, threshold, offset);
		for j = 1:nS
			%Shuffle output
			shuffle = randperm(size(processed.binnedspikes,1));
			processed.binnedspikes = processed.binnedspikes(shuffle, :);
			%Fit a linear model to training data
			data = filters_sp_pos(processed, nK_sp, nK_pos);
			model = MLE_lmfit(data, const);
			%Compute size
			sizes = sqrt(model.b_hat(:,2).^2+model.b_hat(:,3).^2);
			sizes_BC = [sizes_BC; sizes];
		end

		%Manual control 2
		nevpath = [blackrock nevfile2];
		%Load parameters
		eval(paramcode);
		%Preprocess data
		processed = preprocess_smooth(nevpath, binsize, sigma_fr, sigma_trq, threshold, offset);
		for j = 1:nS
			%Shuffle output
			shuffle = randperm(size(processed.binnedspikes,1));
			processed.binnedspikes = processed.binnedspikes(shuffle, :);
			%Fit a linear model to training data
			data = filters_sp_pos(processed, nK_sp, nK_pos);
			model = MLE_lmfit(data, const);
			%Compute size
			sizes = sqrt(model.b_hat(:,2).^2+model.b_hat(:,3).^2);
			sizes_MC2 = [sizes_MC2; sizes];
		end
	else
		display('Cannot find all files, continuing')
	end
end

save('./worksheets/2015_10_22-lineartuningsignificance/shuffedsizes.mat', 'sizes_MC' , 'sizes_BC', 'sizes_MC2');

clf
hist(sizes_MC)
xlabel('Tuning strength')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/sizes_MC.eps')

clf
hist(sizes_BC)
xlabel('Tuning strength')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/sizes_BC.eps')

clf
hist(sizes_MC2)
xlabel('Tuning strength')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/sizes_MC2.eps')

