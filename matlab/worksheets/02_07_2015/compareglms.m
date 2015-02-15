function compareglms(baseglm, otherglms, threshdev)
	if (nargin < 1) baseglm = 'sprc'; end
	if (nargin < 2) otherglms = {'sprc_pos_lv'}; end
	if (nargin < 3) threshdev = 1e20; end

	fn_out = './worksheets/01_17_2015/weeklynevs.mat';
	conds = {'2D Manual Position', 'Dual Control'};
	load(fn_out);
	basedevs = {};
	basetopten = {};
	otherdevs = {};
	othertopten = {};
	diffdevs = {};
	difftopten = {};
	meandevs = {};
	nE = 128;

	%%%%%%%%%%%%%%%%%%
	%Read in base glm%
	%%%%%%%%%%%%%%%%%%
	modeltype = baseglm;
	%Load each file and extract the deviances for each unit
	for i = 1:length(conds)
		condition = conds{i};
		nevs = weeklynevs(condition);
		nW = length(nevs);
		basedevs{i} = zeros(nW, nE);
		basetopten{i} = zeros(nW, 1);
		for j = 1:nW
			recinfo = nevs(j);
			%Get files
			nevfile = ['./blackrock/' recinfo.nevfile];
			matfile = ['./labview/' recinfo.matfile];
			currdate = strrep(recinfo.date, '/', '-');
			curr_fn = ['./worksheets/02_07_2015/glmresults/' modeltype '_' condition '_' currdate '.mat'];
			curr_fn = strrep(curr_fn, ' ', '_');
			if exist(curr_fn, 'file')
				%Add data to matrix for plotting
				display([curr_fn ' exists, reading in deviances.']);
				load(curr_fn);
				nU = size(model.b_hat,1);
				%Extract deviances
				curr_devs = [];
				for k = 1:nU
					curr_devs(k) = model.dev{k};
				end
				curr_devs(curr_devs>threshdev) = 0;
				%Find units used
				units = cellfun(@str2num, processed_mua.unitnames);
				%Add these to devs matrix
				basedevs{i}(j,units) = curr_devs;
				tt = sort(curr_devs, 'descend');
				basetopten{i}(j) = sum(tt(1:8));
			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Read in the other glms to be compared to the base%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for idx = 1:length(otherglms)
		modeltype = otherglms{idx};
		%Load each file and extract the deviances for each unit
		for i = 1:length(conds)
			condition = conds{i};
			nevs = weeklynevs(condition);
			nW = length(nevs);
			otherdevs{idx,i} = zeros(nW, nE);
			othertopten{idx,i} = zeros(nW, 1);
			diffdevs{idx,i} = zeros(nW, nE);
			difftopten{idx,i} = zeros(nW, 1);
			meandevs{idx,i} = zeros(nW, 1);
			for j = 1:nW
				recinfo = nevs(j);
				%Get files
				nevfile = ['./blackrock/' recinfo.nevfile];
				matfile = ['./labview/' recinfo.matfile];
				currdate = strrep(recinfo.date, '/', '-');
				curr_fn = ['./worksheets/02_07_2015/glmresults/' modeltype '_' condition '_' currdate '.mat'];
				curr_fn = strrep(curr_fn, ' ', '_');
				if exist(curr_fn, 'file')
					%Add data to matrix for plotting
					display([curr_fn ' exists, reading in deviances.']);
					load(curr_fn);
					nU = size(model.b_hat,1);
					%Extract deviances
					curr_devs = [];
					for k = 1:nU
						curr_devs(k) = model.dev{k};
					end
					curr_devs(curr_devs>threshdev) = 0;
					%Find units used
					units = cellfun(@str2num, processed_mua.unitnames);
					%Add these to devs matrix
					otherdevs{idx,i}(j,units) = curr_devs;
					diffdevs{idx,i}(j,units) = basedevs{i}(j,units)-curr_devs;
					meandevs{idx,i}(j)=mean(diffdevs{idx,i}(j,units));
					tt = sort(basedevs{i}(j,units)-curr_devs, 'descend');
					othertopten{idx,i}(j) = sum(tt(1:8));
					difftopten{idx,i}(j) = sum(tt(1:8));
				end
			end
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%
	%Plot the differences%
	%%%%%%%%%%%%%%%%%%%%%%
	for i = 1:length(conds)
		clf
		condition = conds{i};
		nW = length(nevs);
		for idx = 1:length(otherglms)
			modeltype = otherglms{idx};
			subplot(3,1,1)
		    colormap(hot);
		    imagesc(diffdevs{idx,i}');
		    title(['Deviance between ' modeltype ' and ' baseglm ' for ' condition ' recordings']);
		    ylabel('Electrode');
		    colorbar;
		    %Look at sum of top 10 performing units for each recording
		    subplot(3,1,2)
		    plot(difftopten{idx,i}, '.')
		    ylabel('Sum of largest 8 deviance changes')
		    %Look at boxplot of all units for each week
			subplot(3,1,3)
			[X,Y] = meshgrid(1:nW, 1:nE);
		    grps = reshape(X, nW*nE,1);
		    rshpdevs = reshape(diffdevs{idx,i}, nW*nE,1);
		    boxplot(rshpdevs(rshpdevs>0), grps(rshpdevs>0));
		    xlabel('Week')
		    ylabel('Deviance')
		    fn_out = ['./worksheets/02_07_2015/plots/' modeltype '-' baseglm '_' condition '.eps'];
		    saveplot(gcf, fn_out, 'eps', [6 15])    
		end
	    clf
	    hold on
	    legstr = {};
	    colormap('lines');
	    cm = colormap;
	    for idx = 1:length(otherglms)
		    plot(meandevs{idx,i}, '.', 'Color', cm(idx,:));
		    legstr = {legstr{:}, otherglms{idx}};
	    end
	    title(['Mean change in deviance for ' condition ' recordings']);
	    legend(legstr);
	    xlabel('Week')
	    ylabel('Deviance')
	    fn_out = ['./worksheets/02_07_2015/plots/' condition '_mean.eps'];
	    saveplot(gcf, fn_out)    		    
	end
end