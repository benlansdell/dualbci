function [GCdevperm, GCpvalperm, GCsigperm, clusters, namesperm] = granger_cluster(GCdev, GCpval, GCsig, unitnames, fn_out, pval);
	%Computes Granger-causality matrices for set of spike and cursor trajectory data. Clusters
	%using MCL algorithm (mcl program)
	%
	%Usage:
	%	[GCdev, GCpval, GCsig] = granger_cluster(GCdev, GCpval, GCsig, unitnames, fn_out);
	%     
	%Input:
	%	GCdev = [nU x nU] square matrix with granger causality measure (change in deviance)
	%	GCpval = square matrix with pvalues corresponding to change in deviance
	%	GCsig = square matrix listing result of significance test (presumed to be equal to pval parameter)
	%	unitnames = [nU x 1] cell array with names of units
	%	fn_out = output filename for eps file
	%	pval = (optional, default = 0.05) default pvalue below which to set edge weight to zero
	%   
	%Output:
	%	GCdev = [nU x nU] matrix where i,j element is the change in deviance when unit i is excluded
	%		from a model of unit j. (i.e. unit i's effect on unit j). Permuted so in clusters
	%	GVpval = [nU x nU] matrix listing p-value of each GCdev value
	%	GVsig = [nU x nU] matrix listing significance of GCdev
	%  
	%Test code:
	%	%Load test preprocessed data
	%	const = 'on';
	%	pval = 0.001;
	%	nK_sp = 6; 
	%	nK_pos = 6;
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_60hz_short24.mat');
	%	fn_out = './worksheets/12_2_2014/plots/testgranger.eps';
	%	processed = pre.processed;
	%	data = filters_sp_pos_network(processed, nK_sp, nK_pos);
	%	[GCdev, GCpval, GCsig] = granger(processed, data, fn_out, pval);
	%	granger_cluster(GCdev, GCpval, GCsig, processed.unitnames, [fn_out '.cluster']);

	if (nargin < 6) pval = 0.001; end

	nU = size(unitnames,2);
	tmpfile = './tmp.adj';
	clusters = {};

	mindev = 1e8;
	inflation = 2.5;
	%Make a file appropriate for mcl program
	%Also: set main diagonal of GCdev to zero (not interested in self-similarity)
	fid = fopen(tmpfile, 'w');
	for i = 1:nU
		for j = 1:nU
			if (i == j)
				GCdev(i,j) = 0;
			%else
			elseif (GCpval(i,j)<pval)
			%if (GCpval(i,j)<pval)
				mindev = min(mindev,GCdev(i,j));
				fprintf(fid, '%d %d %f\n', i, j, GCdev(i,j));
			end
		end
	end
	fclose(fid);
	%Run mcl program
	cmd = ['mcl ./tmp.adj --abc -I ' num2str(inflation) ' -tf ''gq(' num2str(mindev) ')'' -o out.adj']
	stat = system(cmd);

	%Load back in results from mcl program
	fid = fopen('./out.adj');
	tline = fgetl(fid);
	i = 1;
	order = [];
	while ischar(tline)
		cluster = [];
		words = textscan(tline, '%s');
		words = words{1};
		for idx = 1:length(words)
			cluster = [cluster, str2num(words{idx})];
			order = [order, str2num(words{idx})];
		end
		clusters{i} = cluster;
		i = i+1;
	    tline = fgetl(fid);
	end
	fclose(fid);

	%Permute matrices and unitnames by results
	namesperm = unitnames(order);
	GCdevperm = GCdev(order, order);
	GCpvalperm = GCpval(order, order);
	GCsigperm = GCsig(order, order);

	%Plot results
	clf
	subplot(3,1,1)
	colormap(bone);
	imagesc(GCdevperm)
	title('Change in deviance')
	ylabel('Unit')
	xlabel('Unit')
	set(gca,'XTick',1:nU);
	set(gca,'YTick',1:nU);
	set(gca,'XTickLabel',namesperm);
	set(gca,'YTickLabel',namesperm);
	rotateXLabels(gca, 90);
	colorbar
	subplot(3,1,2)
	imagesc(1-GCpvalperm)
	title(['p-value'])
	ylabel('Unit')
	xlabel('Unit')
	set(gca,'XTick',1:nU);
	set(gca,'YTick',1:nU);
	set(gca,'XTickLabel',namesperm);
	set(gca,'YTickLabel',namesperm);
	%rotate x labels
	rotateXLabels(gca, 90);
	colorbar
	subplot(3,1,3)
	imagesc(GCsigperm)
	title(['Significance test @ p<' num2str(pval)])
	ylabel('Unit')
	xlabel('Unit')
	set(gca,'XTick',1:nU);
	set(gca,'YTick',1:nU);
	set(gca,'XTickLabel',namesperm);
	set(gca,'YTickLabel',namesperm);
	rotateXLabels(gca, 90);
	colorbar
	saveplot(gcf, fn_out, 'eps', [6 18]);