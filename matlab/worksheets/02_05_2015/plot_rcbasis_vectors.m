const = 'on';
fn_out = './worksheets/02_05_2015/20130117SpankyUtah001';
nK_sp = 200; 
nK_pos = 5;
dt_pos = 0.2;
dt_sp = 0.002;
%Load test preprocessed data
pre = load('./testdata/test_preprocess_spline.mat');
data = filters_sprc_pos(pre.processed, nK_sp, nK_pos);
%Fit model
model = MLE_glmfit(data, const);
%Plot filters
fn_out = './worksheets/02_05_2015/20130117SpankyUtah001_rc';
plot_filters_rc(model, data, pre.processed, fn_out);
%Compare with usual basis
datasp = filters_sp_pos(pre.processed, nK_sp, nK_pos);
%Fit model
modelsp = MLE_glmfit(datasp, const);
%Plot filters
fn_out = './worksheets/02_05_2015/20130117SpankyUtah001_sp';
plot_filters(modelsp, datasp, pre.processed, fn_out);

nU = size(data.X,1);
%Plot each model_rc filter vs each model_sp filter
for idx=1:nU 
	clf;
	dt_filt = data.k{1,3};
	b_hat = model.b_hat(idx,:);
	b_hat_sp = modelsp.b_hat(idx,:);
	filt = b_hat_sp(1+datasp.sp_hist);
	ttsp = (0:length(filt)-1)*dt_filt*1000;
	%Extract spike history data
	rc_filt = b_hat(1+data.sp_hist);
	%Transform back to original basis
	sp_filt = data.spbasis*rc_filt';
	tt = (0:length(sp_filt)-1)*dt_filt*1000;
	plot(tt, sp_filt, 'b', ttsp, filt, 'r');
	if ischar(pre.processed.unitnames)
		name = pre.processed.unitnames;
	else
		name = pre.processed.unitnames{idx};
	end
	legend('RC basis func', 'Standard basis func.', 'Location', 'SouthWest');
	title(['Unit: ' name ' no. of spikes: ' num2str(sum(pre.processed.binnedspikes(:,idx)))]);
	xlabel('time (ms)');
	%save eps
	saveplot(gcf, [fn_out '_unit_' name '_filter_comp.eps'], 'eps', [12,6]);	
end