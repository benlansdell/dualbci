%Script to generate a bunch of test data
freqlow = 20;
freqhigh = 60;
%freqhigh = 2000;
N = 100000;
binsize = 0.002;
dur = N*binsize;
dt_sp = binsize;
dt_pos = 0.05;
fn_out = './testdata/glm_test1.mat'

%Filters
k_const = [0];
nK_sp = 50;
nK_pos = 6;

t_sp = linspace(0,1,nK_sp);
%k_sp = 1.5*exp(-4*t_sp)-2*exp(-6*t_sp);
k_sp = -2*exp(-6*t_sp);
k_sp = fliplr(k_sp);

t_pos = linspace(0,1,nK_pos);
k_RU = -4*25*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE = 6*25*exp(-5*t_pos);%.*sin(t_pos*2*pi);

b = [k_const, k_sp, k_RU, k_FE];

%Simulate GLM
processed = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU, k_FE, N, binsize);
save(fn_out, 'processed');

%Estimate filters
fn_out = './worksheets/09_15_2014/plots/testdata_large_sp';
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
k = data.k;
nK = size(k,1);
const = 'on';
model = MLE_glmfit(data, const);
stats = model.stats{1};
b_hat = model.b_hat(1,:);
k_const_hat = b_hat(1);
rho = glmval(b', squeeze(data.X(1,:,:)), 'log');
rho_hat = glmval(b_hat', squeeze(data.X(1,:,:)), 'log');

%Plot stuff
nP = 3;
%clf
figure
for j = 1:nK
	%Extract data
	name = k{j,1};
	filt = b_hat(k{j,2}+1);
	filt_true = b(k{j,2}+1);
	se = stats.se(k{j,2}+1)';
	tstat = stats.t(k{j,2}+1);
	pval = stats.p(k{j,2}+1);
	dt_filt = k{j,3};

	%If filter length is zero skip this one
	if length(k{j,2}) < 1
		continue
	end
	%Plot filter plus/minus SE
	subplot(nP, nK+1, j)
	tt = (0:length(filt)-1)*dt_filt*1000;
	hold on
	ymin = min(filt-se)*1.2;
	ymax = max(filt+se)*1.2;
	area(tt, filt+se, ymin, 'FaceColor', [0.8 0.8 0.8])
	area(tt, filt-se, ymin, 'FaceColor', [1 1 1])
	plot(tt, filt);
	plot(tt, filt_true, 'r');
	ylim([ymin ymax]);
	if length(tt) > 1
		xlim([min(tt) max(tt)]);
	end
	title(name);
	%xlabel('time (ms)');
	%Plot tstats
	subplot(nP, nK+1, (nK+1)+j)
	plot(tt, tstat);
	if (j==1)
		ylabel('t statistic');
	end
	%xlabel('time (ms)');
	%p-values
	subplot(nP, nK+1, (nK+1)*2+j)
	hold on
	above = pval > 0.05;
	below = pval <= 0.05;
	plot(tt(above), log(pval(above)), '.b');
	plot(tt(below), log(pval(below)), '.r');
	plot(tt, log(0.05)*ones(length(tt),1), 'k')
	if (j==1)
		ylabel('log(p-val)');
	end
	xlabel('time (ms)');
end
saveplot(gcf, [fn_out '_filters.eps'], 'eps', [12,6]);
saveas(gcf, [fn_out '_filters.fig']);


figure
t_i = 30;
t_f = 33;
%t_f = 40; 
ii = 1:size(data.y,2);
tt = ii*processed.binsize;
plot(tt, rho(ii), tt, rho_hat(ii));
xlim([10 13])
xlabel('time (s)')
ylabel('firing rate (spikes/s)')
legend('Actual rate', 'Estimated rate');
saveplot(gcf, [fn_out '_rho.eps'], 'eps', [12,6]);
saveas(gcf, [fn_out '_rho.fig']);