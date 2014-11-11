binsize = 0.002;
dt_sp = binsize;
maxlag = 4/dt_sp;
nevfile = './testdata/20130117SpankyUtah001.nev';
threshold = 5; offset = 0;
pre = preprocess_spline(nevfile, binsize, threshold, offset);

autocorrRU = xcorr(pre.torque(:,1), maxlag);
autocorrFE = xcorr(pre.torque(:,2), maxlag);
autocorrdRU = xcorr(pre.dtorque(:,1), maxlag);
autocorrdFE = xcorr(pre.dtorque(:,2), maxlag);

figure 
subplot(2,2,1)
tt = ((1:length(autocorrRU))-length(autocorrRU)/2)*dt_sp;
plot(tt,autocorrRU)
xlabel('time(s)')
title('Auto-correlation RU')
subplot(2,2,2)
tt = ((1:length(autocorrFE))-length(autocorrFE)/2)*dt_sp;
plot(tt,autocorrFE)
xlabel('time(s)')
title('Auto-correlation FE')
subplot(2,2,3)
plot(tt,autocorrdFE)
title('Auto-correlation vel RU')
xlabel('time(s)')
subplot(2,2,4)
tt = ((1:length(autocorrdRU))-length(autocorrdRU)/2)*dt_sp;
plot(tt,autocorrdRU)
xlabel('time(s)')
title('Auto-correlation vel FE')
fn_out4 = ['./worksheets/11_11_2014/plots/autocorr_spline.eps'];
saveplot(gcf, fn_out4, 'eps', [6 4]);
