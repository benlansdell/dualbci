files = {'20131031SpankyUtah014', '20131031SpankyUtah016', '20131205SpankyUtah001', '20131205SpankyUtah005' ...
			'20131104SpankyUtah001','20131104SpankyUtah005','20131011SpankyUtah001','20131011SpankyUtah004'};
matfile = {'Spanky_2013-12-05-1400.mat', 'Spanky_2013-11-04-1321.mat', 'Spanky_2013-10-31-1322.mat', 'Spanky_2013-10-11-1440.mat', 'Spanky_2013-09-30-1342.mat'};

files = {'20131031SpankyUtah014'};
matfile = {'Spanky_2013-10-31-1322.mat'};

%Read in spreadsheet
[data, txt] = xlsread('./exptinfo/Results1_excluding2d.xls');

indices = 3:size(data,1);
%Compute tuning angle for MC pos
tunedidx = data(indices,3);

torqueFE1 = data(indices,10);
torqueRU1 = data(indices,11);
torqueFE2 = data(indices,12);
torqueRU2 = data(indices,13);

torqueH1 = atan(torqueRU1./torqueFE1); %theta
torqueS1 = sqrt(torqueRU1.^2 + torqueFE1.^2); %speed
piplus = torqueFE1 < 0 & torqueRU1 > 0;
piminus = torqueFE1 < 0 & torqueRU1 < 0;
torqueH1(piplus) = torqueH1(piplus) + pi;
torqueH1(piminus) = torqueH1(piminus) - pi;
torqueH1(torqueH1<0) = torqueH1(torqueH1<0)+2*pi;

torqueH2 = atan(torqueRU2./torqueFE2); %theta
torqueS2 = sqrt(torqueRU2.^2 + torqueFE2.^2); %speed
piplus = torqueFE2 < 0 & torqueRU2 > 0;
piminus = torqueFE2 < 0 & torqueRU2 < 0;
torqueH2(piplus) = torqueH2(piplus) + pi;
torqueH2(piminus) = torqueH2(piminus) - pi;
torqueH2(torqueH1<0) = torqueH2(torqueH1<0)+2*pi;

%Compute angles used in BC task
taskH1 = (data(indices, 5)==1)*pi+(data(indices, 5)==2)*pi/2;
taskH2 = (data(indices, 5)==2)*0+3*(data(indices, 5)==2)*pi/2;

%Compute difference
deltaH1 = mod(torqueH1 - taskH1, 2*pi);
deltaH2 = mod(torqueH2 - taskH2, 2*pi);

deltaH1(deltaH1>pi) = abs(deltaH1(deltaH1>pi)-2*pi);
deltaH2(deltaH2>pi) = abs(deltaH2(deltaH2>pi)-2*pi);

%Plot performance as function of difference for each unit
performance1BC = data(indices,25);

%subplot(2,1,1)
%plot(deltaH1, performance1BC, '.');
%xlabel('\Delta unit 1')
%ylabel('Targets/minute')

%subplot(2,1,2)
%plot(deltaH2, performance1BC, '.');
%xlabel('\Delta unit 2')
%ylabel('Targets/minute')

%pts that are control = 5
%

clf
colormap(jet);
hold on
scatter(deltaH1(tunedidx==1), deltaH2(tunedidx==1), [], performance1BC(tunedidx==1), '.')
scatter(deltaH1(tunedidx==2), deltaH2(tunedidx==2), [], performance1BC(tunedidx==2), 'o')
scatter(deltaH1(tunedidx==3), deltaH2(tunedidx==3), [], performance1BC(tunedidx==3), 'd')
scatter(deltaH1(tunedidx==4), deltaH2(tunedidx==4), [], performance1BC(tunedidx==4), '+')
scatter(deltaH1(tunedidx==5), deltaH2(tunedidx==5), [], performance1BC(tunedidx==5), '*')
%legend('Tuned', 'Tuned ff', 'Tuned 45', 'Untuned', 'Control', 'location', 'northoutside')
colorbar
xlabel('\Delta unit 1')
ylabel('\Delta unit 2')

saveplot(gcf, './worksheets/05_27_2015/performance_scatter.eps')

%Plot histograms of preferred directions computed by cross correlation
clf
subplot(2,1,1)
hist(torqueH1, 25)
subplot(2,1,2)
hist(torqueH2, 25)
xlabel('Preferred angle')

saveplot(gcf, './worksheets/05_27_2015/preferreddirection_crosscorrelation.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot linear regression preferred direction%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[datalinear, txt] = xlsread('./exptinfo/ResultsTuningPosRegression.xls');

indices = 3:size(datalinear,1);
%Compute tuning angle for MC pos
tunedidx = datalinear(indices,4);
torqueH1 = datalinear(indices,5);
torqueH2 = datalinear(indices,6);

clf
subplot(2,1,1)
hist(torqueH1, 25)
subplot(2,1,2)
hist(torqueH2, 25)
xlabel('Preferred angle')

saveplot(gcf, './worksheets/05_27_2015/preferreddirection_linearregression.eps')

%Plot performance versus change in 
[data, txt] = xlsread('./exptinfo/Results1_excluding2d.xls');

indices = 3:size(data,1);
%Compute tuning angle for MC pos
tunedidx = data(indices,3);

torqueFE1 = data(indices,34);
torqueRU1 = data(indices,35);
torqueFE2 = data(indices,36);
torqueRU2 = data(indices,37);

torqueH1 = atan(torqueRU1./torqueFE1); %theta
torqueS1 = sqrt(torqueRU1.^2 + torqueFE1.^2); %speed
piplus = torqueFE1 < 0 & torqueRU1 > 0;
piminus = torqueFE1 < 0 & torqueRU1 < 0;
torqueH1(piplus) = torqueH1(piplus) + pi;
torqueH1(piminus) = torqueH1(piminus) - pi;
torqueH1(torqueH1<0) = torqueH1(torqueH1<0)+2*pi;

torqueH2 = atan(torqueRU2./torqueFE2); %theta
torqueS2 = sqrt(torqueRU2.^2 + torqueFE2.^2); %speed
piplus = torqueFE2 < 0 & torqueRU2 > 0;
piminus = torqueFE2 < 0 & torqueRU2 < 0;
torqueH2(piplus) = torqueH2(piplus) + pi;
torqueH2(piminus) = torqueH2(piminus) - pi;
torqueH2(torqueH1<0) = torqueH2(torqueH1<0)+2*pi;

%Compute angles used in BC task
taskH1 = (data(indices, 5)==1)*pi+(data(indices, 5)==2)*pi/2;
taskH2 = (data(indices, 5)==2)*0+3*(data(indices, 5)==2)*pi/2;

%Compute difference
deltaH1 = mod(torqueH1 - taskH1, 2*pi)*180/pi;
deltaH2 = mod(torqueH2 - taskH2, 2*pi)*180/pi;

deltaH1(deltaH1>180) = abs(deltaH1(deltaH1>180)-360);
deltaH2(deltaH2>180) = abs(deltaH2(deltaH2>180)-360);

%Plot performance as function of difference for each unit
performance1BC = data(indices,25);

%subplot(2,1,1)
%plot(deltaH1, performance1BC, '.');
%xlabel('\Delta unit 1')
%ylabel('Targets/minute')

%subplot(2,1,2)
%plot(deltaH2, performance1BC, '.');
%xlabel('\Delta unit 2')
%ylabel('Targets/minute')

%pts that are control = 5
%

clf
colormap(jet);
hold on
scatter(deltaH1(tunedidx==1), deltaH2(tunedidx==1), [], performance1BC(tunedidx==1), '.')
scatter(deltaH1(tunedidx==2), deltaH2(tunedidx==2), [], performance1BC(tunedidx==2), 'o')
scatter(deltaH1(tunedidx==3), deltaH2(tunedidx==3), [], performance1BC(tunedidx==3), 'd')
scatter(deltaH1(tunedidx==4), deltaH2(tunedidx==4), [], performance1BC(tunedidx==4), '+')
scatter(deltaH1(tunedidx==5), deltaH2(tunedidx==5), [], performance1BC(tunedidx==5), '*')
%legend('Tuned', 'Tuned ff', 'Tuned 45', 'Untuned', 'Control', 'location', 'northoutside')
colorbar
xlabel('\Delta unit 1')
ylabel('\Delta unit 2')

saveplot(gcf, './worksheets/05_27_2015/performance_scatter_regression.eps')
