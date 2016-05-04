figure

subplot(2,2,1)

x=reshape(normalize_connections(fulltestscan.TE),1,900);
y=reshape(normalize_connections(fulltestscan.real),1,900);
scatter(x,y);
title('Predicted Connections (Covariance) vs actual')
xlabel('Predicted connection (Covariance)')
ylabel('Actual connection')
hold on
plot(0:max(x),0:max(x),'r')

subplot(2,2,2)

x=reshape(normalize_connections(fulltestscan.GLMtest),1,900);
y=reshape(normalize_connections(fulltestscan.real),1,900);
scatter(x,y);
title('Predicted Connections (Granger LM) vs actual')
xlabel('Predicted connection (Granger LM)')
ylabel('Actual connection')
hold on
plot(0:max(x),0:max(x),'r')

subplot(2,2,3)

x=reshape(normalize_connections(fulltestscan.GLMpoissonlin),1,900);
y=reshape(normalize_connections(fulltestscan.real),1,900);
scatter(x,y);
title('Predicted Connections (Granger GLM) vs actual')
xlabel('Predicted connection (Granger GLM)')
ylabel('Actual connection')
hold on
plot(0:max(x),0:max(x),'r')

subplot(2,2,4)

x=reshape(normalize_connections(fulltestscan.TE),1,900);
y=reshape(normalize_connections(fulltestscan.real),1,900);
scatter(x,y);
title('Predicted Connections (TE) vs actual')
xlabel('Predicted connection (TE)')
ylabel('Actual connection')
hold on
plot(0:max(x),0:max(x),'r')