conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tuning strength position%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_data = fetch(exec(conn, ['SELECT fl1.r2, fl2.r2, fl3.r2, fl5.r2, fl1.r2out, fl2.r2out, fl3.r2out, fl5.r2out, et1.`tuning_type` FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
'INNER JOIN `fits` flin3 '...
'ON flin3.`nev file` = et1.`manualrecordingafter`'...
'INNER JOIN `fits_linear` fl3 '...
'ON flin3.id = fl3.id '...
'INNER JOIN `fits` flin5 '...
'ON flin5.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl5 '...
'ON flin5.id = fl5.id '...
'WHERE flin1.modelID = 1 AND flin2.modelID = 1 AND flin3.modelID = 1 AND flin5.modelID = 1 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin5.unit ']));
all_r2 = cell2mat(all_data.Data(:,1:4));
all_r2out = cell2mat(all_data.Data(:,5:8));
tuningtype = cell2mat(all_data.Data(:,9));

corrs(1) = corr(all_r2(:,1), all_r2(:,2));
corrs(2) = corr(all_r2(:,1), all_r2(:,3));
corrs(3) = corr(all_r2(:,3), all_r2(:,2));
corrs(4) = corr(all_r2(:,1), all_r2(:,4));
corrs(5) = corr(all_r2(:,2), all_r2(:,4));

%boxplot(all_r2(:,2), tuningtype)

unrot = tuningtype == 5;
rot = (tuningtype == 1 | tuningtype == 3| tuningtype == 4);
untuned = tuningtype == 4;

corrsrotated(1) = corr(all_r2(rot,1), all_r2(rot,2));
corrsrotated(2) = corr(all_r2(rot,1), all_r2(rot,3));
corrsrotated(3) = corr(all_r2(rot,3), all_r2(rot,2));
corrsrotated(4) = corr(all_r2(rot,1), all_r2(rot,4));
corrsrotated(5) = corr(all_r2(rot,2), all_r2(rot,4));

corrsunrotated(1) = corr(all_r2(unrot,1), all_r2(unrot,2));
corrsunrotated(2) = corr(all_r2(unrot,1), all_r2(unrot,3));
corrsunrotated(3) = corr(all_r2(unrot,3), all_r2(unrot,2));
corrsunrotated(4) = corr(all_r2(unrot,1), all_r2(unrot,4));
corrsunrotated(5) = corr(all_r2(unrot,2), all_r2(unrot,4));

clf
cc = ones(size(tuningtype));
cc(tuningtype == 1 | tuningtype == 3 | tuningtype == 4) = 2;
cc(tuningtype == 5) = 3;
colors = [0 0 0; 1 0 0; 0 0 1];
c = [];
for idx = 1:size(cc, 1);
	c(idx,1:3) = colors(cc(idx),:);
end
colormap(colors);
subplot(2,3,1)
scatter(all_r2(:,1), all_r2(:,2), [], cc)
xlabel('R^2 MC1')
ylabel('R^2 BC1')
title(['corr rotated: ' num2str(corrsrotated(1)) ' corr unrotated: ' num2str(corrsunrotated(1))])
subplot(2,3,2)
scatter(all_r2(:,1), all_r2(:,3), [], cc)
xlabel('R^2 MC1')
ylabel('R^2 MC2')
title(['corr rotated: ' num2str(corrsrotated(2)) ' corr unrotated: ' num2str(corrsunrotated(2))])
subplot(2,3,3)
scatter(all_r2(:,3), all_r2(:,2), [], cc)
ylabel('R^2 BC1')
xlabel('R^2 MC2')
title(['corr rotated: ' num2str(corrsrotated(3)) ' corr unrotated: ' num2str(corrsunrotated(3))])
subplot(2,3,5)
scatter(all_r2(:,1), all_r2(:,4), [], cc)
xlabel('R^2 MC1')
ylabel('R^2 DC')
title(['corr rotated: ' num2str(corrsrotated(4)) ' corr unrotated: ' num2str(corrsunrotated(4))])
subplot(2,3,6)
scatter(all_r2(:,2), all_r2(:,4), [], cc)
xlabel('R^2 BC1')
ylabel('R^2 DC')
title(['corr rotated: ' num2str(corrsrotated(5)) ' corr unrotated: ' num2str(corrsunrotated(5))])

saveplot(gcf, './worksheets/2015_11_03-bcireportplots/tuningstrengthR2R2out.eps', 'eps', [10 6])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tuning strength velocity%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_data = fetch(exec(conn, ['SELECT fl1.r2, fl2.r2, fl3.r2, fl5.r2, fl1.r2out, fl2.r2out, fl3.r2out, fl5.r2out, et1.`tuning_type` FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
'INNER JOIN `fits` flin3 '...
'ON flin3.`nev file` = et1.`manualrecordingafter`'...
'INNER JOIN `fits_linear` fl3 '...
'ON flin3.id = fl3.id '...
'INNER JOIN `fits` flin5 '...
'ON flin5.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl5 '...
'ON flin5.id = fl5.id '...
'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin3.modelID = 30 AND flin5.modelID = 30 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin5.unit ']));
all_r2 = cell2mat(all_data.Data(:,1:4));
all_r2out = cell2mat(all_data.Data(:,5:8));
tuningtype = cell2mat(all_data.Data(:,9));

corrs(1) = corr(all_r2(:,1), all_r2(:,2));
corrs(2) = corr(all_r2(:,1), all_r2(:,3));
corrs(3) = corr(all_r2(:,3), all_r2(:,2));
corrs(4) = corr(all_r2(:,1), all_r2(:,4));
corrs(5) = corr(all_r2(:,2), all_r2(:,4));

unrot = tuningtype == 5;
rot = (tuningtype == 1 | tuningtype == 3| tuningtype == 4);
untuned = tuningtype == 4;

corrsrotated(1) = corr(all_r2(rot,1), all_r2(rot,2));
corrsrotated(2) = corr(all_r2(rot,1), all_r2(rot,3));
corrsrotated(3) = corr(all_r2(rot,3), all_r2(rot,2));
corrsrotated(4) = corr(all_r2(rot,1), all_r2(rot,4));
corrsrotated(5) = corr(all_r2(rot,2), all_r2(rot,4));

corrsunrotated(1) = corr(all_r2(unrot,1), all_r2(unrot,2));
corrsunrotated(2) = corr(all_r2(unrot,1), all_r2(unrot,3));
corrsunrotated(3) = corr(all_r2(unrot,3), all_r2(unrot,2));
corrsunrotated(4) = corr(all_r2(unrot,1), all_r2(unrot,4));
corrsunrotated(5) = corr(all_r2(unrot,2), all_r2(unrot,4));

clf
cc = ones(size(tuningtype));
%Rotated
cc(tuningtype == 1 | tuningtype == 3 | tuningtype == 4) = 2;
%Unrotated
cc(tuningtype == 5) = 3;
colors = [0 0 0; 1 0 0; 0 0 1];
c = [];
for idx = 1:size(cc, 1);
	c(idx,1:3) = colors(cc(idx),:);
end
colormap(colors);
subplot(2,3,1)
scatter(all_r2(:,1), all_r2(:,2), [], cc)
xlabel('R^2 MC1')
ylabel('R^2 BC1')
title(['corr rotated: ' num2str(corrsrotated(1)) ' corr unrotated: ' num2str(corrsunrotated(1))])
subplot(2,3,2)
scatter(all_r2(:,1), all_r2(:,3), [], cc)
xlabel('R^2 MC1')
ylabel('R^2 MC2')
title(['corr rotated: ' num2str(corrsrotated(2)) ' corr unrotated: ' num2str(corrsunrotated(2))])
subplot(2,3,3)
scatter(all_r2(:,3), all_r2(:,2), [], cc)
ylabel('R^2 BC1')
xlabel('R^2 MC2')
title(['corr rotated: ' num2str(corrsrotated(3)) ' corr unrotated: ' num2str(corrsunrotated(3))])
subplot(2,3,5)
scatter(all_r2(:,1), all_r2(:,4), [], cc)
xlabel('R^2 MC1')
ylabel('R^2 DC')
title(['corr rotated: ' num2str(corrsrotated(4)) ' corr unrotated: ' num2str(corrsunrotated(4))])
subplot(2,3,6)
scatter(all_r2(:,2), all_r2(:,4), [], cc)
xlabel('R^2 BC1')
ylabel('R^2 DC')
title(['corr rotated: ' num2str(corrsrotated(5)) ' corr unrotated: ' num2str(corrsunrotated(5))])

saveplot(gcf, './worksheets/2015_11_03-bcireportplots/tuningstrengthVR2R2out.eps', 'eps', [10 6])

%%%%%%%%%%%%%%%%%%%
%Only plot rotated%
%%%%%%%%%%%%%%%%%%%

rotc = cc == 2;
subplot(2,3,1)
scatter(all_r2(rotc,1), all_r2(rotc,2), 'r')
xlabel('R^2 MC1')
ylabel('R^2 BC1')
title(['corr rotatedrot ' num2str(corrsrotated(1))])
subplot(2,3,2)
scatter(all_r2(rotc,1), all_r2(rotc,3), 'r')
xlabel('R^2 MC1')
ylabel('R^2 MC2')
title(['corr rotatedrot ' num2str(corrsrotated(2))])
subplot(2,3,5)
scatter(all_r2(rotc,1), all_r2(rotc,4), 'r')
xlabel('R^2 MC1')
ylabel('R^2 DC')
title(['corr rotatedrot ' num2str(corrsrotated(4))])

saveplot(gcf, './worksheets/2015_11_03-bcireportplots/tuningstrengthVR2R2out_rotated.eps', 'eps', [10 6])


%Tuning angles (velocity)
all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl3.dir, fl5.dir, et1.`tuning_type`, rec1.`successrate`, rec2.`successrate` FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
'INNER JOIN `fits` flin3 '...
'ON flin3.`nev file` = et1.`manualrecordingafter`'...
'INNER JOIN `fits_linear` fl3 '...
'ON flin3.id = fl3.id '...
'INNER JOIN `fits` flin5 '...
'ON flin5.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl5 '...
'ON flin5.id = fl5.id '...
'INNER JOIN `recordings` rec1 '...
'ON rec1.`nev file` = et1.`1DBCrecording` '...
'INNER JOIN `recordings` rec2 '...
'ON rec2.`nev file` = et1.`dualrecording` '...
'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin3.modelID = 30 AND flin5.modelID = 30 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin5.unit ' ...
'AND fl1.r2 > .005 AND fl3.r2 > .005']));
all_r2 = cell2mat(all_data.Data(:,1:4));
tuningtype = cell2mat(all_data.Data(:,5));
performance = 50*cell2mat(all_data.Data(:,6:7))+10;
bciangle = cell2mat(all_data.Data(:,8:9));

corrsdir(1) = corr(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,2)));
corrsdir(2) = corr(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,3)));
corrsdir(3) = corr(all_r2(:,3), translate_angles(all_r2(:,3), all_r2(:,2)));
corrsdir(4) = corr(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,4)));
corrsdir(5) = corr(all_r2(:,2), translate_angles(all_r2(:,2), all_r2(:,4)));

%boxplot(all_r2(:,2), tuningtype)

unrot = tuningtype == 5;
rot = (tuningtype == 1 | tuningtype == 3| tuningtype == 4);
untuned = tuningtype == 4;

corrsrotated(1) = corr(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,2)));
corrsrotated(2) = corr(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,3)));
corrsrotated(3) = corr(all_r2(rot,3), translate_angles(all_r2(rot,3), all_r2(rot,2)));
corrsrotated(4) = corr(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,4)));
corrsrotated(5) = corr(all_r2(rot,2), translate_angles(all_r2(rot,2), all_r2(rot,4)));

corrsunrotated(1) = corr(all_r2(unrot,1), translate_angles(all_r2(unrot,1), all_r2(unrot,2)));
corrsunrotated(2) = corr(all_r2(unrot,1), translate_angles(all_r2(unrot,1), all_r2(unrot,3)));
corrsunrotated(3) = corr(all_r2(unrot,3), translate_angles(all_r2(unrot,3), all_r2(unrot,2)));
corrsunrotated(4) = corr(all_r2(unrot,1), translate_angles(all_r2(unrot,1), all_r2(unrot,4)));
corrsunrotated(5) = corr(all_r2(unrot,2), translate_angles(all_r2(unrot,2), all_r2(unrot,4)));

clf
cc = ones(size(tuningtype));
cc(tuningtype == 1 | tuningtype == 3 | tuningtype == 4) = 2;
cc(tuningtype == 5) = 3;
colors = [0 0 0; 1 0 0; 0 0 1];
c = [];
for idx = 1:size(cc, 1);
	c(idx,1:3) = colors(cc(idx),:);
end
colormap(colors)
subplot(2,3,1)
scatter(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,2)), [], cc)
xlabel('\theta MC1')
ylabel('\theta BC1')
title(['corr rotated: ' num2str(corrsrotated(1)) ' corr unrotated: ' num2str(corrsunrotated(1))])
subplot(2,3,2)
scatter(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,3)), [], cc)
xlabel('\theta MC1')
ylabel('\theta MC2')
title(['corr rotated: ' num2str(corrsrotated(2)) ' corr unrotated: ' num2str(corrsunrotated(2))])
subplot(2,3,3)
scatter(all_r2(:,3), translate_angles(all_r2(:,3), all_r2(:,2)), [], cc)
ylabel('\theta BC1')
xlabel('\theta MC2')
title(['corr rotated: ' num2str(corrsrotated(3)) ' corr unrotated: ' num2str(corrsunrotated(3))])
subplot(2,3,5)
scatter(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,4)), [], cc)
xlabel('\theta MC1')
ylabel('\theta DC')
title(['corr rotated: ' num2str(corrsrotated(4)) ' corr unrotated: ' num2str(corrsunrotated(4))])
subplot(2,3,6)
scatter(all_r2(:,2), translate_angles(all_r2(:,2), all_r2(:,4)), [], cc)
xlabel('\theta BC1')
ylabel('\theta DC')
title(['corr rotated: ' num2str(corrsrotated(5)) ' corr unrotated: ' num2str(corrsunrotated(5))])
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/tuningangleVR2R2out.eps', 'eps', [10 6])


%Tuning angles (position)
all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl3.dir, fl5.dir, et1.`tuning_type`, rec1.`successrate`, rec2.`successrate` FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
'INNER JOIN `fits` flin3 '...
'ON flin3.`nev file` = et1.`manualrecordingafter`'...
'INNER JOIN `fits_linear` fl3 '...
'ON flin3.id = fl3.id '...
'INNER JOIN `fits` flin5 '...
'ON flin5.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl5 '...
'ON flin5.id = fl5.id '...
'INNER JOIN `recordings` rec1 '...
'ON rec1.`nev file` = et1.`1DBCrecording` '...
'INNER JOIN `recordings` rec2 '...
'ON rec2.`nev file` = et1.`dualrecording` '...
'WHERE flin1.modelID = 1 AND flin2.modelID = 1 AND flin3.modelID = 1 AND flin5.modelID = 1 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin5.unit ' ...
'AND fl1.r2 > .005 AND fl3.r2 > .005']));
all_r2 = cell2mat(all_data.Data(:,1:4));
tuningtype = cell2mat(all_data.Data(:,5));
performance = 50*cell2mat(all_data.Data(:,6:7))+10;
bciangle = cell2mat(all_data.Data(:,8:9));


corrsdir(1) = corr(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,2)));
corrsdir(2) = corr(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,3)));
corrsdir(3) = corr(all_r2(:,3), translate_angles(all_r2(:,3), all_r2(:,2)));
corrsdir(4) = corr(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,4)));
corrsdir(5) = corr(all_r2(:,2), translate_angles(all_r2(:,2), all_r2(:,4)));

%boxplot(all_r2(:,2), tuningtype)

unrot = tuningtype == 5;
rot = (tuningtype == 1 | tuningtype == 3| tuningtype == 4);
untuned = tuningtype == 4;

corrsrotated(1) = corr(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,2)));
corrsrotated(2) = corr(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,3)));
corrsrotated(3) = corr(all_r2(rot,3), translate_angles(all_r2(rot,3), all_r2(rot,2)));
corrsrotated(4) = corr(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,4)));
corrsrotated(5) = corr(all_r2(rot,2), translate_angles(all_r2(rot,2), all_r2(rot,4)));

corrsunrotated(1) = corr(all_r2(unrot,1), translate_angles(all_r2(unrot,1), all_r2(unrot,2)));
corrsunrotated(2) = corr(all_r2(unrot,1), translate_angles(all_r2(unrot,1), all_r2(unrot,3)));
corrsunrotated(3) = corr(all_r2(unrot,3), translate_angles(all_r2(unrot,3), all_r2(unrot,2)));
corrsunrotated(4) = corr(all_r2(unrot,1), translate_angles(all_r2(unrot,1), all_r2(unrot,4)));
corrsunrotated(5) = corr(all_r2(unrot,2), translate_angles(all_r2(unrot,2), all_r2(unrot,4)));

figure 
cc = ones(size(tuningtype));
cc(tuningtype == 1 | tuningtype == 3 | tuningtype == 4) = 2;
cc(tuningtype == 5) = 3;
colors = [0 0 0; 1 0 0; 0 0 1];
c = [];
for idx = 1:size(cc, 1);
	c(idx,1:3) = colors(cc(idx),:);
end
colormap(colors);
subplot(2,3,1)
scatter(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,2)), [], cc)
xlabel('\theta MC1')
ylabel('\theta BC1')
title(['corr rotated: ' num2str(corrsrotated(1)) ' corr unrotated: ' num2str(corrsunrotated(1))])
subplot(2,3,2)
scatter(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,3)), [], cc)
xlabel('\theta MC1')
ylabel('\theta MC2')
title(['corr rotated: ' num2str(corrsrotated(2)) ' corr unrotated: ' num2str(corrsunrotated(2))])
subplot(2,3,3)
scatter(all_r2(:,3), translate_angles(all_r2(:,3), all_r2(:,2)), [], cc)
ylabel('\theta BC1')
xlabel('\theta MC2')
title(['corr rotated: ' num2str(corrsrotated(3)) ' corr unrotated: ' num2str(corrsunrotated(3))])
subplot(2,3,5)
scatter(all_r2(:,1), translate_angles(all_r2(:,1), all_r2(:,4)), [], cc)
xlabel('\theta MC1')
ylabel('\theta DC')
title(['corr rotated: ' num2str(corrsrotated(4)) ' corr unrotated: ' num2str(corrsunrotated(4))])
subplot(2,3,6)
scatter(all_r2(:,2), translate_angles(all_r2(:,2), all_r2(:,4)), [], cc)
xlabel('\theta BC1')
ylabel('\theta DC')
title(['corr rotated: ' num2str(corrsrotated(5)) ' corr unrotated: ' num2str(corrsunrotated(5))])
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/tuninganglePR2R2out.eps', 'eps', [10 6])
