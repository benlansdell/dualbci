%Check all electrodes are good...
%Load nev file and plot number of all 144 channels activity

fn = './data/20130117SpankyUtah001.nev';

NEV = openNEV(fn);

nelectrodes = 144;
electrode = zeros(nelectrodes,1);

for e=NEV.Data.Spikes.Electrode
	electrode(e) = electrode(e) + 1;
end

plot(1:nelectrodes, electrode);
xlabel('Electrode #');
ylabel('# spikes');
title('')
saveplot(gcf, './worksheets/diagnostics/electrodecheck.eps');
