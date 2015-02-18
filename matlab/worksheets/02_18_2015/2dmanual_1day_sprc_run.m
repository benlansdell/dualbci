matfile_in = './expts/2dmanualpos_1day.mat';
settings = setupExperiment('sprc_def');
exptname = '2dmanualpos_sprc_def';
expts = runExperiment(matfile_in, settings, exptname);
save(['./expts/' exptname '.mat'], expts);