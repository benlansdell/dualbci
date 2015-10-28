matfile_in = './expts/2dmanualpos_1day.mat';
settings = setupExperiment('sprc_pos_def');

%For all expts
load(matfile_in);

for expt = expts
	data = load([settings.matfiledir expt.matfile_start]);
	display(['Mat file: ' expt.matfile_start]);
	display(['Labview sample rate: ' num2str(data.data.sampleRate)])
	data = load([settings.matfiledir expt.matfile_end]);
	display(['Mat file: ' expt.matfile_end]);
	display(['Labview sample rate: ' num2str(data.data.sampleRate)])
end

%Matfiles whose sample rate is NOT 60Hz:
%Spanky_2013-02-27-1259.mat: 59Hz
%Spanky_2013-03-19-1408.mat: 48Hz
%Spanky_2013-03-20-1330.mat: 48Hz
%Spanky_2013-03-21-1308.mat: 48Hz

%Why???