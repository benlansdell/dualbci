function settings = setupExperiment(expttype)
	%Fit GLM on start and end recordings of all experiments located.
    %       See source code for details on each model.
    %
    %Usage:
    %   settings = setupExperiment(expttype)
    %
    %Input:
    %   expttype = one of:
    %       'sprc_def'
    %
    %Output:
    %	settings = structure array listing experiment settings:
    %		.binsize = bin size
    %       .duration = duration to truncate all recordings to
    %       .threshold = minimum firing rate below which to ignore unit
    %       .offset = temporal offset to apply between torque and spiking signals
    %       .const = whether to fit a constant term to models
    %       .modelname = expttype
    %       .process = function handle to appropriate preprocess function
    %       .filters = function handle to appropriate filter function for model
    %       .matfiledir = location of mat files
    %       .nevfiledir = location of nev files
    %
    %Test code:
    %   settings = setupExperiment('sprc_def');

    switch expttype
        %Model with only spiking history, using raised cosine basis
        case 'sprc_def'
            %Preprocess settings
            settings.binsize = 0.002;
            settings.duration = 360; 
            settings.threshold = 3;
            settings.offset = 0;
            settings.const = 'on';
            settings.process = @(nev, mat) preprocess_spline(nev, settings.binsize, settings.threshold, settings.offset);
            %Filter settings
            settings.dt_sp = settings.binsize;
            settings.nK_sp = 100;
            settings.nK_pos = 0;
            settings.filters = @(proc) filters_sprc_pos(proc, settings.nK_sp, settings.nK_pos);
        otherwise
            error('See help setupExperiment for valid settings names');
    end
    settings.modelname = expttype;
    settings.matfiledir = '/home/lansdell/fairbanks/labview/Moritz Primate MAT Files 2012-2013/';
    settings.nevfiledir = '/home/lansdell/fairbanks/blackrock/';
end