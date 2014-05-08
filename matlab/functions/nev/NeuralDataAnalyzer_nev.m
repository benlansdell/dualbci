function NeuralDataAnalyzer_nev(nevfile, fn_out, unit, torque_axis, time_range)
	%NeuralDataAnalyzer_nev Function to prepare data for SangWook's NeuralDataAnalyzer program. Reads in nev file and extracts spike times within a given range, 
	%		for a specified unit, and a specified torque axis.	
	%
	%		Usage:
	%			NeuralDataAnalyzer_nev(nevfile, fn_out, unit, torque_axis, timerange)
	%
	% 		Input:
	%			nevfile = file to process
	%			fn_out = file to output diagnostic plots to, if desired
	%			unit = unit name to process. In the format 'Electrode.SortCode'
	%			torque_axis = 1 or 2, 1 = FE, 2 = RU
	%			time_range = [tmin, tmax] range of times from which to extract spikes
	%		
	%		Output:
	%			none
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
	%			fn_out = './worksheets/NeuralDataAnalyzer/20130117SpankyUtah001_0_100_16.2';
	%			unit = '16.2';
	%			torque_axis = 1;
	%			time_range = [0 100];
	%			NeuralDataAnalyzer_nev(nevfile, fn_out, unit, torque_axis, time_range);
	
	%%%%%%%%%%%%%%%%%%%%%%
	%Extract spiking data%
	%%%%%%%%%%%%%%%%%%%%%%
	nunits = 5;
	nE = 128;
	nU = nunits*nE;
	ns3file = [nevfile(1:end-3) 'ns3'];
	NEV = openNEV(nevfile);
	%Find the duration and sample rate of the nev file recording
	nevsamplerate = NEV.MetaTags.TimeRes;
	%Convert spike times into doubles
	spiketimes = double(NEV.Data.Spikes.TimeStamp);
    elecs = cell(1,nU);
    spikemuas = struct('times', elecs);
    unitnames = cell(1,nU);
    averate = zeros(1,nU);
    for idx=1:nU
        spikemuas(idx).times = [0];    
    end
    for i=1:length(spiketimes)
       	E = NEV.Data.Spikes.Electrode(i);
       	u = NEV.Data.Spikes.Unit(i);
       	U = single((E-1)*nunits)+single(u)+1;
       	sptime = spiketimes(i);
       	if (sptime/nevsamplerate) > time_range(1) & (sptime/nevsamplerate) < time_range(2)
	       	spikemuas(U).times = [spikemuas(U).times; sptime];
    	end
       	unitnames{U} = [num2str(E) '.' num2str(u)];
    end
    %Find unit we want, otherwise throw error
    nu = find(strcmp(unit, unitnames));
    if isempty(nu)
    	error('Unit not found in nev file');
    end
    %Save
    spikes = spikemuas(nu).times;
    save([fn_out '_sptimes.mat'], 'spikes');

	%%%%%%%%%%%%%%%%%%%%%
	%Process torque data%
	%%%%%%%%%%%%%%%%%%%%%
	clear torque;
	NS3 = openNSx(ns3file, 'read', 'c:138:139');
	nsxtorque = double(NS3.Data);
	nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
	%Switch sign of FE axis for coordinate consistency
	nsxtorque(2,:)=-nsxtorque(2,:);
    for j=1:2
        %Scale from uint16 value to proportion
        nsxtorque(j,:) = nsxtorque(j,:)/(2^15);
        %Subtract mean
        nsxtorque(j,:) = nsxtorque(j,:)-mean(nsxtorque(j,:));
    end
    %Extract axis we want
    torque = nsxtorque(torque_axis,:);
    %Only look at times in specified range
    tt = (1:length(torque))/nsxsamplerate;
    withinrange = (tt > time_range(1)) & (tt < time_range(2));
    torque = torque(withinrange);
    %Resample at rate of spiking data
    torque=resample(torque,nevsamplerate,nsxsamplerate);
    save([fn_out '_torque.mat'], 'torque');