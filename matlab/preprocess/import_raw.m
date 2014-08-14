function t_out = import_raw(t_in, chans)
	%import_raw       Import raw extra-cellular recording data from ns3 file. Includes Blackrock electrode data, and torque data
	%
	% Usage:
	%           t_out = import_raw(t_in, chans)
	%
	% Input:
	%           t_in = a trial structure produced by import_trials().
	%			chans = (optional, default = '') if specified will only import channels requested
	%
	% Output:
	%			t_out = trial structure with channel data sampled at labview rate
	%
	% Examples:
	%           trials = import_trials('Spanky_2013-01-17-1325.mat');
	%			t = trials(117);
	%			%Read in torque channels only
	%			t = import_raw(t, 'c:138:139');

	if (nargin < 2)
		chans = '';
	end

	if length(t_in.ns3file) > 0
		if length(chans) > 0
			NS3=openNSx(t_in.ns3file, 'read', chans);
		else
			NS3 = openNSx(t_in.ns3file, 'read');
		end
	else
		t_out = t_in;
		return;
	end
	
	t_in.ns3samplerate = NS3.MetaTags.SamplingFreq;
	labviewsamplerate = 60;
	flank = t_in.flank;
	ns3data_raw = double(NS3.Data);
	ns3data = resample(double(NS3.Data'), labviewsamplerate, t_in.ns3samplerate)';
	%Only find data within trial span
	s = size(ns3data);
	nsamples = s(2);
	s = size(ns3data_raw);
	nsamples_raw = s(2);
	times = t_in.offset + (0:(nsamples-1))/labviewsamplerate;
	times_raw = t_in.offset + (0:(nsamples_raw-1))/t_in.ns3samplerate;
	withintrial = (times > t_in.starttime) & (times < t_in.endtime);
	withintrial_raw = (times_raw > t_in.starttime) & (times_raw < t_in.endtime);
	withintrial_flank = (times > t_in.starttime-flank) & (times < t_in.endtime+flank);
	withintrial_raw_flank = (times_raw > t_in.starttime-flank) & (times_raw < t_in.endtime+flank);
	t_in.ns3data = ns3data(:,withintrial);
	t_in.ns3data_raw = ns3data_raw(:,withintrial_raw);
	t_in.ns3data_flank = ns3data(:,withintrial_flank);
	t_in.ns3data_raw_flank = ns3data_raw(:,withintrial_raw_flank);
	t_out = t_in;
end
