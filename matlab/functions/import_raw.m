function t_out = import_raw(t_in)
        %import_raw       Import raw extra-cellular recording data from ns3 file. Includes Blackrock electrode data, and torque data
        %
        % Usage:
        %                       t_out = import_raw(t_in)
        %
        % Input:
        %                       t_in = a trial structure produced by import_trials().
	%
	% Output:
	%			t_out = trial structure with channel data sampled at labview rate
        %
        % Examples:
        %                       trials = import_trials('Spanky_2013-01-17-1325.mat');
	%			t = trials(117);
	%			t = import_raw(t);
        if length(t_in.ns3file) > 0
                NS3=openNSx(t_in.ns3file, 'read');
        else
                return;
        end

        t_in.ns3samplerate = NS3.MetaTags.SamplingFreq;
        labviewsamplerate = 60;

	t_in.ns3data_raw = double(NS3.Data);
	t_in.ns3data = resample(double(NS3.Data'), labviewsamplerate, t_in.ns3samplerate)';
	%Only find data within trial span
	s = size(t_in.ns3data);
	nsamples = s(2);
	s = size(t_in.ns3data_raw);
	nsamples_raw = s(2);
	times = t_in.offset + (0:(nsamples-1))/labviewsamplerate;
	times_raw = t_in.offset + (0:(nsamples_raw-1))/t_in.ns3samplerate;
	withintrial = (times > t_in.starttime) & (times < t_in.endtime);
	withintrial_raw = (times_raw > t_in.starttime) & (times_raw < t_in.endtime);
	t_in.ns3data = t_in.ns3data(:,withintrial);
	t_in.ns3data_raw = t_in.ns3data_raw(:,withintrial_raw);
        t_out = t_in;
end
