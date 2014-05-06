function [S, f] = plot_lfp_ns3(ns3s, chans, fn)
        %plot_lfp_ns3       Plot spectral content of nsx file for specified channels. Output to eps file
        %
        % Usage:
        %                       plot_lfp_ns3(ns3s,chans,fn)
        %
        % Input:
        %                       ns3s = cell array of input ns3 files
        %                       chans = [chmin chmax] specifies the range of channels to read from
        %                       fn = output file name for plot
        %
        % Examples:
        %                       fn = './worksheets/diagnostics/plots/test_lfp_ns3.eps';
        %                       chans = [138, 139];
        %                       ns3s = {'./testdata/20130117SpankyUtah005.ns3'};
	%			plot_lfp_ns3(ns3s, chans, fn);

        close all;
	fig = figure;

        if (nargin < 3)
                throw(MException('Argin:MoreExpected', 'More input arguments expected'));
        end
        nsxlfps = [];
        totaltime = 0;

        %Load data from each file
        for idx = 1:length(ns3s)
                ns3file = ns3s{idx};     
                %Load torque data from NS3 file
                NS3 = openNSx(ns3file, 'read', ['c:' num2str(chans(1)) ':' num2str(chans(2))]);
                nsxlfp = double(NS3.Data)';
                size(nsxlfp);
                nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
                duration = double(NS3.MetaTags.DataDurationSec);
                %Concatenate to previously loaded files
                totaltime = totaltime + duration;
                nsxlfps = [nsxlfps; nsxlfp];
        end

        %Plot the data
        params.Fs = nsxsamplerate;
        [S, f] =  mtspectrumc(nsxlfps, params);
        for idx = 1:length(chans)
               plot_vector(S(:,idx),f);
               title('Frequency content of channel ' num2str(chans(idx)))
	       %Write file
	       saveplot(gcf, [fn '_ch_' num2str(chans(idx)) '.eps']);
        end
end
