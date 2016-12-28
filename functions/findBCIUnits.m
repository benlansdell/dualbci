function [units, durs] = findBCIUnits(data, mapping, activechannels, offset, dur)
    units = {};
    durs = [];
    %If manual just return an empty list, otherwise, find which were used based on data.pvmap
    if isempty(strfind(mapping, 'Manual'))
        %Find the options that were used during this nev file
        optidx1st = 0;
        optidx2nd = 0;
        for idx = 1:length(data.pvmap.time)
            if data.pvmap.time(idx)<offset
                optidx1st = idx;
            end
            if data.pvmap.time(idx)<offset+dur
                optidx2nd = idx;
            end
        end
        optidx1st = max(optidx1st,1);
        optidx2nd = max(optidx2nd,1);
        %If the options changed within the recording then we need to take more care...
        if optidx1st ~= optidx2nd
            for idx = optidx1st:optidx2nd
                enabledchans = data.pvmap.enabled(idx,1:4);
                units{idx-optidx1st+1} = floor(activechannels(enabledchans==1));
                %1st time
                if idx == optidx1st
                    segdur = data.pvmap.time(idx+1)-offset;
                %Last time
                elseif idx == optidx2nd
                    segdur = offset+dur-data.pvmap.time(idx);
                %Other times
                else
                    segdur = data.pvmap.time(idx+1)-data.pvmap.time(idx);
                end
                durs(idx-optidx1st+1) = segdur;
            end
        %Otherwise we don't
        else
            enabledchans = data.pvmap.enabled(optidx1st,1:4);
            units = {floor(activechannels(enabledchans==1))};
            durs = [dur];
        end
    end
end
