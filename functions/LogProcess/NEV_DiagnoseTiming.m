% This pile of scripts is used to hand-fix timing mismatch errors between
% NEV and labview recordings

% This file assumes a data object exists that we want to fix

%% Overiew of timing information for current data

SR = data.sampleRate;
fprintf('\nMatch Summary:\n')
fprintf('NEV\tCH #1\tCH#2\tCH#3\tCH#4\tMM:SS Toffset\n')
Toffsets = cell2mat({data.nev.Toffset}');
num2str([(1:length(data.nev))', ...
        cell2mat({data.nev.chans}'), ...
        floor(Toffsets(:,1)/SR/60), ...
        mod(Toffsets(:,1)/SR,60)],...
        '%d\t%2.1f\t%2.1f\t%2.1f\t%2.1f\t%d:%2.3f\n')

    ok = {'---FAIL---','ok'};
fprintf('Segment Info (**could be same segment matched to multiple NEVs**):\n')
fprintf('NEV\tTime (s)\tMap\t\tMatch ok?\n')
for i=1:length(data.nev)
    fprintf('%d\t',i);
    if ~isempty(data.nev(i).segInfo)
        seg = data.nev(i).segInfo;
        durationSec = seg.tpause-seg.trun;
        fprintf('%d:%2.f\t%s\t%s\n',floor(durationSec/60),mod(round(durationSec),60),seg.map,ok{data.nev(i).ok+1})
    else
        fprintf('\n');
    end
end

%% Brute force parallel search of timing offset with guessed channels

