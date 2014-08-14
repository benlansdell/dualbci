function [spike_time,logli] = STCModel(Stim,Neuron,FR)

global stc_scale_factor ppms

if isempty(stc_scale_factor)
    ii = 2;
    while ii > 1

        len = length(Stim);
        pred_spt_len = (len/(ppms*1000))*FR;

        if ii == 2 && isempty(stc_scale_factor)
            stc_scale_factor = 1;
        end

        [spike_time,logli] = subSTCModel(Stim,Neuron,FR,stc_scale_factor);

        if isempty(spike_time)
            stc_scale_factor = stc_scale_factor.*2;
        elseif ~isempty(spike_time) && ii == 2
            esc = 1;
            while esc
                stc_scale_factor = (stc_scale_factor*(pred_spt_len/length(spike_time))); % divided by 10 for mean divided by 1 for std
                [spike_time,logli] = subSTCModel(Stim,Neuron,FR,stc_scale_factor);
                if length(spike_time) >= (pred_spt_len - round(pred_spt_len*.1)) && ...
                        length(spike_time) <= (pred_spt_len + round(pred_spt_len*.1))
                    esc = 0;
                end
            end
            ii = ii - 1;
        end
    end
else
    [spike_time,logli] = subSTCModel(Stim,Neuron,FR,stc_scale_factor);
end

end

function [spike_time,logli] = subSTCModel(Stim,Neuron,FR,stc_scale_factor)

etol = 1e-100;

spike_time = [];

Stim = Stim(:);

Mode1filt = Neuron.Mode1;
Mode2filt = Neuron.Mode2;

Mode1filtresp = sameconv(Stim,Mode1filt);
Mode2filtresp = sameconv(Stim,Mode2filt);

Mode1filtresp = Mode1filtresp - min(Mode1filtresp);
Mode2filtresp = Mode2filtresp - min(Mode2filtresp);

Mode1filtresp = 1 + round(Mode1filtresp.*100);
Mode2filtresp = 1 + round(Mode2filtresp.*100);

Mode1threshold = (Neuron.Mode1Thresh);
Mode2threshold = (Neuron.Mode2Thresh);

Mode1threshold = interp1(1:length(Mode1threshold),Mode1threshold,...
    linspace(1,length(Mode1threshold),max(Mode1filtresp)));
Mode2threshold = interp1(1:length(Mode2threshold),Mode2threshold,...
    linspace(1,length(Mode2threshold),max(Mode2filtresp)));

STCthreshold = stc_scale_factor.*(Mode1threshold'*(Mode2threshold))./FR;

lambda = zeros(1,length(Mode1filtresp));
t = 1;
potential = 0;
isiNext = exprnd(1);
Trm1 = 0;
Trm2 = 0;
rrsum = 0;
logli = 0;
while t <= length(Stim)-1
    %integrate lambda
    lambda(t) = STCthreshold(Mode1filtresp(t),Mode2filtresp(t));
    rr = lambda(t);
    if rr < etol
        rr = etol;
    end
    potential = potential + lambda(t);

    % if driving a spike
    if potential > isiNext
        %save spike
        spike_time = [spike_time,t];
        %draw next firing potential
        isiNext = isiNext + exprnd(1);
        Trm1 = rr*5;
        Trm2 = -log(rr)/2.5;
        rrsum = Trm1 + Trm2;
    else
        Trm1 = rr*5;
        rrsum = Trm1;
    end

    logli = logli + rrsum;
    %continue
    t = t+1;
end
end