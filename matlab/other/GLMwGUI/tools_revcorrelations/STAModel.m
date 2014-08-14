function [spike_time,logli] = STAModel(Stim,Neuron,FR)

global sta_scale_factor ppms

if isempty(sta_scale_factor)
    ii = 2;
    while ii > 1
        
        len = length(Stim);
        pred_spt_len = (len/(ppms*1000))*FR;
        
        if ii == 2 && isempty(sta_scale_factor)
            sta_scale_factor = 1;
        end
        
        [spike_time,logli] = subSTAModel(Stim,Neuron,FR,sta_scale_factor);
        
        if isempty(spike_time)
            sta_scale_factor = sta_scale_factor.*2;
        elseif ~isempty(spike_time) && ii == 2
            esc = 1;
            while esc
                sta_scale_factor = (sta_scale_factor*(pred_spt_len/length(spike_time)));
                [spike_time,logli] = subSTAModel(Stim,Neuron,FR,sta_scale_factor);
                if length(spike_time) >= (pred_spt_len - round(pred_spt_len*.1)) && ...
                        length(spike_time) <= (pred_spt_len + round(pred_spt_len*.1))
                    esc = 0;
                end
            end
            ii = ii - 1;
        end
    end
else
    [spike_time,logli] = subSTAModel(Stim,Neuron,FR,sta_scale_factor);
end

end

function [spike_time,logli] = subSTAModel(Stim,Neuron,FR,sta_scale_factor)

spike_time = [];

Stim = Stim(:);

STAfilt = Neuron.STA;
STAthreshold = sta_scale_factor.*(Neuron.Thresh)./FR;

filtresp = sameconv(Stim,STAfilt);
filtresp = filtresp - min(filtresp);
filtresp = 1 + round(filtresp.*100); % scaling filtered response

STAthreshold = interp1(1:length(STAthreshold),STAthreshold,linspace(1,...
    length(STAthreshold),max(filtresp)));

lambda = zeros(1,length(filtresp));
t = 1;
potential = 0;
isiNext = exprnd(1);
while t <= length(Stim)-1
    %integrate lambda
    lambda(t) = STAthreshold(filtresp(t));
    potential = potential + lambda(t);

    %if driving a spike
    if potential > isiNext
        %save spike
        spike_time = [spike_time,t];
        %draw next firing potential
        isiNext = isiNext + exprnd(1);
    end

    %continue
    t = t+1;
end

rr = STAthreshold(filtresp);
etol = 1e-100;
iiz = find(rr<=etol);
rr(iiz) = etol;
Trm1 = sum(rr)*0.1;
Trm2 = -sum(log(rr(spike_time)));
logli = Trm1 + Trm2;
end