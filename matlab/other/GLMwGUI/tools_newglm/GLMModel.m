function [spike_time,logli] = GLMModel(Stim,Neuron,FR)

global glm_scale_factor ppms

if isempty(glm_scale_factor)
    ii = 2;
    while ii > 1
        len = length(Stim);
        pred_spt_len = (len/(ppms*1000))*FR;
        
        if ii == 2 && isempty(glm_scale_factor)
            glm_scale_factor = 1;
        end
        
        [spike_time,logli] = subGLMModel(Stim,Neuron,FR,glm_scale_factor);
        
        if isempty(spike_time)
            glm_scale_factor = glm_scale_factor.*2;
        elseif ~isempty(spike_time) && ii == 2
            esc = 1;
            while esc
                glm_scale_factor = glm_scale_factor*(pred_spt_len/length(spike_time));
                [spike_time,logli] = subGLMModel(Stim,Neuron,FR,glm_scale_factor);
                if length(spike_time) >= (pred_spt_len - round(pred_spt_len*.15)) && ...
                        length(spike_time) <= (pred_spt_len + round(pred_spt_len*.3))
                    esc = 0;
                end
            end
            ii = ii - 1;
        end
    end
else
    [spike_time,logli] = subGLMModel(Stim,Neuron,FR,glm_scale_factor);
end

end

function [spike_time,logli] = subGLMModel(Stim,Neuron,FR,glm_scale_factor)

etol = 1e-100;
Stim = Stim(:)';

Stimfilt = Neuron.Stimfilt(:)';
Histfilt = Neuron.Histfilt(:)';

Histfilt = interp1(1:length(Histfilt),Histfilt,linspace(1,length(Histfilt),1000));
neuthresh = @(x) exp(x);

inputCurrent = conv(Stim, Stimfilt, 'full');
historyCurrent = zeros(1,length(inputCurrent));
lambda = zeros(1,length(inputCurrent));
spike = false(size(Stim));

% generate spikes, for all timebins
t = 1;
potential = 0;
isiNext = exprnd(1);
while t <= length(Stim)-1
    %integrate lambda
    lambda(t) = 1e-10.*glm_scale_factor*FR.*(neuthresh(inputCurrent(t) + historyCurrent(t)));
    potential = potential + lambda(t);

    %if driving a spike
    if potential > isiNext
        %save spike
        spike(t) = true;

        %draw next firing potential
        isiNext = isiNext + exprnd(1);

        %add history filter to history current
        future = min(t+length(Histfilt)-1,length(inputCurrent))-t;
        historyCurrent(t+1:t+future) = historyCurrent(t+1:t+future) +...
            Histfilt(1:future);
    end

    %continue
    t = t+1;
end

glmspk = spike;
spike_time = [];
for i = 1:length(glmspk)
    if glmspk(i)
        spike_time = [spike_time,i];
    end
end

iiz = find(lambda <= etol);
lambda(iiz) = etol;
        
Trm1 = sum(lambda)*.1;
Trm2 = -sum(log(lambda(spike_time)));
logli = (Trm1+Trm2)*0.01;
end