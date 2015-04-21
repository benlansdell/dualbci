%Find all units common between start and end recording
mat_sprc_pos = load('./expts/2dmanualpos_sprc_pos_lv_def.mat');
mat_sprc = load('./expts/2dmanualpos_sprc_def.mat');
expts_sprc_pos = mat_sprc_pos.expts;
expts_sprc = mat_sprc.expts;
nE = length(expts_sprc);
nC = length(conds);

DP = zeros(nU, nE, 4+nC);

%Note their change in deviance between sprc and sprc_pos_lv
for idx = 1:nE
%Do analysis per unit:
%For each unit, find all experiments containing that unit

%Append this experiment's info (deviances for each case, and the seconds of recording in between)

end


%Plot:
%Scatter plot showing DeltaD_p before (xaxis) and after (yaxis) for all expts

%Color this plot in different ways to highlight (hopefully) interesting features

%Color by minutes under manual recording, minutes under brain recording, minutes used to drive BCI

%Fit LM:

%Make the same plots but with all units' information aggregated