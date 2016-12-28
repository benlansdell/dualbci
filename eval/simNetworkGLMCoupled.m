function tsp = simNetworkGLMCoupled(nU)
	% simNetworkGLMCoupled.m
	%
	% Test code for simulating and fitting a coupled GLM (nU = number of neurons,
	% default 30)
	%
	% Returns spike times
	if (nargin < 1)	nU = 30; end

	global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
	RefreshRate = 100;
		
	%Randomly generate coupling matrix
	%Excitation terms 10%
	%Inhibition terms 20%
	pE = 0.1;
	pI = 0.2;
	excitation = rand(nU);
	inhibition = rand(nU);
	
	coupling = zeros(nU);
	coupling(excitation < pE) = 1;
	coupling(inhibition < pI) = -1;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% 1 Set parameters and display for GLM  %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	ggsims = {};
	DTsim = .01; % Bin size for simulating model & computing likelihood (in units of stimulus frames)
	nkt = 20;    % Number of time bins in filter;
	ttk = [-nkt+1:0]'; 
	for idx = 1:nU
		ggsims{idx} = makeSimStruct_GLM(nkt,DTsim);  % Create GLM struct with default params
		k = ggsims{idx}.k;  % Stimulus filter with random offsets
		offset = floor(-4*rand(1));
		filtersize = 2*rand(1)-1;
		ggsims{idx}.k = shift(k,offset)*filtersize;
	end
	ggsim = makeSimStruct_GLMcpl(ggsims{:});
	
	% Make some coupling kernels
	for i = 1:nU
		for j = 1:nU
			[iht,ihbas,ihbasis] = makeBasis_PostSpike(ggsim.ihbasprs,DTsim);
			%Excitatory connections
			if coupling(i,j) == 1
				hhcpl = ihbasis*[.6;.47;.25;0;0]*2*rand(1);
			%Inhibitory connections
			elseif coupling(i,j) == -1
				hhcpl = ihbasis*[-1;-1;0;0;.25]*2*rand(1);
			%No connection
			else
				hhcpl = zeros(size(iht));
			end
			ggsim.ih(:,i,j) = hhcpl; % 2nd cell coupling to first
		end
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% 2 Make GWN stimulus & simulate the glm model response %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	slen = 50; % Stimulus length (frames) & width (# pixels)
	swid = size(ggsim.k,2);
	Stim = 2*randn(slen,swid);  % Gaussian white noise stimulus
	[tsp, vmem,Ispk] = simGLM(ggsim, Stim); % Simulate GLM response
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%% 3 Generate some training data %%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	slen = 2500;  % Stimulus length (frames);  More samples gives better fit
	Stim = round(rand(slen,swid))*4-2;  %  Run model on long, binary stimulus
	[tsp,vmem,ispk] = simGLM(ggsim,Stim);  % run model
end