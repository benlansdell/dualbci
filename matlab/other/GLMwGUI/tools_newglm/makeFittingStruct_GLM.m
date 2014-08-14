function gg = makeFittingStruct_GLM(sta,DTsim,timeafterspk,ppms)
% gg = makeFittingStruct_GLM(sta,DTsim,glmstruct,cellnumToFit);
%
% Initialize parameter structure for fitting of GLM model,
% normal parametrization of stim kernel
%
% Inputs:  sta = initial guess at kernel (or use all zeros if unknown)

% Set up structure
gg.k = [];
gg.dc = 0;
gg.ih = [];
gg.iht = [];
gg.ihbas = [];
gg.ihbas2 = [];
gg.ihbasprs = [];
gg.kt = [];
gg.ktbas = [];
gg.kbasprs = [];
gg.tsp = [];
gg.tspi = [];
gg.dt = DTsim;
gg.ih2 = [];
gg.ihbas2 = [];
gg.ihbasprs2 = [];
gg.tsp2 = [];
gg.couplednums = [];

% === Make temporal basis for stimulus filter =======================
[nkt,nkx] = size(sta);
% % ----- Set up temporal basis for stimulus kernel -----------
kbasprs.neye = 5; % Number of "identity" basis vectors near time of spike;
kbasprs.ncos = 5; % Number of raised-cosine vectors to use  
kbasprs.kpeaks = [0 round(nkt/3)];  % Position of first and last bump (relative to identity bumps)
kbasprs.b = 3; % Offset for nonlinear scaling (larger -> more linear)
ktbas = makeBasis_StimKernel(kbasprs,nkt);
gg.ktbas = ktbas;
gg.kbasprs = kbasprs;

% ======================================================================
% Set up basis for post-spike kernel

global flagval Basispars
if isempty(flagval)
    Basispars = [8,0.9,4];
end
ihbasprs.ncols = Basispars(1);  % Number of basis vectors for post-spike kernel
% SangWook's change 8/13
ihbasprs.hpeaks = [DTsim*ppms*Basispars(3) DTsim*ppms*timeafterspk*Basispars(2)];  % Peak location for first and last vectors

ihbasprs.b = .4;  % How nonlinear to make spacings
ihbasprs.absref = DTsim*ppms; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim,timeafterspk,ppms);
gg.iht = iht;
gg.ihbas = ihbas;
gg.ihbasprs = ihbasprs;
gg.ih = zeros(size(ihbas,2),1);

% % ==================================================================
% set up initial K params
gg.kt = inv(gg.ktbas'*gg.ktbas)*gg.ktbas'*sta;
gg.k = gg.ktbas*gg.kt;