function data = filters_sprc_relpos_target_lv(processed, nK_sp, nK_pos, nK_tar, dt_sp, dt_pos, dt_tar)
	%Prepare spike and torque data for GLM which includes spike history and cursor position (x_1, x_2) filters:
	%Spike history filter is saved as a rasied cosine function coefficients
	%
	%	y(i) ~ Pn(g(eta_i))
	%
	%where
	%
	%	eta_i = \sum y(i-j) k_sp(i) + \sum x_1(i+j) k_1(j) + \sum x_2(i+j) k_2(j)
	%
	%Usage:
	%	data = filters_sprc_relpos_lv(processed, nK_sp, nK_pos, nK_tar, dt_sp, dt_pos, dt_tar)
	%     
	%Input:
	%	processed = structure output from one of the preprocess functions.
	%	nK_sp = number of timebins used for spike history filter
	%	nK_pos = number of timebins used for relative position filters 
	%	nK_tar = number of timebins used for target position filters 
	%	dt_sp = (optional, default = binsize in processed structure) step size of spike history filter
	%		in seconds. Must be a multiple of the data's binsize.
	%	dt_pos = (optional, default = binze in processed structure) step size of relative position filter
	%		in seconds. Must be a multiple of the data's binsize
	%	dt_tar = (optional, default = binsize) step size of target position filter
	%   
	%Output:
	%	data is a structure containing the following fields:
	%		y = [nU x nB] array where y_ij is the number of spikes at time bin j for unit i.
	%		X = [nU x nB x nK] array where X_ijk is the value of covariate k, at time bin j, for unit i
	%			Note: nK = nK_sp + 2*nK_pos
	%		k = Names of each filter, a [n x 2] cell array in which each row is of the form ['filter j name', [idxj1 idxj2 ...]]
	%			Note: The second column lists indices in 1:nK to which the label applies
	%		target = trimmed target position data
	%		torque = torque data trimmed in the same way X and y are. 
	%			Note: truncated at start and end because spike and cursor trajectory are not defined for first 
	%			and last nK_sp and nK_pos timebins respectively.
	%		dtorque = trimmed dtorque
	%		ddtorque = trimmed ddtorque
	%		cursor = trimmed cursor position 
	%		dcursor = trimmed diff of cursor position
	%		ddcursor = trimmed diff of diff of cursor position	
	%  
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	nK_sp = 50; 
	%	nK_pos = 10;
	%	nK_tar = 1;
	%	dt_sp = 0.002;
	%	dt_pos = 0.05;
	%	dt_tar = 0.002;
	%	data = filters_sprc_relpos_target_lv(pre.processed, nK_sp, nK_pos, nK_tar, dt_sp, dt_pos, dt_tar);

	if (nargin < 5) dt_sp = processed.binsize; end
	if (nargin < 6) dt_pos = processed.binsize; end
	if (nargin < 7)	dt_tar = processed.binsize; end

	%Check dt's specified are valid
	assert(rem(dt_sp,processed.binsize)==0, 'Invalid dt_sp. Must be a multiple of binsize');
	assert(rem(dt_pos,processed.binsize)==0, 'Invalid dt_tar. Must be a multiple of binsize');
	assert(rem(dt_tar,processed.binsize)==0, 'Invalid dt_tar. Must be a multiple of binsize');
	steps_sp = dt_sp/processed.binsize;
	steps_pos = dt_pos/processed.binsize;
	steps_tar = dt_tar/processed.binsize;

	%Check if target vector is all zero, if it is then this recording contains no trial information
	%and the use of the 'target' model is inappropriate
	if length(unique(processed.target)) == 1
		error('BCIGLM:filters_sprc_pos_target:noTargetInformation', 'No target information during recording available during recording, use different GLM model');
	end

	nU = size(processed.binnedspikes,2);
	nB = size(processed.binnedspikes,1);
	nK = nK_sp + 2*nK_pos;

	T = nK_sp*dt_sp;
	[rcbasis, spbasis, nK_rc] = makeRCBasis(dt_sp, T);
	nK = nK_rc + 2*nK_pos + 2*nK_tar;

	data.X = zeros(nU, nB, nK);
	data.k = cell(5,3);
	data.k{1,1} = 'spike history'; 
	data.k{1,2} = 1:nK_rc;
	data.k{1,3} = dt_sp;
	data.k{2,1} = 'RU rel. pos'; 
	data.k{2,2} = (1:nK_pos) + nK_rc;
	data.k{2,3} = dt_pos;
	data.k{3,1} = 'FE rel. pos'; 
	data.k{3,2} = (1:nK_pos) + nK_rc + nK_pos;
	data.k{3,3} = dt_pos;
	data.k{4,1} = 'RU tar. pos'; 
	data.k{4,2} = (1:nK_tar) + nK_rc + 2*nK_pos;
	data.k{4,3} = dt_tar;
	data.k{5,1} = 'FE tar. pos'; 
	data.k{5,2} = (1:nK_tar) + nK_rc + 2*nK_pos + nK_tar;
	data.k{5,3} = dt_tar;
	%Record specifically which indices are spike history indices for model simulation
	data.sp_hist = data.k{1,2};

	startbin = (nK_sp*steps_sp+1);
	endlen = max(nK_pos*steps_pos, nK_tar*steps_tar);
	endbin = (nB-endlen);

	%For each unit, add data to X array
	for idx=1:nU 
		%Make stimulus vector at each timebin
		for j = startbin:endbin
			%(past) spike history
			shist = project_rc(processed.binnedspikes(j-nK_sp*steps_sp:steps_sp:j-steps_sp, idx), rcbasis);
			%(future) torque trajectory
			torqueRU = processed.cursor(j:steps_pos:(j+(nK_pos-1)*steps_pos),1);
			torqueFE = processed.cursor(j:steps_pos:(j+(nK_pos-1)*steps_pos),2);
			%(future) target position
			targetRU = processed.target(j:steps_pos:(j+(nK_pos-1)*steps_pos),1);
			targetFE = processed.target(j:steps_pos:(j+(nK_pos-1)*steps_pos),2);			
			%(future) relative cursor position to target
			reltorqueRU = torqueRU-targetRU;
			reltorqueFE = torqueFE-targetFE;			
			%Extra target-based filter
			tarRU = processed.target(j:steps_tar:(j+(nK_tar-1)*steps_tar),1);
			tarFE = processed.target(j:steps_tar:(j+(nK_tar-1)*steps_tar),2);			
			%Form stim vector
			data.X(idx,j,:) = [shist' reltorqueRU' reltorqueFE' tarRU' tarFE'];
		end
	end
	%Truncate to exclude start and end of recording where spike history 
	%and cursor trajectory aren't well defined
	data.X = data.X(:,startbin:endbin,:); %(nkt+1:end-nkt,:);
	data.y = processed.binnedspikes(startbin:endbin, :)';
	%Truncate other data for comparison, too
	data.target = processed.target(startbin:endbin,:);
	data.torque = processed.torque(startbin:endbin,:); 
	data.dtorque = processed.dtorque(startbin:endbin,:);
	data.ddtorque = processed.ddtorque(startbin:endbin,:);
	data.cursor = processed.cursor(startbin:endbin,:); 
	data.dcursor = processed.dcursor(startbin:endbin,:);
	data.ddcursor = processed.ddcursor(startbin:endbin,:);	

	data.rcbasis = rcbasis;
	data.spbasis = spbasis;
end

function sphistory_rc = project_rc(sphistory, rcbasis)
	sphistory_rc = rcbasis*sphistory;
end

function [rcbasis, spbasis, nK_rc] = makeRCBasis(dt, T)
	%Define basis of raised cosine functions
	%Hard-wired log-scale
	a = 15;
	nTotal = floor(T/dt);
	%Create basis function function
	B = @(t, a, psi, phi) iif(a*log(t-psi)>phi-pi & a*log(t-psi)<phi+pi, 0.5*(cos((a*log(t-psi)-phi))+1), ...
		true, zeros(size(t)));
	tt = 0:dt:(T-dt);
	nT = length(tt);
	%Compute phi0 and psi for desired basis
	phi1 = a*log(dt*(1-exp(-pi/a))^(-1));
	phi0 = phi1-pi;
	psi = -exp(phi0/a);
	%Compute each function for phi_i
	phi = phi0;
	%Compute matrix of basis vectors
	rc = zeros(nT, nTotal);
	for j = 1:nTotal
		for i = 1:nT
			rc(i, j) = B(tt(i), a, psi, phi);
		end
		%Normalize each column
		if norm(rc(:,j))>0
			rc(:,j) = rc(:,j)/norm(rc(:,j));
		end
		phi = phi + pi;
	end
	%Truncate to just the non-zero columns
	nK_rc = rank(rc);
	rc = rc(:,1:nK_rc);
	%Flip so time near spike is best resolved
	rc = flipud(rc);
	%Compute pseudo inverse of this matrix
	rcbasis = inv(rc'*rc)*rc';
	spbasis = rc;
end