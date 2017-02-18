function [processed, wghts] = resampleProcessed(processed_orig, nsamp, reftorqdist)
	%Brain control, etc, data
	processed = processed_orig;
	nTar = size(processed_orig.torque,1);
	%Generate nsamp samples from reftorqdist (Manual control)
	load(reftorqdist);
	nBins = 30;
	binnedreftorq = zeros(nBins, nBins);
	maxt = max(max(torque), max(processed_orig.torque));
	mint = min(min(torque), min(processed_orig.torque));
	maxtFE = maxt(1);
	mintFE = mint(1);
	maxtRU = maxt(2);
	mintRU = mint(2);

	for idx = 1:size(torque,1)
		%compute FE, RU bin
		idxFE = max(min(ceil(nBins*(torque(idx,1)-mintFE)/(maxtFE - mintFE)), nBins), 1);
		idxRU = max(min(ceil(nBins*(torque(idx,2)-mintRU)/(maxtRU - mintRU)), nBins), 1);
		binnedreftorq(idxFE, idxRU) = binnedreftorq(idxFE, idxRU) + 1;
	end

	binnedtargtorq = zeros(nBins, nBins);
	%maxt = max(processed_orig.torque);
	%mint = min(processed_orig.torque);
	%maxtFE = maxt(1);
	%mintFE = mint(1);
	%maxtRU = maxt(2);
	%mintRU = mint(2);

	for idx = 1:nTar
		%compute FE, RU bin
		idxFE = max(min(ceil(nBins*(processed_orig.torque(idx,1)-mintFE)/(maxtFE - mintFE)), nBins), 1);
		idxRU = max(min(ceil(nBins*(processed_orig.torque(idx,2)-mintRU)/(maxtRU - mintRU)), nBins), 1);
		binnedtargtorq(idxFE, idxRU) = binnedtargtorq(idxFE, idxRU) + 1;
	end

	G = fspecial('gaussian', [6 6], 1);
	filtertargettorq = imfilter(binnedtargtorq/sum(sum(binnedtargtorq)), G, 'same');
	filterreftorq = imfilter(binnedreftorq/sum(sum(binnedreftorq)), G, 'same');

	wghts = zeros(nTar,1);
	for idx = 1:nTar
		idxFE = max(min(ceil(nBins*(processed_orig.torque(idx,1)-mintFE)/(maxtFE - mintFE)), nBins), 1);
		idxRU = max(min(ceil(nBins*(processed_orig.torque(idx,2)-mintRU)/(maxtRU - mintRU)), nBins), 1);
		wghts(idx) = filterreftorq(idxFE,idxRU)/filtertargettorq(idxFE, idxRU);
	end

	%Now take random sample from
	rsidx = randsample(nTar,nsamp,true,wghts);
	processed.binnedspikes = processed.binnedspikes(rsidx,:);
	processed.rates = processed.rates(rsidx,:);
	processed.torque = processed.torque(rsidx,:);
	processed.dtorque = processed.dtorque(rsidx,:);
	processed.ddtorque = processed.ddtorque(rsidx,:);

end