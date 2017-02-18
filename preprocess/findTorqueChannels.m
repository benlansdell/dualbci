function chans = findTorqueChannels(nevfile)
	nsxfile = [nevfile(1:end-3), 'ns3'];
	nsx = openNSx(nsxfile, 'read', 't:0:1');
	nChannels = size(nsx.Data, 1);
	fe = nChannels-2;
	ru = nChannels-1;
	chans = [fe, ru];
end