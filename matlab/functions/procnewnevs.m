function procnewnevs 
	global matpath
	global nevpath
	global metaData
	matpath = '../unprocessedmat'; 
	nevpath = './blackrock';
	load([matpath '/Spanky_LogMetaData2.mat']);
	anp = dir([matpath '/*2014*.mat']);
	
	for idx = 1:length(anp)
		fn = anp(idx).name
		%system(['cp ' matpath '/' fn ' ~/projects/bci/unprocessedmat']);
		load([matpath '/' fn]);
		NEV_FindTimesForMAT1b(data)
	end
end