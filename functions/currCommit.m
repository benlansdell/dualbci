function comm = currCommit()
	%Try 10 times to get a good answer, otherwise return NULL
	nTries = 1;
	comm = 'asdfasdfasdfasdfasdfasdfasdf';
	while length(comm) > 7 & nTries <= 10
		comm = git('rev-parse --short HEAD');
		nTries = nTries + 1;
	end
	if length(comm) > 7
		comm = 'NULL';
	end
end