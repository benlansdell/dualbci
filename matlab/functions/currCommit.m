function comm = currCommit()
	comm = git('rev-parse --short HEAD');
end