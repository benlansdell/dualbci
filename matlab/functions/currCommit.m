function comm = currCommit()
	comm = git('rev-parse --short HEAD');
	%Filler to stop git() function from including output of next line's 
	%execution (who knows why???)
	a = 0;
	b = 0;
	c = 0;
end