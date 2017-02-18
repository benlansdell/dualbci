function fileList = GetUnprocessedFiles( fn_in)
	fileList = {};
	c = 1;
	fh = fopen(fn_in);
	tline = fgetl(fh);
	while ischar(tline)
		fileList{c} = tline;
		tline = fgetl(fh);
		c = c + 1;
	end
	fclose(fh);
end