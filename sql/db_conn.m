function conn = db_conn(c)
	if nargin < 1
		conn = database('',databaseuser,databasepwd,'com.mysql.jdbc.Driver', ...
			databaseurl);
		return
	end		
	a = exec(c, 'SELECT * FROM recordings');
	if strcmp(a.Message, 'Invalid connection.')
		conn = database('',databaseuser,databasepwd,'com.mysql.jdbc.Driver', ...
			databaseurl);
	else
		conn = c;
	end
end