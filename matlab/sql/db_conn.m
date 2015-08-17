function conn = db_conn(c)
	if nargin < 1
		conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
			'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
		return
	end		
	a = exec(c, 'SELECT * FROM recordings')
	if strcmp(a.Message, 'Invalid connection.')
		conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
			'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
	else
		conn = c;
	end
end