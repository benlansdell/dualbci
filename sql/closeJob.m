function closejob(conn, id)
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	status = 1;
	tablename = 'Jobs';
	colnames = {'`end date`', 'status'};
	sqldata = { stamp, status};
	whereclause = ['WHERE id = ' num2str(id)];
	update(conn,tablename,colnames,sqldata,whereclause);
end