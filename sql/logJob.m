function id = logJob(conn, scriptname, scriptdesc)
	host = hostname();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	comm = currCommit();
	tablename = 'Jobs';
	status = 0;
	fitcols = {'name', 'description', 'computer', '`start date`', 'commit', 'status'};
	sqldata = { scriptname, scriptdesc, host, stamp, comm, status};
	datainsert(conn,tablename,fitcols,sqldata);
	id = fetch(exec(conn, 'SELECT LAST_INSERT_ID()'));
	id = id.Data{1};
end