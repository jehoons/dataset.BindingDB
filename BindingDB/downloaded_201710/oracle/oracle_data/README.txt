This directory contains two files:
  README.txt               --- the procedure used by BindingDB to import a dump
  BindingDB_All_Oracle.dmp --- oracle dump by expdp

  The Oracle version is 12.1

  I used impdp to install like "impdp parfile=imp_all.txt"

--------imp_all.txt--------
Userid=BIND_2000_05/xxxxxxxxxxxx
DIRECTORY=DUMP_DIR
dumpfile=BindingDB_All_Oracle.dmp
logfile=bind_imp.log
--------imp_all.txt--------

  In your system, A DUMP_DIR should be made to hold BindingDB_All_Oracle.dmp
and bind_imp.log.
